"""Command-line orchestration of preflight and PyExoCross validation."""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
import sys
from typing import Any

import numpy as np
import pandas as pd
import yaml

from .calculations import (
    calculate_lifetimes,
    compare_partition_functions,
    reference_partition_function,
)
from .exomol_validator import (
    first_excited_energy,
    validate_pf,
    validate_states,
    validate_transitions,
)
from .metadata import parse_ion
from .pyexocross_adapter import (
    PyExoCrossAdapter,
    PyExoCrossExecutionError,
    PyExoCrossUnavailable,
    pyexocross_info,
    stage_local_files,
)
from .reporting import normalize_stick_output, save_partition_plots, write_report

LOGGER = logging.getLogger("validation")


def _load_config(path: Path | None) -> dict[str, Any]:
    if path is None:
        return {}
    with path.open(encoding="utf-8") as handle:
        value = yaml.safe_load(handle)
    if not isinstance(value, dict):
        raise ValueError("configuration root must be a mapping")
    return value


def _write_issues(path: Path, issues: list[dict[str, object]]) -> None:
    if issues:
        pd.DataFrame(issues).to_csv(path, index=False)
    else:
        pd.DataFrame(columns=["severity", "message"]).to_csv(path, index=False)


def _temperatures(
    pf: pd.DataFrame, requested: list[float] | None
) -> np.ndarray:
    if requested:
        result = np.asarray(sorted(set(requested)), dtype=float)
    else:
        available = pf["temperature"].to_numpy(float)
        # Data-aware low/mid/high smoke temperatures, selected from the actual
        # file range and rounded to the nearest existing integer grid point.
        fractions = np.array([0.10, 0.50, 0.90])
        indices = np.unique(np.round(fractions * (len(available) - 1)).astype(int))
        result = np.unique(np.round(available[indices]))
    if len(result) < 3:
        raise ValueError("at least three distinct temperatures are required")
    if result.min() < pf["temperature"].min() or result.max() > pf["temperature"].max():
        raise ValueError("requested temperatures must lie inside the .pf range")
    return result


def _status(results: list[dict[str, object]], runtime_errors: list[str]) -> str:
    if runtime_errors or any(not bool(result["passed"]) for result in results):
        return "FAIL"
    if any(result["warnings"] for result in results):
        return "PASS_WITH_WARNINGS"
    return "PASS"


def _to_wavenumber_range(
    values: tuple[float, float], unit: str
) -> tuple[float, float]:
    """Convert an ordered wavenumber/vacuum-wavelength interval to cm-1."""
    low, high = map(float, values)
    if low <= 0 or high <= low:
        raise ValueError("--range requires positive MIN < MAX")
    if unit == "cm-1":
        return low, high
    factor = 1e7 if unit == "nm" else 1e4
    return factor / high, factor / low


def validate_one(
    *,
    ion: str,
    states_path: Path,
    trans_path: Path,
    pf_path: Path,
    output_dir: Path,
    charge: int | None,
    temperatures: list[float] | None,
    wn_range: tuple[float, float],
    skip_pf: bool,
    skip_lifetime: bool,
    skip_spectrum: bool,
    cross_section: bool,
    strict: bool,
    config: dict[str, Any],
) -> str:
    """Run one complete validation and return PASS/PASS_WITH_WARNINGS/FAIL."""
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = output_dir / "run.log"
    file_handler = logging.FileHandler(log_path, mode="w", encoding="utf-8")
    file_handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    LOGGER.addHandler(file_handler)
    runtime_errors: list[str] = []
    try:
        metadata = parse_ion(ion, states_path, charge)
        states_result, states = validate_states(
            states_path,
            enforce_g_relation=bool(config.get("tolerances", {}).get("enforce_g_relation", True)),
        )
        trans_result = validate_transitions(
            trans_path,
            states,
            tolerance_cm1=float(config.get("tolerances", {}).get("wavenumber_absolute_cm1", 0.5)),
            relative_tolerance=float(config.get("tolerances", {}).get("wavenumber_relative", 1e-5)),
            chunk_size=int(config.get("chunk_size", 250_000)),
        )
        ground_g = float(states.loc[states["energy"].idxmin(), "g"]) if not states.empty else None
        pf_result, pf = validate_pf(
            pf_path,
            ground_degeneracy=ground_g,
            monotonic_rtol=float(config.get("tolerances", {}).get("pf_monotonic_rtol", 1e-8)),
            first_excited_energy_cm1=first_excited_energy(states),
        )
        _write_issues(output_dir / "states_issues.csv", states_result.issues)
        _write_issues(output_dir / "transitions_issues.csv", trans_result.issues)
        selected_temperatures = _temperatures(pf, temperatures)
        LOGGER.info("selected temperatures K: %s", selected_temperatures.tolist())

        staged = stage_local_files(
            states_path, trans_path, pf_path, metadata,
            output_dir / "_pyexocross_input",
            states_frame=states,
            transition_count=trans_result.metrics.get("rows"),
            max_wavenumber=trans_result.metrics.get("max_wavenumber_cm-1"),
            max_temperature=pf_result.metrics.get("temperature_max_K"),
        )
        adapter = PyExoCrossAdapter(
            staged, output_dir / "pyexocross_raw", log_path
        )
        pf_metrics: dict[str, object] = {"skipped": skip_pf}
        if not skip_pf:
            q_low_level = adapter.low_level_partition_function(selected_temperatures)
            q_public, raw_pf = adapter.calculate_partition_function(selected_temperatures)
            reference = reference_partition_function(
                states["energy"].to_numpy(float),
                states["g"].to_numpy(int),
                selected_temperatures,
            )
            comparison, comparison_metrics = compare_partition_functions(
                pf, selected_temperatures, q_public, reference
            )
            comparison.to_csv(output_dir / "pf_comparison.csv", index=False)
            save_partition_plots(comparison, output_dir)
            pf_metrics = {
                **comparison_metrics,
                "raw_pyexocross_pf": [str(path) for path in raw_pf],
                "low_level_public_api_max_abs_difference": float(
                    np.max(np.abs(q_low_level - q_public))
                ),
                "validation_class": "internal/circular consistency unless .pf provenance is independent",
            }
        else:
            pd.DataFrame().to_csv(output_dir / "pf_comparison.csv", index=False)

        lifetime_metrics: dict[str, object] = {"skipped": skip_lifetime}
        if not skip_lifetime:
            lifetime = calculate_lifetimes(
                states, trans_path, int(config.get("chunk_size", 250_000))
            )
            lifetime.to_csv(output_dir / "lifetime_comparison.csv", index=False)
            raw_lifetime = adapter.calculate_lifetimes()
            compared = lifetime["relative_error"].dropna()
            lifetime_metrics = {
                "state_count": len(lifetime),
                "compared_count": len(compared),
                "no_decay_channel_count": int(
                    (lifetime["status"] == "no_valid_downward_transition").sum()
                ),
                "median_absolute_relative_error": (
                    float(compared.abs().median()) if len(compared) else None
                ),
                "raw_pyexocross_files": [str(path) for path in raw_lifetime],
                "unit": "s",
            }
        else:
            pd.DataFrame().to_csv(output_dir / "lifetime_comparison.csv", index=False)

        spectrum_metrics: dict[str, object] = {"skipped": skip_spectrum}
        if not skip_spectrum:
            raw_spectrum = adapter.calculate_stick_spectrum(
                selected_temperatures.tolist(), wn_range[0], wn_range[1],
                str(config.get("spectrum", {}).get("mode", "Ab")),
            )
            spectrum_metrics = normalize_stick_output(
                raw_spectrum, output_dir, states_path
            )
            if not spectrum_metrics.get("generated"):
                runtime_errors.append("PyExoCross ran but stick output could not be normalized")
        else:
            pd.DataFrame().to_csv(output_dir / "strongest_lines.csv", index=False)

        cross_metrics: dict[str, object] = {"requested": cross_section}
        if cross_section:
            cross_cfg = config.get("cross_section", {})
            raw_cross = adapter.calculate_cross_section(
                float(selected_temperatures[1]),
                float(cross_cfg.get("pressure_bar", 0.0)),
                wn_range[0],
                wn_range[1],
                str(cross_cfg.get("profile", "Gaussian")),
                float(cross_cfg.get("bin_size_cm-1", 0.2)),
                cross_cfg.get("cutoff_cm-1"),
            )
            cross_metrics.update({
                "status": "exploratory / format-validation result",
                "raw_files": [str(path) for path in raw_cross],
                "output_unit": "cm2/molecule",
                "precision_claim": False,
            })

        results = [
            states_result.as_dict(), trans_result.as_dict(), pf_result.as_dict()
        ]
        overall = _status(results, runtime_errors)
        if strict and overall == "PASS_WITH_WARNINGS":
            overall = "FAIL"
        summary: dict[str, Any] = {
            "overall_status": overall,
            "metadata": metadata.as_dict(),
            "inputs": {
                "states": str(states_path.resolve()),
                "trans": str(trans_path.resolve()),
                "pf": str(pf_path.resolve()),
            },
            "pyexocross": {**pyexocross_info(), "calls": adapter.calls},
            "selected_temperatures_K": selected_temperatures.tolist(),
            "wavenumber_range_cm-1": list(wn_range),
            "states_validation": states_result.as_dict(),
            "transitions_validation": trans_result.as_dict(),
            "pf_validation": pf_result.as_dict(),
            "partition_comparison": pf_metrics,
            "lifetime_comparison": lifetime_metrics,
            "spectrum": spectrum_metrics,
            "cross_section": cross_metrics,
            "runtime_errors": runtime_errors,
        }
        write_report(output_dir, summary, {
            "Schema checks": json.dumps({
                "states": states_result.as_dict(),
                "transitions": trans_result.as_dict(),
                "partition_function": pf_result.as_dict(),
            }, indent=2),
            "Partition-function comparison": json.dumps(pf_metrics, indent=2),
            "Lifetime comparison": json.dumps(lifetime_metrics, indent=2),
            "LTE stick spectrum": json.dumps(spectrum_metrics, indent=2),
            "Cross-section": json.dumps(cross_metrics, indent=2),
            "Warnings": "\n".join(
                f"- {message}" for result in results for message in result["warnings"]
            ) or "None.",
        })
        return overall
    except (
        OSError, ValueError, PyExoCrossUnavailable, PyExoCrossExecutionError
    ) as exc:
        LOGGER.exception("validation failed")
        failure = {
            "overall_status": "FAIL",
            "inputs": {
                "states": str(states_path), "trans": str(trans_path), "pf": str(pf_path),
            },
            "error": f"{type(exc).__name__}: {exc}",
        }
        (output_dir / "validation_summary.json").write_text(
            json.dumps(failure, indent=2), encoding="utf-8"
        )
        (output_dir / "validation_report.md").write_text(
            f"# Validation failed\n\n{failure['error']}\n", encoding="utf-8"
        )
        return "FAIL"
    finally:
        LOGGER.removeHandler(file_handler)
        file_handler.close()


def _discover(input_root: Path) -> list[tuple[str, Path, Path, Path]]:
    datasets: list[tuple[str, Path, Path, Path]] = []
    for states in input_root.rglob("*.states"):
        base = states.name[:-7]
        trans = states.with_name(base + ".trans")
        pf = states.with_name(base + ".pf")
        if not (trans.exists() and pf.exists()):
            continue
        ion = states.parent.name if "-" in states.parent.name else base.split("__")[0].replace("_", "-")
        datasets.append((ion, states, trans, pf))
    return sorted(datasets)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Validate Kurucz ExoAtom files with PyExoCross."
    )
    parser.add_argument("--ion", help="spectroscopic ion, e.g. Sc-III")
    parser.add_argument("--charge", type=int, help="explicit charge cross-check")
    parser.add_argument("--states", type=Path)
    parser.add_argument("--trans", type=Path)
    parser.add_argument("--pf", type=Path)
    parser.add_argument("--input-root", type=Path)
    parser.add_argument("--all-ions", action="store_true")
    parser.add_argument("--temperature", type=float, action="append")
    parser.add_argument("--temperatures", nargs="+", type=float)
    parser.add_argument("--range", nargs=2, type=float, metavar=("MIN", "MAX"), default=(1000.0, 100000.0))
    parser.add_argument(
        "--range-unit", choices=["cm-1", "nm", "um"], default="cm-1",
        help="input range unit; wavelength ranges are vacuum wavelengths",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--skip-spectrum", action="store_true")
    parser.add_argument("--skip-lifetime", action="store_true")
    parser.add_argument("--skip-pf", action="store_true")
    parser.add_argument("--cross-section", action="store_true")
    parser.add_argument("--strict", action="store_true")
    parser.add_argument("--config", type=Path)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    try:
        config = _load_config(args.config)
        wn_range = _to_wavenumber_range(tuple(args.range), args.range_unit)
        requested_temperatures = args.temperatures or args.temperature
        if args.all_ions:
            if args.input_root is None:
                raise ValueError("--all-ions requires --input-root")
            datasets = _discover(args.input_root)
            if not datasets:
                raise ValueError("no complete .states/.trans/.pf datasets found")
        else:
            missing = [
                name for name in ["ion", "states", "trans", "pf"]
                if getattr(args, name) is None
            ]
            if missing:
                raise ValueError("single-ion mode requires --" + ", --".join(missing))
            datasets = [(args.ion, args.states, args.trans, args.pf)]
        statuses = []
        for ion, states, trans, pf in datasets:
            ion_output = args.output / ion if args.all_ions else args.output
            statuses.append(validate_one(
                ion=ion,
                states_path=states,
                trans_path=trans,
                pf_path=pf,
                output_dir=ion_output,
                charge=args.charge,
                temperatures=requested_temperatures,
                wn_range=wn_range,
                skip_pf=args.skip_pf,
                skip_lifetime=args.skip_lifetime,
                skip_spectrum=args.skip_spectrum,
                cross_section=args.cross_section,
                strict=args.strict,
                config=config,
            ))
        LOGGER.info("completed %d ion(s): %s", len(statuses), statuses)
        return 1 if "FAIL" in statuses else 0
    except (OSError, ValueError, yaml.YAMLError) as exc:
        LOGGER.error("%s", exc)
        return 2


if __name__ == "__main__":
    sys.exit(main())
