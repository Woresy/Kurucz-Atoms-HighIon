"""Generic higher-ionization Kurucz/NIST/PyExoCross validation CLI."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import shutil
import sys

import numpy as np
import pandas as pd

from .cross_source import (
    compare_partition_sources,
    compare_transitions,
    match_states,
    read_nist_states,
)
from .cross_source_plots import generate_cross_source_plots
from .discovery import Availability, discover_complete_kurucz_ions, discover_ion_files
from .exomol_validator import (
    ValidationResult,
    first_excited_energy,
    iter_transitions,
    read_pf,
    read_states,
    validate_pf,
    validate_states,
)
from .metadata import parse_ion
from .pyexocross_adapter import PyExoCrossAdapter, pyexocross_info, stage_local_files
from .reporting import normalize_stick_output


DEFAULT_KURUCZ = Path("Kurucz-Nist-Overlap-data")
DEFAULT_NIST = Path("Nist-temp-data")
DEFAULT_OUTPUT = Path("validation/reports")


def _json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, allow_nan=True), encoding="utf-8")


def _transition_format(path: Path, states: pd.DataFrame) -> ValidationResult:
    """Bounded-memory transition format/reference validation for huge lists."""
    result = ValidationResult()
    ids = set(states["id"].dropna().astype(int))
    total = missing = invalid_a = invalid_wn = 0
    minimum_wavenumber = math.inf
    maximum_wavenumber = -math.inf
    for chunk in iter_transitions(path):
        total += len(chunk)
        missing += int((~chunk.upper_id.isin(ids) | ~chunk.lower_id.isin(ids)).sum())
        invalid_a += int((~np.isfinite(chunk.A) | (chunk.A <= 0)).sum())
        invalid_wn += int((~np.isfinite(chunk.wavenumber) | (chunk.wavenumber <= 0)).sum())
        positive = chunk.loc[
            np.isfinite(chunk.wavenumber) & (chunk.wavenumber > 0),
            "wavenumber",
        ]
        if len(positive):
            minimum_wavenumber = min(minimum_wavenumber, float(positive.min()))
            maximum_wavenumber = max(maximum_wavenumber, float(positive.max()))
    result.metrics.update({
        "rows": total,
        "referenced_state_missing_count": missing,
        "invalid_A_count": invalid_a,
        "invalid_wavenumber_count": invalid_wn,
        "minimum_wavenumber_cm-1": (
            minimum_wavenumber
            if math.isfinite(minimum_wavenumber) else math.nan
        ),
        "maximum_wavenumber_cm-1": (
            maximum_wavenumber
            if math.isfinite(maximum_wavenumber) else math.nan
        ),
    })
    if missing:
        result.error(f"{missing} transitions reference missing states")
    if invalid_a:
        result.error(f"{invalid_a} invalid Einstein A coefficients")
    if invalid_wn:
        result.error(f"{invalid_wn} invalid wavenumbers")
    return result


def _write_nist_compat(states: pd.DataFrame, path: Path) -> None:
    """Create a disposable ten-column ExoAtom view for PyExoCross 1.1.9."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for row in states.itertuples(index=False):
            lande = row.lande_g if np.isfinite(row.lande_g) else 0.0
            configuration = ".".join(str(row.configuration).split())
            term = ".".join(str(row.term).split())
            handle.write(
                f"{int(row.id):12d} {row.energy:16.6f} {int(row.g):6d} "
                f"{row.J:7.1f} {row.uncertainty:12.6f} {float('inf'):12g} "
                f"{lande:10.6f} {configuration:<24s} {term:<12s} NI\n"
            )


def _run_pyexocross(
    ion: str,
    source: str,
    states_path: Path,
    trans_path: Path,
    pf_path: Path,
    states: pd.DataFrame,
    output: Path,
    temperatures: list[float],
    wavelength_range_nm: tuple[float, float],
    parity_states_path: Path | None = None,
) -> dict[str, object]:
    source_dir = output / f"pyexocross_{source.lower()}"
    metadata = parse_ion(ion, states_path, source_database=source)
    staged = stage_local_files(
        states_path, trans_path, pf_path, metadata,
        source_dir / "_input", dataset=source,
    )
    adapter = PyExoCrossAdapter(staged, source_dir / "raw", source_dir / "run.log")
    temps = np.asarray(temperatures, dtype=float)
    q, raw_pf = adapter.calculate_partition_function(temps)
    pf_frame = read_pf(pf_path)
    q_file = np.interp(temps, pf_frame.temperature, pf_frame.Q)
    comparison = pd.DataFrame({
        "temperature_K": temps, f"Q_{source.lower()}_file": q_file,
        f"Q_{source.lower()}_pyexocross": q,
        "relative_difference": (q - q_file) / q_file,
    })
    comparison.to_csv(output / f"pyexocross_{source.lower()}_pf.csv", index=False)
    low_nm, high_nm = wavelength_range_nm
    raw_spectrum = adapter.calculate_stick_spectrum(
        temperatures, 1e7 / high_nm, 1e7 / low_nm
    )
    normalized = normalize_stick_output(
        raw_spectrum, source_dir, parity_states_path or states_path
    )
    spectrum = source_dir / "stick_spectrum.csv"
    target = output / f"pyexocross_{source.lower()}_spectrum.csv"
    if spectrum.exists():
        if target.exists():
            target.unlink()
        try:
            target.hardlink_to(spectrum)
        except OSError:
            shutil.copy2(spectrum, target)
    return {
        "status": "SUCCESS",
        "version": pyexocross_info()["version"],
        "pf_raw": str(raw_pf),
        "pf_relative_RMS": float(np.sqrt(np.mean(np.square(comparison.relative_difference)))),
        "spectrum": normalized,
        "calls": adapter.calls,
    }


def validate_ion(
    ion: str,
    kurucz_root: Path,
    nist_root: Path,
    output_root: Path,
    temperatures: list[float],
    wavelength_range_nm: tuple[float, float] | None,
    skip_pyexocross: bool = False,
) -> dict[str, object]:
    resolved = discover_ion_files(ion, kurucz_root, nist_root)
    output = output_root / ion.replace("_", "-")
    output.mkdir(parents=True, exist_ok=True)
    kf, nf = resolved.kurucz_files, resolved.nist_files
    if kf.states is None:
        raise FileNotFoundError(f"Kurucz states not available for {ion}")

    k_states_result, k_states = validate_states(kf.states)
    _json(output / "kurucz_states_validation.json", k_states_result.as_dict())
    k_trans_result = _transition_format(kf.trans, k_states) if kf.trans else ValidationResult(False, ["NOT_AVAILABLE"])
    _json(output / "kurucz_trans_validation.json", k_trans_result.as_dict())
    k_pf_result, _ = validate_pf(
        kf.pf,
        float(k_states.loc[k_states.energy.idxmin(), "g"]),
        first_excited_energy_cm1=first_excited_energy(k_states),
    ) if kf.pf else (ValidationResult(False, ["NOT_AVAILABLE"]), pd.DataFrame())
    _json(output / "kurucz_pf_validation.json", k_pf_result.as_dict())

    n_states = pd.DataFrame()
    n_states_result = ValidationResult()
    if nf.states:
        try:
            n_states = read_nist_states(nf.states)
            n_states_result.metrics = {"rows": len(n_states), "schema": "NIST ExoAtom states"}
            if n_states.empty:
                n_states_result.error("NIST states file is empty")
        except (OSError, ValueError) as exc:
            n_states_result.error(f"NIST states parse failed: {exc}")
    else:
        n_states_result.warning("NOT_AVAILABLE")
    _json(output / "nist_states_validation.json", n_states_result.as_dict())
    n_trans_result = _transition_format(nf.trans, n_states) if nf.trans and not n_states.empty else ValidationResult()
    if not nf.trans:
        n_trans_result.warning("NOT_AVAILABLE")
    _json(output / "nist_trans_validation.json", n_trans_result.as_dict())
    if nf.pf:
        n_pf_result, _ = validate_pf(
            nf.pf,
            float(n_states.loc[n_states.energy.idxmin(), "g"])
            if not n_states.empty else None,
            first_excited_energy_cm1=first_excited_energy(n_states),
        )
    else:
        n_pf_result = ValidationResult()
        n_pf_result.warning("NOT_AVAILABLE")
    _json(output / "nist_pf_validation.json", n_pf_result.as_dict())

    if wavelength_range_nm is None:
        minima = [
            result.metrics.get("minimum_wavenumber_cm-1")
            for result in (k_trans_result, n_trans_result)
        ]
        maxima = [
            result.metrics.get("maximum_wavenumber_cm-1")
            for result in (k_trans_result, n_trans_result)
        ]
        minima = [float(value) for value in minima if value is not None and np.isfinite(value)]
        maxima = [float(value) for value in maxima if value is not None and np.isfinite(value)]
        if not minima or not maxima:
            raise ValueError(f"cannot determine a full wavelength range for {ion}")
        # Add a tiny numerical margin so transitions exactly at either endpoint
        # are retained by external implementations using strict inequalities.
        min_wn = min(minima) * (1 - 1e-10)
        max_wn = max(maxima) * (1 + 1e-10)
        effective_wavelength_range_nm = (1e7 / max_wn, 1e7 / min_wn)
    else:
        effective_wavelength_range_nm = wavelength_range_nm
        if (
            effective_wavelength_range_nm[0] <= 0
            or effective_wavelength_range_nm[1]
            <= effective_wavelength_range_nm[0]
        ):
            raise ValueError(
                "--wavelength-range-nm requires positive MIN < MAX"
            )

    summary: dict[str, object] = {
        "ion": ion.replace("_", "-"), "metadata": resolved.as_dict(),
        "element": resolved.element_symbol,
        "stage": resolved.spectroscopic_stage,
        "charge": resolved.charge,
        "availability": {
            "kurucz": kf.availability.value, "nist": nf.availability.value,
            "kurucz_components": kf.as_dict(), "nist_components": nf.as_dict(),
            "component_status": {
                "kurucz_states": (
                    "AVAILABLE" if k_states_result.passed else "INVALID"
                ),
                "kurucz_trans": (
                    "NOT_AVAILABLE" if kf.trans is None
                    else "AVAILABLE" if k_trans_result.passed else "INVALID"
                ),
                "kurucz_pf": (
                    "NOT_AVAILABLE" if kf.pf is None
                    else "AVAILABLE" if k_pf_result.passed else "INVALID"
                ),
                "nist_states": (
                    "NOT_AVAILABLE" if nf.states is None
                    else "AVAILABLE" if n_states_result.passed else "INVALID"
                ),
                "nist_trans": (
                    "NOT_AVAILABLE" if nf.trans is None
                    else "AVAILABLE" if n_trans_result.passed else "INVALID"
                ),
                "nist_pf": (
                    "NOT_AVAILABLE" if nf.pf is None
                    else "AVAILABLE" if n_pf_result.passed else "INVALID"
                ),
            },
        },
        "classification_notes": {
            "kurucz_only": "DATABASE_COVERAGE_DIFFERENCE or PREDICTED_KURUCZ_DATA; not automatically incorrect",
            "nist_only": "DATABASE_COVERAGE_DIFFERENCE or OBSERVED_NIST_DATA; not automatically a conversion error",
            "pyexocross_success": "compatibility evidence, not proof of scientific correctness",
        },
        "kurucz_states": len(k_states),
        "nist_states": len(n_states),
        "kurucz_transitions": k_trans_result.metrics.get("rows", math.nan),
        "nist_transitions": n_trans_result.metrics.get("rows", math.nan),
        "wavelength_range_nm": list(effective_wavelength_range_nm),
        "wavelength_range_mode": (
            "ALL_AVAILABLE_TRANSITIONS"
            if wavelength_range_nm is None else "USER_SELECTED"
        ),
    }

    if not n_states.empty:
        state_matches = match_states(k_states, n_states)
        state_matches.to_csv(output / "kurucz_vs_nist_states.csv", index=False)
        matched_state_rows = state_matches[
            state_matches["classification"] != "UNMATCHED"
        ]
        de = (
            matched_state_rows["delta_energy_cm-1"]
            if len(matched_state_rows) else pd.Series(dtype=float)
        )
        summary.update({
            "matched_states": len(matched_state_rows),
            "state_match_rate": len(matched_state_rows) / len(n_states),
            "energy_MAE_cm-1": float(de.abs().mean()) if len(de) else math.nan,
            "energy_RMS_cm-1": float(np.sqrt(np.mean(np.square(de)))) if len(de) else math.nan,
        })
        if nf.trans and kf.trans:
            trans_cmp, trans_metrics = compare_transitions(
                kf.trans, nf.trans, state_matches
            )
            trans_cmp.to_csv(output / "kurucz_vs_nist_transitions.csv", index=False)
            summary.update(trans_metrics)
            if len(trans_cmp):
                summary.update({
                    "wavelength_MAE_nm": float(trans_cmp.delta_wavelength_nm.abs().mean()),
                    "wavenumber_MAE_cm-1": float(trans_cmp["delta_wavenumber_cm-1"].abs().mean()),
                    "median_abs_delta_logA": float(trans_cmp.delta_logA.abs().median()),
                })
                trans_cmp.nlargest(50, "A_nist_s-1").to_csv(
                    output / "strongest_line_comparison.csv", index=False
                )
            else:
                pd.DataFrame().to_csv(
                    output / "strongest_line_comparison.csv", index=False
                )
    if not (output / "strongest_line_comparison.csv").exists():
        pd.DataFrame().to_csv(output / "strongest_line_comparison.csv", index=False)

    pd.DataFrame(columns=[
        "state_match_classification", "kurucz_id", "nist_id",
        "lifetime_kurucz_s", "lifetime_nist_s", "relative_difference",
        "status",
    ]).to_csv(output / "lifetime_comparison.csv", index=False)
    summary["lifetime_comparison_status"] = (
        "NOT_AVAILABLE: NIST states export contains no lifetime column"
    )

    if kf.pf and nf.pf:
        pf_cmp, pf_metrics = compare_partition_sources(kf.pf, nf.pf)
        pf_cmp.to_csv(output / "kurucz_vs_nist_pf.csv", index=False)
        summary.update(pf_metrics)
    else:
        pd.DataFrame().to_csv(output / "kurucz_vs_nist_pf.csv", index=False)
        summary["pf_relative_RMS"] = math.nan
        summary["pf_comparison_status"] = "NOT_AVAILABLE"

    for source, files, states in (("Kurucz", kf, k_states), ("NIST", nf, n_states)):
        key = f"pyexocross_{source.lower()}_status"
        if skip_pyexocross:
            summary[key] = "SKIPPED"
        elif files.availability is not Availability.AVAILABLE:
            summary[key] = "NOT_AVAILABLE"
        else:
            try:
                states_path = files.states
                if source == "NIST":
                    states_path = output / "_nist_compat" / files.states.name
                    _write_nist_compat(states, states_path)
                px = _run_pyexocross(
                    ion, source, states_path, files.trans, files.pf, states,
                    output, temperatures, effective_wavelength_range_nm,
                    # The ExoAtom-compatible copy has no parity column, so
                    # transition classification must read the original file.
                    parity_states_path=files.states,
                )
                summary[key] = "SUCCESS"
                summary[f"pyexocross_{source.lower()}"] = px
            except Exception as exc:
                summary[key] = "PYEXOCROSS_COMPATIBILITY_ISSUE"
                summary[f"pyexocross_{source.lower()}_error"] = f"{type(exc).__name__}: {exc}"

    ks = output / "pyexocross_kurucz_spectrum.csv"
    ns = output / "pyexocross_nist_spectrum.csv"
    if ks.exists() and ns.exists():
        kspec = pd.read_csv(ks, usecols=["vacuum_wavelength_nm"])
        nspec = pd.read_csv(ns, usecols=["vacuum_wavelength_nm"])
        kw = np.sort(kspec.vacuum_wavelength_nm.to_numpy(float))
        nw = nspec.vacuum_wavelength_nm.to_numpy(float)
        indices = np.searchsorted(kw, nw)
        left = np.abs(nw - kw[np.maximum(indices - 1, 0)])
        right = np.abs(nw - kw[np.minimum(indices, len(kw) - 1)])
        nearest = np.minimum(left, right)
        summary["spectrum_line_match_tolerance_nm"] = 0.01
        summary["spectrum_line_match_rate"] = float(np.mean(nearest <= 0.01))
        summary["spectrum_nist_line_count"] = len(nw)
        summary["spectrum_kurucz_line_count"] = len(kw)
    else:
        summary["spectrum_line_match_rate"] = math.nan

    kpf = output / "pyexocross_kurucz_pf.csv"
    npf = output / "pyexocross_nist_pf.csv"
    if kpf.exists() and npf.exists():
        kp = pd.read_csv(kpf)
        npx = pd.read_csv(npf)
        combined_pf = kp.merge(npx, on="temperature_K", suffixes=("_kurucz", "_nist"))
        combined_pf.to_csv(output / "pyexocross_pf_comparison.csv", index=False)
    else:
        pd.DataFrame().to_csv(output / "pyexocross_pf_comparison.csv", index=False)

    summary["plots"] = generate_cross_source_plots(
        output, ion.replace("_", "-")
    )

    failures = not k_states_result.passed or not k_trans_result.passed or not k_pf_result.passed
    summary["overall_status"] = "FAIL" if failures else "PASS_WITH_WARNINGS"
    _json(output / "validation_summary.json", summary)
    _write_report(output / "validation_report.md", summary)
    return summary


def _write_report(path: Path, summary: dict[str, object]) -> None:
    metadata = summary["metadata"]
    lines = [
        f"# Generic higher-ionization validation: {summary['ion']}", "",
        f"Overall status: **{summary['overall_status']}**", "",
        "## Identity", "",
        f"- Element: {metadata['element_symbol']} (Z={metadata['atomic_number']})",
        f"- Stage: {metadata['spectroscopic_stage']} ({metadata['roman_stage']})",
        f"- Charge: +{metadata['charge']}",
        f"- Bound electrons: {metadata['electron_count']}", "",
        "## Availability", "",
        f"- Kurucz: {summary['availability']['kurucz']}",
        f"- NIST: {summary['availability']['nist']}", "",
        "## Quantitative results", "",
    ]
    keys = [
        "kurucz_states", "nist_states", "matched_states", "state_match_rate",
        "energy_MAE_cm-1", "energy_RMS_cm-1", "kurucz_transitions",
        "nist_transitions", "matched_transitions", "transition_match_rate",
        "wavelength_MAE_nm", "wavenumber_MAE_cm-1",
        "median_abs_delta_logA", "pf_relative_RMS",
        "pyexocross_kurucz_status", "pyexocross_nist_status",
    ]
    lines.extend(f"- {key}: {summary.get(key, 'NOT_AVAILABLE')}" for key in keys)
    lines += [
        "", "## Interpretation rules", "",
        "- Format failures are distinguished from conversion, internal-consistency, database-disagreement, coverage, and PyExoCross compatibility issues.",
        "- Kurucz-only predicted data and NIST-only observed data are coverage differences, not automatic errors.",
        "- PyExoCross success demonstrates executable compatibility, not scientific correctness.",
        "", "## Remaining assumptions", "",
        "- Text labels are normalized generically; principal-shell digits and punctuation are ignored for configuration comparison.",
        "- Missing parity in Kurucz data weakens, rather than prevents, a match.",
        "- Transition matching is performed only through the one-to-one matched-state mapping.",
    ]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--ion")
    group.add_argument("--ions", nargs="+")
    group.add_argument("--all-ions", action="store_true")
    parser.add_argument("--kurucz-root", type=Path, default=DEFAULT_KURUCZ)
    parser.add_argument("--nist-root", type=Path, default=DEFAULT_NIST)
    parser.add_argument("--output-root", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--temperatures", nargs="+", type=float, default=[1000, 5000, 10000])
    wavelength = parser.add_mutually_exclusive_group()
    wavelength.add_argument(
        "--wavelength-range-nm", nargs=2, type=float, metavar=("MIN", "MAX")
    )
    wavelength.add_argument(
        "--all-wavelengths", action="store_true",
        help="derive the full positive wavelength range from all available "
             "transitions (this is the default; the flag states it explicitly)",
    )
    parser.add_argument("--skip-pyexocross", action="store_true")
    return parser


def write_batch_summary(output_root: Path) -> Path:
    """Rebuild summary.csv from every completed ion report under a root."""
    rows: list[dict[str, object]] = []
    for summary_path in sorted(output_root.glob("*/validation_summary.json")):
        try:
            item = json.loads(summary_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if not isinstance(item, dict) or "ion" not in item:
            continue
        rows.append({
            key: value for key, value in item.items()
            if not isinstance(value, dict)
        })
    path = output_root / "summary.csv"
    pd.DataFrame(rows).sort_values("ion").to_csv(path, index=False)
    return path


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    ions = (
        discover_complete_kurucz_ions(args.kurucz_root) if args.all_ions
        else args.ions or [args.ion]
    )
    summaries = []
    for ion in ions:
        # Default to the full range covered by the transitions themselves. A
        # fixed window is nearly always the wrong one: NIST publishes only
        # evaluated lines, and for most ions none of them fall in any given
        # narrow band, which leaves nothing to compare against Kurucz.
        wavelength_range = (
            tuple(args.wavelength_range_nm) if args.wavelength_range_nm else None
        )
        summaries.append(validate_ion(
            ion, args.kurucz_root, args.nist_root, args.output_root,
            args.temperatures, wavelength_range,
            args.skip_pyexocross,
        ))
    write_batch_summary(args.output_root)
    return 1 if any(item["overall_status"] == "FAIL" for item in summaries) else 0


if __name__ == "__main__":
    sys.exit(main())
