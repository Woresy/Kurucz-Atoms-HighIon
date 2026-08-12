"""Streaming parsers and preflight validation for ExoAtom text files."""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from pathlib import Path
from typing import Iterator

import numpy as np
import pandas as pd


STATES_COLUMNS = [
    "id", "energy", "g", "J", "uncertainty", "lifetime",
    "lande_g", "configuration", "term", "source_flag",
]
TRANS_COLUMNS = ["upper_id", "lower_id", "A", "wavenumber"]
PF_COLUMNS = ["temperature", "Q"]


@dataclass
class ValidationResult:
    """Machine-readable result shared by every preflight validator."""

    passed: bool = True
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    metrics: dict[str, object] = field(default_factory=dict)
    issues: list[dict[str, object]] = field(default_factory=list)

    def error(self, message: str, **issue: object) -> None:
        self.passed = False
        self.errors.append(message)
        if issue:
            self.issues.append({"severity": "error", "message": message, **issue})

    def warning(self, message: str, **issue: object) -> None:
        self.warnings.append(message)
        if issue:
            self.issues.append({"severity": "warning", "message": message, **issue})

    def as_dict(self) -> dict[str, object]:
        return {
            "passed": self.passed,
            "errors": self.errors,
            "warnings": self.warnings,
            "metrics": self.metrics,
        }


def read_states(path: Path) -> pd.DataFrame:
    """Read the project's headerless 10-column states schema."""
    records: list[list[object]] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, 1):
            parts = line.split()
            if not parts:
                continue
            if len(parts) < 10:
                raise ValueError(f"{path}:{line_number}: expected >=10 fields")
            configuration = " ".join(parts[7:-2])
            records.append(parts[:7] + [configuration, parts[-2], parts[-1]])
    frame = pd.DataFrame(records, columns=STATES_COLUMNS)
    for column in ["energy", "J", "uncertainty", "lifetime", "lande_g"]:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    for column in ["id", "g"]:
        numeric = pd.to_numeric(frame[column], errors="coerce")
        frame[column] = numeric.astype("Int64")
    return frame


def iter_transitions(path: Path, chunk_size: int = 250_000) -> Iterator[pd.DataFrame]:
    """Yield transitions in bounded-memory chunks."""
    yield from pd.read_csv(
        path,
        sep=r"\s+",
        header=None,
        names=TRANS_COLUMNS,
        chunksize=chunk_size,
        dtype={"upper_id": "Int64", "lower_id": "Int64", "A": float, "wavenumber": float},
    )


def read_pf(path: Path) -> pd.DataFrame:
    """Read a headerless temperature/partition-function file."""
    return pd.read_csv(
        path, sep=r"\s+", header=None, names=PF_COLUMNS,
        dtype={"temperature": float, "Q": float},
    )


def validate_states(path: Path, enforce_g_relation: bool = True) -> tuple[ValidationResult, pd.DataFrame]:
    """Validate IDs, energies, quantum numbers and optional columns."""
    result = ValidationResult()
    try:
        states = read_states(path)
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        result.error(f"states parse failed: {exc}")
        return result, pd.DataFrame(columns=STATES_COLUMNS)
    result.metrics.update(
        rows=len(states),
        schema="id E[cm^-1] g J uncertainty[cm^-1] lifetime[s] lande_g configuration term flag",
        delimiter="whitespace/fixed-width compatible",
        header=False,
    )
    if states.empty:
        result.error("states file is empty")
        return result, states
    ids = states["id"]
    if ids.isna().any():
        result.error("state ID contains non-integer or missing values")
    valid_ids = ids.dropna().astype(int)
    duplicates = int(valid_ids.duplicated().sum())
    if duplicates:
        result.error(f"{duplicates} duplicate state IDs")
    expected = np.arange(1, len(valid_ids) + 1)
    continuous = len(valid_ids) == len(states) and np.array_equal(np.sort(valid_ids), expected)
    result.metrics["ids_contiguous_from_one"] = continuous
    if not continuous:
        result.error("state IDs are not continuous from 1")
    finite_energy = np.isfinite(states["energy"])
    if not finite_energy.all():
        result.error(f"{int((~finite_energy).sum())} non-finite energies")
    negative = int((states["energy"] < 0).sum())
    if negative:
        result.error(f"{negative} negative state energies")
    ground = float(states["energy"].min())
    result.metrics["ground_energy_cm-1"] = ground
    if abs(ground) > 1e-3:
        result.warning(f"ground-state energy is {ground:g} cm^-1, not near zero")
    sorted_energy = bool(states["energy"].is_monotonic_increasing)
    result.metrics["energy_sorted"] = sorted_energy
    if not sorted_energy:
        result.warning("state energies are not non-decreasing")
    invalid_g = states["g"].isna() | (states["g"] <= 0)
    if invalid_g.any():
        result.error(f"{int(invalid_g.sum())} non-positive/invalid degeneracies")
    doubled_j = 2 * states["J"]
    legal_j = np.isfinite(states["J"]) & (states["J"] >= 0) & np.isclose(doubled_j, np.round(doubled_j), atol=1e-8)
    if not legal_j.all():
        result.error(f"{int((~legal_j).sum())} J values are not non-negative integer/half-integer")
    if enforce_g_relation:
        comparable = legal_j & states["g"].notna()
        mismatch = comparable & (states["g"].astype(float) != (2 * states["J"] + 1))
        result.metrics["g_equals_2J_plus_1"] = bool(not mismatch.any())
        if mismatch.any():
            result.error(f"{int(mismatch.sum())} rows violate project-defined g = 2J + 1")
    bad_unc = (~np.isfinite(states["uncertainty"])) | (states["uncertainty"] < 0)
    if bad_unc.any():
        result.error(f"{int(bad_unc.sum())} invalid uncertainties")
    lifetime = states["lifetime"]
    missing_lifetime = lifetime.isna()
    bad_lifetime = (~missing_lifetime) & (~np.isinf(lifetime)) & (lifetime <= 0)
    result.metrics["missing_lifetime_count"] = int(missing_lifetime.sum())
    result.metrics["nonpositive_lifetime_count"] = int(bad_lifetime.sum())
    if bad_lifetime.any():
        result.warning(
            f"{int(bad_lifetime.sum())} non-positive lifetimes; treated as unavailable"
        )
    blank_text = (
        states["configuration"].astype(str).str.strip().eq("")
        | states["term"].astype(str).str.strip().eq("")
    )
    if blank_text.any():
        result.error(f"{int(blank_text.sum())} blank configuration/term values")
    allowed_flags = {"NI", "CA"}
    unknown_flags = sorted(set(states["source_flag"]) - allowed_flags)
    if unknown_flags:
        result.warning(f"unrecognised predicted/measured flags: {unknown_flags}")
    result.metrics["energy_unit"] = "cm^-1 (project schema)"
    result.metrics["lifetime_unit"] = "s"
    return result, states


def validate_transitions(
    path: Path,
    states: pd.DataFrame,
    tolerance_cm1: float = 0.5,
    relative_tolerance: float = 1e-5,
    chunk_size: int = 250_000,
) -> ValidationResult:
    """Validate references and energy/wavenumber consistency by chunks."""
    result = ValidationResult()
    if states.empty or states["id"].isna().any():
        result.error("valid states are required before transition validation")
        return result
    energy_by_id = states.set_index("id")["energy"]
    valid_ids = set(int(value) for value in states["id"])
    total = missing_ref = invalid_a = invalid_wn = self_count = reversed_energy = 0
    sum_sq = max_error = 0.0
    max_wavenumber = 0.0
    above = 0
    duplicate_count = 0
    conflicting_a = 0
    seen: dict[tuple[int, int], float] = {}
    try:
        chunks = iter_transitions(path, chunk_size)
        for chunk_number, chunk in enumerate(chunks, 1):
            total += len(chunk)
            present = chunk["upper_id"].isin(valid_ids) & chunk["lower_id"].isin(valid_ids)
            missing_ref += int((~present).sum())
            finite_a = np.isfinite(chunk["A"]) & (chunk["A"] > 0)
            finite_wn = np.isfinite(chunk["wavenumber"]) & (chunk["wavenumber"] > 0)
            invalid_a += int((~finite_a).sum())
            invalid_wn += int((~finite_wn).sum())
            self_count += int((chunk["upper_id"] == chunk["lower_id"]).sum())
            usable = chunk[present & finite_a & finite_wn].copy()
            if usable.empty:
                continue
            max_wavenumber = max(max_wavenumber, float(usable["wavenumber"].max()))
            upper_e = usable["upper_id"].map(energy_by_id)
            lower_e = usable["lower_id"].map(energy_by_id)
            reversed_energy += int((upper_e <= lower_e).sum())
            calculated = (upper_e - lower_e).abs()
            abs_error = (usable["wavenumber"] - calculated).abs()
            rel_error = abs_error / calculated.replace(0, np.nan)
            sum_sq += float(np.square(abs_error).sum())
            max_error = max(max_error, float(abs_error.max()))
            bad = (abs_error > tolerance_cm1) & (rel_error > relative_tolerance)
            above += int(bad.sum())
            for row_index in usable.index[bad][:100]:
                result.issues.append({
                    "severity": "error",
                    "message": "transition wavenumber mismatch",
                    "row": int(row_index),
                    "absolute_error_cm-1": float(abs_error.loc[row_index]),
                    "relative_error": float(rel_error.loc[row_index]),
                })
            for upper, lower, a_value in usable[["upper_id", "lower_id", "A"]].itertuples(index=False):
                key = (int(upper), int(lower))
                previous = seen.get(key)
                if previous is not None:
                    duplicate_count += 1
                    scale = max(abs(previous), abs(float(a_value)), 1.0)
                    if abs(previous - float(a_value)) / scale > 1e-6:
                        conflicting_a += 1
                else:
                    seen[key] = float(a_value)
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        result.error(f"transitions parse failed: {exc}")
        return result
    result.metrics.update({
        "rows": total,
        "referenced_state_missing_count": missing_ref,
        "invalid_A_count": invalid_a,
        "invalid_wavenumber_count": invalid_wn,
        "self_transition_count": self_count,
        "upper_energy_not_higher_count": reversed_energy,
        "duplicate_transition_count": duplicate_count,
        "conflicting_A_count": conflicting_a,
        "wavenumber_rms_error_cm-1": math.sqrt(sum_sq / total) if total else math.nan,
        "wavenumber_max_error_cm-1": max_error if total else math.nan,
        "max_wavenumber_cm-1": max_wavenumber if total else math.nan,
        "wavenumber_above_tolerance_count": above,
        "A_unit": "s^-1",
        "wavenumber_unit": "cm^-1",
    })
    for count, label in [
        (missing_ref, "transitions reference missing states"),
        (invalid_a, "invalid Einstein A coefficients"),
        (invalid_wn, "invalid wavenumbers"),
        (self_count, "self-transitions"),
        (reversed_energy, "transitions whose declared upper state is not higher in energy"),
        (above, "wavenumbers exceed consistency tolerance"),
        (conflicting_a, "duplicate transitions have conflicting A coefficients"),
    ]:
        if count:
            result.error(f"{count} {label}")
    if duplicate_count and not conflicting_a:
        result.warning(f"{duplicate_count} duplicate upper/lower transition pairs")
    return result


def first_excited_energy(states: pd.DataFrame) -> float | None:
    """Return the lowest energy strictly above the ground state, if any."""
    if states.empty:
        return None
    energies = states["energy"].to_numpy(float)
    above_ground = energies[energies > energies.min()]
    return float(above_ground.min()) if above_ground.size else None


def validate_pf(
    path: Path,
    ground_degeneracy: float | None = None,
    monotonic_rtol: float = 1e-8,
    first_excited_energy_cm1: float | None = None,
) -> tuple[ValidationResult, pd.DataFrame]:
    """Validate the temperature grid and partition-function values."""
    result = ValidationResult()
    try:
        pf = read_pf(path)
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        result.error(f"partition-function parse failed: {exc}")
        return result, pd.DataFrame(columns=PF_COLUMNS)
    result.metrics.update(rows=len(pf), temperature_unit="K", Q_unit="dimensionless")
    if pf.empty:
        result.error("partition-function file is empty")
        return result, pf
    invalid_t = (~np.isfinite(pf["temperature"])) | (pf["temperature"] <= 0)
    invalid_q = (~np.isfinite(pf["Q"])) | (pf["Q"] <= 0)
    if invalid_t.any():
        result.error(f"{int(invalid_t.sum())} invalid temperatures")
    if invalid_q.any():
        result.error(f"{int(invalid_q.sum())} invalid partition functions")
    duplicate_t = int(pf["temperature"].duplicated().sum())
    if duplicate_t:
        result.error(f"{duplicate_t} duplicate temperatures")
    if not pf["temperature"].is_monotonic_increasing:
        result.error("temperatures are not strictly increasing")
    delta_q = pf["Q"].diff().iloc[1:]
    scale = np.maximum(pf["Q"].iloc[:-1].to_numpy(), 1.0)
    decreasing = delta_q.to_numpy() < -monotonic_rtol * scale
    if decreasing.any():
        result.error(f"{int(decreasing.sum())} significant decreases in Q(T)")
    spacing = pf["temperature"].diff().dropna().to_numpy()
    uniform = bool(len(spacing) < 2 or np.allclose(spacing, spacing[0]))
    result.metrics.update(
        temperature_min_K=float(pf["temperature"].min()),
        temperature_max_K=float(pf["temperature"].max()),
        uniform_temperature_grid=uniform,
    )
    if not uniform:
        result.warning("temperature grid is non-uniform (valid, interpolation is required)")
    if ground_degeneracy is not None:
        ratio = float(pf["Q"].iloc[0] / ground_degeneracy)
        result.metrics["low_temperature_Q_over_ground_g"] = ratio
        # Q(T) sums positive terms of which the ground state contributes exactly
        # g_ground, so Q below g_ground is unphysical at any temperature.
        if ratio < 0.8:
            result.warning("lowest-temperature Q is below the ground-state degeneracy")
        elif ratio > 1.5:
            # Q only collapses onto g_ground once the first excited level is
            # thermally inaccessible. Every ion whose ground term carries fine
            # structure a few tens of cm^-1 up keeps a legitimately larger ratio
            # at 100 K, so this half of the check means nothing until that
            # level's Boltzmann population is negligible.
            from .calculations import C2_CM_K  # deferred: calculations imports this module

            lowest_temperature = float(pf["temperature"].iloc[0])
            population = (
                math.exp(-C2_CM_K * first_excited_energy_cm1 / lowest_temperature)
                if first_excited_energy_cm1 is not None
                else None
            )
            if population is None:
                result.warning(
                    "lowest-temperature Q exceeds the ground-state degeneracy "
                    "(no state list supplied to check for low-lying levels)"
                )
            else:
                result.metrics["first_excited_population_at_lowest_T"] = population
                if population < 0.01:
                    result.warning(
                        "lowest-temperature Q exceeds the ground-state degeneracy "
                        "without a thermally populated low-lying level"
                    )
    if float(pf["temperature"].max()) < 10_000:
        result.warning("temperature coverage may be low for a highly ionized species")
    return result, pf
