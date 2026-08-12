"""Independent reference calculations and comparisons."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd

from .exomol_validator import iter_transitions

# Second radiation constant c2 = h*c/k_B in cm*K. Since the 2019 SI revision h,
# c, and k_B are all exact, so this value is exact rather than measured. The
# previously used 1.438776877 was the same constant truncated to the CODATA-2018
# published digits; that truncation is a 3.5e-10 relative error, which Boltzmann
# weighting amplifies to ~1e-9 in Q and tripped the adapter's 1e-9 agreement
# assertion against PyExoCross, which carries the full-precision value.
C2_CM_K = 1.4387768775039338


def reference_partition_function(
    energies_cm1: np.ndarray, degeneracies: np.ndarray, temperatures_k: np.ndarray
) -> np.ndarray:
    """Compute Q(T) independently for sanity checking."""
    energy = np.asarray(energies_cm1, dtype=float)
    g = np.asarray(degeneracies, dtype=float)
    temperatures = np.asarray(temperatures_k, dtype=float)
    if np.any(temperatures <= 0):
        raise ValueError("temperatures must be positive")
    return np.array(
        [np.sum(g * np.exp(-C2_CM_K * energy / temperature)) for temperature in temperatures]
    )


def interpolate_in_bounds(
    source_x: np.ndarray, source_y: np.ndarray, target_x: np.ndarray
) -> np.ndarray:
    """Linearly interpolate, rejecting silent extrapolation."""
    x = np.asarray(source_x, dtype=float)
    targets = np.asarray(target_x, dtype=float)
    if targets.size and (targets.min() < x.min() or targets.max() > x.max()):
        raise ValueError("target temperature lies outside the source range")
    return np.interp(targets, x, np.asarray(source_y, dtype=float))


def compare_partition_functions(
    pf_file: pd.DataFrame,
    temperatures: np.ndarray,
    q_pyexocross: np.ndarray,
    q_reference: np.ndarray,
) -> tuple[pd.DataFrame, dict[str, float]]:
    """Compare file, PyExoCross, and independent reference values."""
    q_file = interpolate_in_bounds(
        pf_file["temperature"].to_numpy(), pf_file["Q"].to_numpy(), temperatures
    )
    result = pd.DataFrame({
        "temperature_K": temperatures,
        "Q_file": q_file,
        "Q_pyexocross": q_pyexocross,
        "Q_reference": q_reference,
    })
    result["absolute_difference"] = result["Q_pyexocross"] - result["Q_file"]
    result["relative_difference"] = (
        result["absolute_difference"] / result["Q_file"]
    )
    result["log_ratio"] = np.log(result["Q_pyexocross"] / result["Q_file"])
    rel = result["relative_difference"].to_numpy()
    metrics = {
        "mean_relative_error": float(np.mean(np.abs(rel))),
        "rms_relative_error": float(math.sqrt(np.mean(np.square(rel)))),
        "maximum_relative_error": float(np.max(np.abs(rel))),
        "pyexocross_reference_max_relative_error": float(np.max(
            np.abs(q_pyexocross - q_reference) / np.maximum(np.abs(q_reference), 1e-300)
        )),
    }
    return result, metrics


def calculate_lifetimes(
    states: pd.DataFrame, trans_path: Path, chunk_size: int = 250_000
) -> pd.DataFrame:
    """Calculate 1/sum(A) from every spontaneous upper-state channel."""
    sums: dict[int, float] = {}
    channels: dict[int, int] = {}
    for chunk in iter_transitions(trans_path, chunk_size):
        valid = chunk[np.isfinite(chunk["A"]) & (chunk["A"] > 0)]
        grouped = valid.groupby("upper_id")["A"].agg(["sum", "count"])
        for state_id, row in grouped.iterrows():
            key = int(state_id)
            sums[key] = sums.get(key, 0.0) + float(row["sum"])
            channels[key] = channels.get(key, 0) + int(row["count"])
    output = states[["id", "energy", "lifetime"]].copy()
    output["sum_A_s-1"] = output["id"].map(sums).fillna(0.0)
    output["decay_channels"] = output["id"].map(channels).fillna(0).astype(int)
    output["calculated_lifetime_s"] = np.where(
        output["sum_A_s-1"] > 0, 1.0 / output["sum_A_s-1"], np.nan
    )
    original_valid = np.isfinite(output["lifetime"]) & (output["lifetime"] > 0)
    output["relative_error"] = np.where(
        original_valid,
        (output["calculated_lifetime_s"] - output["lifetime"]) / output["lifetime"],
        np.nan,
    )
    output["status"] = np.where(
        output["sum_A_s-1"] <= 0,
        "no_valid_downward_transition",
        np.where(original_valid, "compared", "original_lifetime_missing"),
    )
    return output
