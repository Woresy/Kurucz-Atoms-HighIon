"""Generic Kurucz/NIST parsing, matching, and quantitative comparison."""

from __future__ import annotations

import math
import re
from pathlib import Path

import numpy as np
import pandas as pd

from .exomol_validator import TRANS_COLUMNS, iter_transitions, read_pf, read_states


NIST_COLUMNS = [
    "id", "energy", "g", "J", "uncertainty", "lande_g",
    "configuration", "term", "parity",
]


def read_nist_states(path: Path) -> pd.DataFrame:
    """Parse both NIST state variants (with or without Landé g)."""
    rows: list[list[object]] = []
    with path.open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            tokens = line.split()
            if not tokens:
                continue
            numeric = 0
            for token in tokens:
                try:
                    float(token)
                    numeric += 1
                except ValueError:
                    break
            if numeric not in (5, 6) or len(tokens) < numeric + 3:
                raise ValueError(f"{path}:{number}: unsupported NIST states row")
            lead = tokens[:numeric]
            if numeric == 5:
                lead.append(np.nan)
            rows.append(lead + [" ".join(tokens[numeric:-2]), tokens[-2], tokens[-1]])
    frame = pd.DataFrame(rows, columns=NIST_COLUMNS)
    for column in ["energy", "g", "J", "uncertainty", "lande_g"]:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    frame["id"] = pd.to_numeric(frame["id"], errors="raise").astype(int)
    frame["lifetime"] = np.nan
    frame["source_flag"] = "NIST"
    return frame


def normalize_label(value: object) -> str:
    return re.sub(r"[^a-z0-9]", "", str(value).lower())


def normalize_configuration(value: object) -> str:
    value = normalize_label(value)
    return re.sub(r"(?<![a-z])\d+(?=[spdfgh])", "", value)


def _parity_from_configuration(value: object) -> str:
    text = str(value).lower()
    if "°" in text or text.endswith("o") or "*" in text:
        return "-"
    return ""


def prepare_states(frame: pd.DataFrame, source: str) -> pd.DataFrame:
    out = frame.copy()
    out["configuration_original"] = out["configuration"].astype(str)
    out["term_original"] = out["term"].astype(str)
    out["configuration_norm"] = out["configuration"].map(normalize_configuration)
    out["term_norm"] = out["term"].map(normalize_label)
    if "parity" not in out:
        out["parity"] = out["configuration"].map(_parity_from_configuration)
    out["parity_norm"] = out["parity"].astype(str).str.strip().replace({"+": "+", "-": "-"})
    out["source"] = source
    return out


def match_states(
    kurucz: pd.DataFrame,
    nist: pd.DataFrame,
    energy_tolerance_cm1: float = 5.0,
) -> pd.DataFrame:
    """One-to-one quantum-aware matching independent of database state IDs."""
    k = prepare_states(kurucz, "Kurucz")
    n = prepare_states(nist, "NIST")
    candidates: list[tuple[float, int, int, int, bool, bool, bool]] = []
    by_j = {float(j): sub for j, sub in k.groupby("J", dropna=True)}
    for ni, row in n.iterrows():
        sub = by_j.get(float(row["J"]))
        if sub is None:
            continue
        delta = (sub["energy"] - row["energy"]).abs()
        for ki in delta.nsmallest(8).index:
            de = float(delta.at[ki])
            if de > energy_tolerance_cm1:
                continue
            kr = k.loc[ki]
            config = bool(row["configuration_norm"] and row["configuration_norm"] == kr["configuration_norm"])
            term = bool(row["term_norm"] and row["term_norm"] == kr["term_norm"])
            parity = bool(
                not row["parity_norm"] or not kr["parity_norm"]
                or row["parity_norm"] == kr["parity_norm"]
            )
            score = de + (0 if term else 20) + (0 if config else 10) + (0 if parity else 50)
            candidates.append((score, ni, ki, de, config, term, parity))
    candidates.sort()
    used_n: set[int] = set()
    used_k: set[int] = set()
    rows: list[dict[str, object]] = []
    for _, ni, ki, de, config, term, parity in candidates:
        if ni in used_n or ki in used_k:
            continue
        used_n.add(ni)
        used_k.add(ki)
        nr, kr = n.loc[ni], k.loc[ki]
        exact = config and term and parity
        classification = (
            "EXACT_QUANTUM_MATCH" if exact and de <= 0.1
            else "HIGH_CONFIDENCE" if term and parity
            else "AMBIGUOUS"
        )
        rows.append({
            "nist_id": int(nr["id"]), "kurucz_id": int(kr["id"]),
            "classification": classification,
            "energy_nist_cm-1": nr["energy"], "energy_kurucz_cm-1": kr["energy"],
            "delta_energy_cm-1": kr["energy"] - nr["energy"],
            "J": nr["J"], "configuration_match": config, "term_match": term,
            "parity_consistent": parity,
            "configuration_nist": nr["configuration_original"],
            "configuration_kurucz": kr["configuration_original"],
            "term_nist": nr["term_original"], "term_kurucz": kr["term_original"],
        })
    for ni, nr in n.iterrows():
        if ni in used_n:
            continue
        rows.append({
            "nist_id": int(nr["id"]), "kurucz_id": pd.NA,
            "classification": "UNMATCHED",
            "energy_nist_cm-1": nr["energy"], "energy_kurucz_cm-1": np.nan,
            "delta_energy_cm-1": np.nan, "J": nr["J"],
            "configuration_match": False, "term_match": False,
            "parity_consistent": False,
            "configuration_nist": nr["configuration_original"],
            "configuration_kurucz": "",
            "term_nist": nr["term_original"], "term_kurucz": "",
        })
    return pd.DataFrame(rows)


def compare_transitions(
    kurucz_path: Path,
    nist_path: Path,
    state_matches: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, float]]:
    """Compare transitions through the matched-state mapping, streaming Kurucz."""
    matched_states = state_matches[
        state_matches["classification"] != "UNMATCHED"
    ]
    mapping = dict(zip(matched_states["nist_id"], matched_states["kurucz_id"]))
    nist = pd.read_csv(nist_path, sep=r"\s+", header=None, names=TRANS_COLUMNS)
    nist["ku"] = nist["upper_id"].map(mapping)
    nist["kl"] = nist["lower_id"].map(mapping)
    eligible = nist[nist["ku"].notna() & nist["kl"].notna()].copy()
    eligible["key"] = [
        (min(int(a), int(b)), max(int(a), int(b)))
        for a, b in zip(eligible["ku"], eligible["kl"])
    ]
    wanted = set(eligible["key"])
    strongest: dict[tuple[int, int], tuple[float, float]] = {}
    total_kurucz = 0
    for chunk in iter_transitions(kurucz_path):
        total_kurucz += len(chunk)
        for upper, lower, aval, wn in chunk.itertuples(index=False, name=None):
            key = (min(int(upper), int(lower)), max(int(upper), int(lower)))
            if key in wanted and (
                key not in strongest or float(aval) > strongest[key][0]
            ):
                strongest[key] = (float(aval), float(wn))
    rows = []
    for row in eligible.itertuples(index=False):
        found = strongest.get(row.key)
        if found is None:
            continue
        ak, wnk = found
        rows.append({
            "nist_upper_id": row.upper_id, "nist_lower_id": row.lower_id,
            "kurucz_upper_id": row.ku, "kurucz_lower_id": row.kl,
            "wavenumber_nist_cm-1": row.wavenumber,
            "wavenumber_kurucz_cm-1": wnk,
            "delta_wavenumber_cm-1": wnk - row.wavenumber,
            "wavelength_nist_nm": 1e7 / row.wavenumber,
            "wavelength_kurucz_nm": 1e7 / wnk,
            "delta_wavelength_nm": 1e7 / wnk - 1e7 / row.wavenumber,
            "A_nist_s-1": row.A, "A_kurucz_s-1": ak,
            "delta_logA": math.log10(ak) - math.log10(row.A),
        })
    result = pd.DataFrame(rows)
    metrics = {
        "kurucz_transitions": total_kurucz,
        "nist_transitions": len(nist),
        "eligible_nist_transitions": len(eligible),
        "matched_transitions": len(result),
        "transition_match_rate": len(result) / len(nist) if len(nist) else math.nan,
    }
    return result, metrics


def compare_partition_sources(kurucz_path: Path, nist_path: Path) -> tuple[pd.DataFrame, dict[str, float]]:
    kurucz, nist = read_pf(kurucz_path), read_pf(nist_path)
    low = max(kurucz.temperature.min(), nist.temperature.min())
    high = min(kurucz.temperature.max(), nist.temperature.max())
    grid = kurucz[kurucz.temperature.between(low, high)].copy()
    grid = grid.rename(columns={"Q": "Q_kurucz"})
    grid["Q_nist"] = np.interp(grid.temperature, nist.temperature, nist.Q)
    grid["relative_difference"] = (grid.Q_kurucz - grid.Q_nist) / grid.Q_nist
    rms = float(np.sqrt(np.mean(np.square(grid.relative_difference)))) if len(grid) else math.nan
    return grid, {"pf_relative_RMS": rms, "temperature_overlap_min_K": low, "temperature_overlap_max_K": high}
