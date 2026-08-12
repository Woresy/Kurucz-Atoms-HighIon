#!/usr/bin/env python3
"""Compare Kurucz-derived ExoMol files against NIST-derived ExoMol files.

Compares three quantities for one species:
  1. Energy levels (.states)  -- match by statistical weight g + nearest energy,
                                 one-to-one greedy assignment within --etol.
  2. Transitions  (.trans)    -- map NIST (upper, lower) level pairs onto matched
                                 Kurucz level IDs and compare Einstein A values.
  3. Partition functions (.pf) -- interpolate NIST Q(T) onto the Kurucz T grid
                                 over the overlapping range and compare.

Outputs per-species CSVs, plots and a summary.txt under <out>/<species>/.

Usage:
  python3 compare_nist.py --species Fe_I --nist-dir Nist-temp-data
  # kurucz dir defaults to Kurucz-data/<Fe-I>/exomol, output to reports/nist-compare
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

KURUCZ_STATES_COLS = ["id", "E", "g", "J", "unc", "lifetime", "lande", "config", "term", "tag"]
NIST_STATES_COLS = ["id", "E", "g", "J", "unc", "lande", "config", "term", "parity"]
TRANS_COLS = ["id1", "id2", "A", "wn"]
PF_COLS = ["T", "Q"]


def find_file(root: Path, filename: str) -> Path:
    """Locate filename directly under root or in any subdirectory."""
    direct = root / filename
    if direct.exists():
        return direct
    hits = sorted(root.rglob(filename))
    if not hits:
        raise FileNotFoundError(f"{filename} not found under {root}")
    return hits[0]


def read_states(path: Path, columns: list[str]) -> pd.DataFrame:
    """Whitespace-parse a states file whose config field may itself contain spaces.

    Leading columns (up to config) are numeric and trailing columns (after
    config) are single tokens, so everything left in the middle is the config.
    """
    n_lead = columns.index("config")
    n_trail = len(columns) - n_lead - 1
    rows = []
    with path.open() as handle:
        for line in handle:
            tokens = line.split()
            if not tokens:
                continue
            rows.append(tokens[:n_lead] + [" ".join(tokens[n_lead:-n_trail])] + tokens[-n_trail:])
    df = pd.DataFrame(rows, columns=columns)
    for col in columns[:n_lead]:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df["id"] = df["id"].astype(int)
    df["g"] = df["g"].astype(int)
    return df


def read_nist_states(path: Path) -> pd.DataFrame:
    """Read a NIST states file, tolerating the variant without a lande column.

    Some NIST exports lack lande entirely; detect by counting leading numeric
    tokens on the first data line (6 with lande, 5 without).
    """
    with path.open() as handle:
        for line in handle:
            tokens = line.split()
            if tokens:
                break
    n_numeric = 0
    for tok in tokens:
        try:
            float(tok)
        except ValueError:
            break
        n_numeric += 1
    columns = NIST_STATES_COLS if n_numeric >= 6 else [c for c in NIST_STATES_COLS if c != "lande"]
    df = read_states(path, columns)
    if "lande" not in df.columns:
        df["lande"] = np.nan
    return df


def fix_inconsistent_g(states: pd.DataFrame, name: str, log: list[str]) -> pd.DataFrame:
    """If the g column widely disagrees with 2J+1, rebuild it from J.

    Matching keys on g, so a corrupt g column (seen in some NIST exports)
    would silently kill the level match; J is the more trustworthy column.
    """
    bad = (states["g"] - (2 * states["J"] + 1)).abs() > 0.01
    if bad.mean() > 0.05:
        log.append(f"WARNING: {name} g column != 2J+1 for {int(bad.sum())}/{len(states)} "
                   f"levels; using g = 2J+1 derived from J for matching")
        states = states.copy()
        states["g"] = (2 * states["J"] + 1).round().astype(int)
    return states


def read_trans(path: Path) -> pd.DataFrame:
    return pd.read_csv(
        path, sep=r"\s+", names=TRANS_COLS,
        dtype={"id1": np.int64, "id2": np.int64, "A": np.float64, "wn": np.float64},
    )


def read_pf(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, sep=r"\s+", names=PF_COLS, dtype=np.float64)


def orient_trans(trans: pd.DataFrame, states: pd.DataFrame, name: str, log: list[str]) -> pd.DataFrame:
    """Return trans with columns upper/lower ordered so E[upper] > E[lower].

    The Kurucz pipeline writes endpoints in raw-file order, which is not
    guaranteed to be upper-first, so orientation is decided per row by energy.
    """
    e_by_id = np.full(states["id"].max() + 1, np.nan)
    e_by_id[states["id"].to_numpy()] = states["E"].to_numpy()
    e1 = e_by_id[trans["id1"].to_numpy()]
    e2 = e_by_id[trans["id2"].to_numpy()]
    first_is_upper = e1 >= e2
    out = trans.copy()
    out["upper"] = np.where(first_is_upper, trans["id1"], trans["id2"])
    out["lower"] = np.where(first_is_upper, trans["id2"], trans["id1"])
    ediff = np.abs(e1 - e2)
    dev = np.abs(ediff - trans["wn"].to_numpy())
    log.append(
        f"{name}: {len(out)} transitions; first-column-is-upper in "
        f"{first_is_upper.mean() * 100:.1f}% of rows; "
        f"median |(E_u-E_l) - wn| = {np.nanmedian(dev):.4g} cm-1, max = {np.nanmax(dev):.4g} cm-1"
    )
    return out[["upper", "lower", "A", "wn"]]


def match_levels(nist: pd.DataFrame, kurucz: pd.DataFrame, etol: float):
    """Greedy one-to-one matching: same g, smallest |dE| first, within etol."""
    candidates = []
    kur_by_g = {g: sub for g, sub in kurucz.groupby("g")}
    for row in nist.itertuples():
        sub = kur_by_g.get(row.g)
        if sub is None:
            continue
        de = (sub["E"] - row.E).abs()
        # keep the 3 nearest candidates so a NIST level can fall back if its
        # nearest Kurucz level is claimed by a closer NIST level
        for kidx in de.nsmallest(3).index:
            candidates.append((de[kidx], row.Index, kidx))
    candidates.sort()
    used_nist, used_kur, pairs = set(), set(), []
    nearest_de = {}
    for de, nidx, kidx in candidates:
        nearest_de.setdefault(nidx, de)
        if de > etol or nidx in used_nist or kidx in used_kur:
            continue
        used_nist.add(nidx)
        used_kur.add(kidx)
        pairs.append((nidx, kidx))
    return pairs, nearest_de


def config_consistent(nist_config: str, kurucz_config: str) -> bool:
    n = str(nist_config).replace(".", "").lower()
    k = str(kurucz_config).replace(".", "").lower()
    return n.endswith(k) or k.endswith(n)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species", required=True, help="e.g. Fe_I")
    parser.add_argument("--nist-dir", required=True, type=Path)
    parser.add_argument("--kurucz-dir", type=Path, default=None,
                        help="default: Kurucz-data/<species with - >/exomol")
    parser.add_argument("--out", type=Path, default=Path("reports/nist-compare"))
    parser.add_argument("--etol", type=float, default=2.0,
                        help="level-match energy tolerance in cm-1 (default 2.0)")
    args = parser.parse_args()

    species = args.species
    kurucz_dir = args.kurucz_dir or Path("Kurucz-data") / species.replace("_", "-") / "exomol"
    out_dir = args.out / species
    out_dir.mkdir(parents=True, exist_ok=True)
    log: list[str] = [f"Species: {species}", f"Kurucz dir: {kurucz_dir}", f"NIST dir: {args.nist_dir}",
                      f"Energy tolerance: {args.etol} cm-1", ""]

    kur_states = read_states(kurucz_dir / f"{species}__Kurucz.states", KURUCZ_STATES_COLS)
    nist_states = read_nist_states(find_file(args.nist_dir, f"{species}__NIST.states"))
    log.append(f"Levels: Kurucz {len(kur_states)}, NIST {len(nist_states)}")
    kur_states = fix_inconsistent_g(kur_states, "Kurucz", log)
    nist_states = fix_inconsistent_g(nist_states, "NIST", log)

    # ---------------- 1. energy levels ----------------
    pairs, nearest_de = match_levels(nist_states, kur_states, args.etol)
    matched = pd.DataFrame({
        "nist_id": [nist_states.at[n, "id"] for n, _ in pairs],
        "kurucz_id": [kur_states.at[k, "id"] for _, k in pairs],
        "g": [nist_states.at[n, "g"] for n, _ in pairs],
        "E_nist": [nist_states.at[n, "E"] for n, _ in pairs],
        "E_kurucz": [kur_states.at[k, "E"] for _, k in pairs],
        "term_nist": [nist_states.at[n, "term"] for n, _ in pairs],
        "term_kurucz": [kur_states.at[k, "term"] for _, k in pairs],
        "config_nist": [nist_states.at[n, "config"] for n, _ in pairs],
        "config_kurucz": [kur_states.at[k, "config"] for _, k in pairs],
        "lande_nist": [nist_states.at[n, "lande"] for n, _ in pairs],
        "lande_kurucz": [kur_states.at[k, "lande"] for _, k in pairs],
    }).sort_values("E_nist")
    matched["dE"] = matched["E_kurucz"] - matched["E_nist"]
    matched["term_match"] = (
        matched["term_nist"].str.strip().str.lower()
        == matched["term_kurucz"].str.strip().str.lower()
    )
    matched["config_match"] = [
        config_consistent(a, b) for a, b in zip(matched["config_nist"], matched["config_kurucz"])
    ]
    matched.to_csv(out_dir / "levels_matched.csv", index=False)

    matched_nidx = {n for n, _ in pairs}
    unmatched = nist_states.loc[[i for i in nist_states.index if i not in matched_nidx]].copy()
    unmatched["nearest_dE_same_g"] = [nearest_de.get(i, np.nan) for i in unmatched.index]
    unmatched.to_csv(out_dir / "levels_nist_unmatched.csv", index=False)

    abs_de = matched["dE"].abs()
    lande_diff = (matched["lande_kurucz"] - matched["lande_nist"]).abs()
    log += [
        "",
        "== Energy levels ==",
        f"matched {len(matched)}/{len(nist_states)} NIST levels "
        f"({len(matched) / len(nist_states) * 100:.1f}%) within {args.etol} cm-1",
        f"|dE|: median {abs_de.median():.4g}, mean {abs_de.mean():.4g}, "
        f"rms {np.sqrt((matched['dE'] ** 2).mean()):.4g}, max {abs_de.max():.4g} cm-1",
        f"exact energy agreement (|dE| <= 0.001): {(abs_de <= 0.001).mean() * 100:.1f}%",
        f"term label agreement: {matched['term_match'].mean() * 100:.1f}%",
        f"config consistency (suffix match): {matched['config_match'].mean() * 100:.1f}%",
        f"Lande g |diff| median: {lande_diff.median():.4g} (where both defined)",
        f"unmatched NIST levels: {len(unmatched)} (see levels_nist_unmatched.csv)",
    ]

    # ---------------- 2. transitions ----------------
    kur_trans = orient_trans(read_trans(kurucz_dir / f"{species}__Kurucz.trans"), kur_states, "Kurucz trans", log)
    nist_trans = orient_trans(read_trans(find_file(args.nist_dir, f"{species}__NIST.trans")), nist_states, "NIST trans", log)

    dup = kur_trans.duplicated(subset=["upper", "lower"]).sum()
    if dup:
        log.append(f"note: Kurucz trans has {dup} duplicate (upper,lower) pairs; keeping strongest A")
        kur_trans = kur_trans.sort_values("A", ascending=False).drop_duplicates(["upper", "lower"])

    nist2kur = dict(zip(matched["nist_id"], matched["kurucz_id"]))
    nt = nist_trans.copy()
    nt["kur_upper"] = nt["upper"].map(nist2kur)
    nt["kur_lower"] = nt["lower"].map(nist2kur)
    endpoint_ok = nt["kur_upper"].notna() & nt["kur_lower"].notna()

    merged = nt[endpoint_ok].merge(
        kur_trans.rename(columns={"upper": "kur_upper", "lower": "kur_lower",
                                  "A": "A_kurucz", "wn": "wn_kurucz"}),
        on=["kur_upper", "kur_lower"], how="left",
    ).rename(columns={"upper": "nist_upper", "lower": "nist_lower",
                      "A": "A_nist", "wn": "wn_nist"})
    found = merged["A_kurucz"].notna()
    merged["dlogA"] = np.log10(merged["A_kurucz"] / merged["A_nist"])
    merged["dwn"] = merged["wn_kurucz"] - merged["wn_nist"]
    merged.to_csv(out_dir / "trans_matched.csv", index=False)

    miss = pd.concat([
        nt[~endpoint_ok].assign(reason="endpoint level unmatched"),
        merged[~found][["nist_upper", "nist_lower", "A_nist", "wn_nist"]]
        .rename(columns={"nist_upper": "upper", "nist_lower": "lower",
                         "A_nist": "A", "wn_nist": "wn"})
        .assign(reason="no Kurucz transition for matched levels"),
    ])
    miss.to_csv(out_dir / "trans_nist_unmatched.csv", index=False)

    ok = merged[found]
    log += [
        "",
        "== Transitions (Einstein A) ==",
        f"NIST lines: {len(nist_trans)}; endpoints matched to Kurucz levels: {int(endpoint_ok.sum())}; "
        f"found in Kurucz trans: {len(ok)} ({len(ok) / len(nist_trans) * 100:.1f}% of all NIST lines)",
    ]
    if len(ok):
        log += [
            f"dlogA = log10(A_Kurucz/A_NIST): median {ok['dlogA'].median():+.3f}, "
            f"mean {ok['dlogA'].mean():+.3f}, std {ok['dlogA'].std():.3f}",
            f"|dlogA| <= 0.05 (12%): {(ok['dlogA'].abs() <= 0.05).mean() * 100:.1f}%   "
            f"<= 0.30 (factor 2): {(ok['dlogA'].abs() <= 0.30).mean() * 100:.1f}%   "
            f"<= 1.00 (factor 10): {(ok['dlogA'].abs() <= 1.00).mean() * 100:.1f}%",
            f"wavenumber |diff| median: {ok['dwn'].abs().median():.4g} cm-1",
        ]
    log.append(f"NIST lines without Kurucz counterpart: {len(miss)} (see trans_nist_unmatched.csv)")

    # ---------------- 3. partition function ----------------
    kur_pf = read_pf(kurucz_dir / f"{species}__Kurucz.pf")
    nist_pf = read_pf(find_file(args.nist_dir, f"{species}__NIST.pf"))
    t_lo = max(kur_pf["T"].min(), nist_pf["T"].min())
    t_hi = min(kur_pf["T"].max(), nist_pf["T"].max())
    grid = kur_pf[(kur_pf["T"] >= t_lo) & (kur_pf["T"] <= t_hi)].copy()
    grid["Q_nist"] = np.interp(grid["T"], nist_pf["T"], nist_pf["Q"])
    grid = grid.rename(columns={"Q": "Q_kurucz"})
    grid["ratio"] = grid["Q_kurucz"] / grid["Q_nist"]
    grid.to_csv(out_dir / "pf_compare.csv", index=False)

    log += ["", "== Partition function ==",
            f"overlap: {t_lo:.0f}-{t_hi:.0f} K ({len(grid)} Kurucz grid points)"]
    for t in (200, 500, 1000, 2000, 3000, 4000, 5000, 6000):
        if t_lo <= t <= t_hi:
            r = np.interp(t, grid["T"], grid["ratio"])
            log.append(f"  Q_Kurucz/Q_NIST at {t:>5d} K: {r:.4f}")
    worst = grid.loc[(grid["ratio"] - 1).abs().idxmax()]
    log.append(f"  largest deviation: {worst['ratio']:.4f} at {worst['T']:.0f} K")

    # ---------------- plots ----------------
    fig, ax = plt.subplots(figsize=(6, 4))
    clipped = matched["dE"].clip(-args.etol, args.etol)
    ax.hist(clipped, bins=80, color="#4878a8")
    ax.set_xlabel(r"$E_{Kurucz} - E_{NIST}$ (cm$^{-1}$)")
    ax.set_ylabel("levels")
    ax.set_title(f"{species}: level energy differences ({len(matched)} matched)")
    fig.tight_layout()
    fig.savefig(out_dir / "levels_dE_hist.png", dpi=150)
    plt.close(fig)

    if len(ok):
        fig, ax = plt.subplots(figsize=(5.5, 5.5))
        ax.scatter(np.log10(ok["A_nist"]), np.log10(ok["A_kurucz"]), s=8, alpha=0.4, color="#4878a8")
        lims = [min(np.log10(ok["A_nist"]).min(), np.log10(ok["A_kurucz"]).min()) - 0.5,
                max(np.log10(ok["A_nist"]).max(), np.log10(ok["A_kurucz"]).max()) + 0.5]
        ax.plot(lims, lims, "k--", lw=1)
        ax.set_xlim(lims), ax.set_ylim(lims)
        ax.set_xlabel(r"log$_{10}$ A$_{NIST}$ (s$^{-1}$)")
        ax.set_ylabel(r"log$_{10}$ A$_{Kurucz}$ (s$^{-1}$)")
        ax.set_title(f"{species}: Einstein A comparison ({len(ok)} lines)")
        fig.tight_layout()
        fig.savefig(out_dir / "trans_logA_scatter.png", dpi=150)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(6, 4))
        ax.hist(ok["dlogA"].clip(-2, 2), bins=80, color="#4878a8")
        ax.axvline(0, color="k", lw=1, ls="--")
        ax.set_xlabel(r"log$_{10}$(A$_{Kurucz}$ / A$_{NIST}$)")
        ax.set_ylabel("lines")
        ax.set_title(f"{species}: A-value ratio distribution")
        fig.tight_layout()
        fig.savefig(out_dir / "trans_dlogA_hist.png", dpi=150)
        plt.close(fig)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 6), sharex=True,
                                   gridspec_kw={"height_ratios": [2, 1]})
    ax1.plot(kur_pf["T"], kur_pf["Q"], label="Kurucz", color="#4878a8")
    ax1.plot(nist_pf["T"], nist_pf["Q"], label="NIST levels", color="#c44e52", ls="--")
    ax1.set_ylabel("Q(T)")
    ax1.set_yscale("log")
    ax1.legend()
    ax1.set_title(f"{species}: partition function")
    ax2.plot(grid["T"], grid["ratio"], color="#4878a8")
    ax2.axhline(1, color="k", lw=1, ls="--")
    ax2.set_xlabel("T (K)")
    ax2.set_ylabel(r"Q$_{Kurucz}$ / Q$_{NIST}$")
    fig.tight_layout()
    fig.savefig(out_dir / "pf_compare.png", dpi=150)
    plt.close(fig)

    summary = "\n".join(log)
    (out_dir / "summary.txt").write_text(summary + "\n")
    print(summary)
    print(f"\nOutputs written to {out_dir}/")


if __name__ == "__main__":
    main()
