#!/usr/bin/env python3
"""Test whether Einstein A derived from .lines reproduces the A published in .agafgf.

Some ions publish no full .agafgf, so their transitions are rebuilt from the .lines
file, which carries log(gf) rather than log(A). That is a change of provenance, and
it has to be shown to be a change of provenance only: the derived coefficient must
agree with the published one wherever both exist.

    A_ul = 6.6702e15 * gf / (g_u * lambda_A^2) = 0.66702 * gf * wn^2 / g_u

with g_u = 2J_u + 1 taken from whichever endpoint has the larger |E|. This script
pairs the two sources on their endpoint energies and J values, then reports the
distribution of |log10(A_lines / A_agafgf)|.

The comparison is meaningful only within the rounding of the source fields: Kurucz
stores log(gf) to three decimals, so agreement cannot be shown to better than about
1e-3 dex regardless of how correct the conversion is.

Usage:
  python3 check_a_equivalence.py [--ions Ni-X Mn-III] [--out reports/a-equivalence]
"""

from __future__ import annotations

import argparse
import csv
import math
import statistics
import subprocess
from pathlib import Path

from process_kurucz_atom import (
    BASE_URL,
    ELEMENTS,
    discover_ion,
    is_gz,
    iter_agafgf_records,
    parse_lines_transition,
    roman,
)

STAGES = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI"]


def fetch(code: str, name: str, timeout: int) -> str:
    """Download one Kurucz source file, transparently decompressing .gz members."""
    url = f"{BASE_URL}/{code}/{name}"
    raw = subprocess.run(["curl", "-sfL", url, "-m", str(timeout)], capture_output=True).stdout
    if is_gz(name):
        import gzip
        try:
            raw = gzip.decompress(raw)
        except (OSError, EOFError):
            return ""
    return raw.decode("latin-1", "replace")


def stream(code: str, name: str, timeout: int):
    """Yield one Kurucz source file line by line without holding it in memory."""
    url = f"{BASE_URL}/{code}/{name}"
    command = ["curl", "-sfL", url, "-m", str(timeout)]
    if is_gz(name):
        process = subprocess.Popen(f"{' '.join(command)} | gzip -dc", shell=True,
                                   stdout=subprocess.PIPE, text=True, errors="replace")
    else:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, text=True, errors="replace")
    try:
        yield from process.stdout
    finally:
        process.stdout.close()
        process.wait()


def key(e1: float, j1: float, e2: float, j2: float) -> tuple:
    """Endpoint-pair identity, orientation-independent and tolerant of |E| sign."""
    a = (round(abs(e1), 3), round(j1, 1))
    b = (round(abs(e2), 3), round(j2, 1))
    return tuple(sorted([a, b]))


def compare(ion: str, timeout: int) -> tuple[dict, list[dict]]:
    element, stage = ion.split("-", 1)
    charge = STAGES.index(stage)
    discovery = discover_ion(element, charge, timeout, BASE_URL)

    published: dict[tuple, float] = {}
    for source in discovery.trans_sources:
        for e1, j1, _l1, e2, j2, _l2, _wn, a_value in iter_agafgf_records(fetch(discovery.code, source, timeout)):
            if a_value > 0:
                published[key(e1, j1, e2, j2)] = a_value

    # The .lines member can run to gigabytes while the published set is a few tens of
    # thousands of rows, so stream it and keep only the pairs actually under test.
    derived: dict[tuple, float] = {}
    for source in discovery.lines_sources:
        for line in stream(discovery.code, source, timeout):
            record = parse_lines_transition(line)
            if record is None:
                continue
            e1, j1, _l1, e2, j2, _l2, _wn, a_value = record
            if a_value <= 0:
                continue
            pair = key(e1, j1, e2, j2)
            if pair in published:
                derived[pair] = a_value

    rows = []
    for pair, a_published in published.items():
        a_derived = derived.get(pair)
        if a_derived is None:
            continue
        rows.append({
            "ion": ion,
            "E_lower": pair[0][0], "J_lower": pair[0][1],
            "E_upper": pair[1][0], "J_upper": pair[1][1],
            "A_agafgf_s-1": f"{a_published:.6e}",
            "A_from_lines_s-1": f"{a_derived:.6e}",
            "abs_delta_log10": f"{abs(math.log10(a_derived / a_published)):.3e}",
        })

    deltas = [float(r["abs_delta_log10"]) for r in rows]
    summary = {
        "ion": ion,
        "code": discovery.code,
        "agafgf_sources": ";".join(discovery.trans_sources),
        "lines_sources": ";".join(discovery.lines_sources),
        "published_transitions": len(published),
        "lines_transitions": len(derived),
        "matched": len(rows),
        "match_rate": f"{len(rows) / len(published):.4f}" if published else "",
        "median_abs_delta_log10": f"{statistics.median(deltas):.3e}" if deltas else "",
        "max_abs_delta_log10": f"{max(deltas):.3e}" if deltas else "",
        "within_0.01_dex": f"{sum(d <= 0.01 for d in deltas) / len(deltas):.4f}" if deltas else "",
        "within_0.001_dex": f"{sum(d <= 0.001 for d in deltas) / len(deltas):.4f}" if deltas else "",
    }
    return summary, rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ions", nargs="+", default=["Ni-X", "Mn-III"])
    parser.add_argument("--out", default="reports/a-equivalence")
    parser.add_argument("--timeout", type=int, default=600)
    args = parser.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    summaries, details = [], []
    for ion in args.ions:
        element = ion.split("-", 1)[0]
        if element not in ELEMENTS:
            print(f"{ion}: unknown element")
            continue
        print(f"comparing {ion} ...")
        summary, rows = compare(ion, args.timeout)
        summaries.append(summary)
        details.extend(rows)
        print(f"  agafgf={summary['published_transitions']}  lines={summary['lines_transitions']}  "
              f"matched={summary['matched']}  median|dlog10|={summary['median_abs_delta_log10']}  "
              f"max={summary['max_abs_delta_log10']}")

    if summaries:
        with (out_dir / "summary.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(summaries[0].keys()))
            writer.writeheader()
            writer.writerows(summaries)
    if details:
        with (out_dir / "matched-transitions.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(details[0].keys()))
            writer.writeheader()
            writer.writerows(details)
    print(f"\nwrote {out_dir}/summary.csv and {out_dir}/matched-transitions.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
