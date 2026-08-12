#!/usr/bin/env python3
"""Apply the wavenumber consistency guard to every delivered transition file.

The guard is a pipeline invariant, so the interesting question is not whether it
works on a chosen example but what it reports when applied to the dataset as a
whole. For each species this recomputes

    D = | |E_u - E_l| - wavenumber |

for every transition from the delivered .states and .trans, and reports the
distribution and the number of records exceeding the production threshold.

A species passing this check is not thereby correct -- the guard tests one identity
and nothing else -- but a species failing it cannot be, because the two quantities
are not independent.

Usage:
  python3 sweep_guards.py [--min-stage 3] [--tolerance 0.5] [--out reports/guard-sweep]
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import statistics
from pathlib import Path

STAGES = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI"]
TREES = ("Kurucz-Only-data", "Kurucz-Nist-Overlap-data")


def load_energies(path: Path) -> dict[int, float]:
    energies = {}
    with path.open() as handle:
        for line in handle:
            fields = line.split(None, 2)
            if len(fields) >= 2:
                try:
                    energies[int(fields[0])] = float(fields[1])
                except ValueError:
                    continue
    return energies


def sweep(states: Path, trans: Path, tolerance: float) -> dict:
    energies = load_energies(states)
    total = violations = unresolved = 0
    worst = 0.0
    residuals = []
    keep_every = 97  # sample the distribution; the counters stay exact
    with trans.open() as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 4:
                continue
            total += 1
            try:
                upper, lower, wavenumber = int(fields[0]), int(fields[1]), float(fields[3])
            except ValueError:
                continue
            e_upper, e_lower = energies.get(upper), energies.get(lower)
            if e_upper is None or e_lower is None:
                unresolved += 1
                continue
            residual = abs(abs(e_upper - e_lower) - wavenumber)
            if residual > worst:
                worst = residual
            if residual > tolerance:
                violations += 1
            if total % keep_every == 0:
                residuals.append(residual)
    return {
        "transitions": total,
        "unresolved_state_reference": unresolved,
        "violations": violations,
        "violation_rate": f"{violations / total:.8f}" if total else "",
        "max_residual_cm-1": f"{worst:.6g}",
        "median_residual_cm-1": f"{statistics.median(residuals):.6g}" if residuals else "",
        "p99_residual_cm-1": f"{sorted(residuals)[int(len(residuals) * 0.99)]:.6g}" if len(residuals) > 100 else "",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-stage", type=int, default=3)
    parser.add_argument("--tolerance", type=float, default=0.5)
    parser.add_argument("--out", default="reports/guard-sweep")
    args = parser.parse_args()

    targets = []
    for tree in TREES:
        for path in sorted(glob.glob(f"{tree}/*/*.trans")):
            directory = Path(path).parent
            name = directory.name
            if "-" not in name:
                continue
            stage = name.split("-", 1)[1]
            if stage not in STAGES or STAGES.index(stage) + 1 < args.min_stage:
                continue
            states = directory / Path(path).name.replace(".trans", ".states")
            if states.exists():
                targets.append((name, states, Path(path)))
    print(f"sweeping {len(targets)} species at stage {args.min_stage}+")

    rows = []
    for done, (name, states, trans) in enumerate(sorted(targets), 1):
        row = {"ion": name, "element": name.split("-")[0], "stage": name.split("-")[1]}
        row.update(sweep(states, trans, args.tolerance))
        rows.append(row)
        if row["violations"] or row["unresolved_state_reference"]:
            print(f"  {name}: {row['violations']} violations, "
                  f"{row['unresolved_state_reference']} unresolved, max {row['max_residual_cm-1']}")
        if done % 25 == 0:
            print(f"  [{done}/{len(targets)}]")

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    target = out_dir / "wavenumber-guard-sweep.csv"
    with target.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    transitions = sum(r["transitions"] for r in rows)
    violations = sum(r["violations"] for r in rows)
    unresolved = sum(r["unresolved_state_reference"] for r in rows)
    failing = [r for r in rows if r["violations"]]
    worst = max(float(r["max_residual_cm-1"]) for r in rows if r["max_residual_cm-1"])
    medians = sorted(float(r["median_residual_cm-1"]) for r in rows if r["median_residual_cm-1"])

    print(f"\nwrote {target}")
    print(f"species                     : {len(rows)}")
    print(f"transitions checked         : {transitions:,}")
    print(f"unresolved state references : {unresolved:,}")
    print(f"records above {args.tolerance} cm^-1     : {violations:,} "
          f"({violations / transitions:.3%}) in {len(failing)} species")
    print(f"largest residual anywhere   : {worst:.6g} cm^-1")
    if medians:
        print(f"per-species median residual : {medians[len(medians) // 2]:.6g} cm^-1 "
              f"(range {medians[0]:.3g} to {medians[-1]:.6g})")
    for row in sorted(failing, key=lambda r: -r["violations"])[:15]:
        print(f"  {row['ion']:<9} {row['violations']:>10,} / {row['transitions']:>12,} "
              f"max {row['max_residual_cm-1']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
