#!/usr/bin/env python3
"""Tabulate declared versus delivered transition counts for every processed ion.

Every Kurucz level file states, in its first record, how many transitions its own
computation produced. That figure is an authoritative completeness reference which
costs nothing to obtain and needs no external data set: an ion whose delivered
.trans carries far fewer rows than its own source declares is reading a subset
product, not losing rows to filtering.

The pipeline already applies this check at run time (see
``check_transition_completeness``), but only for ions it is actively processing.
This script reconstructs the same comparison for an existing dataset without
rebuilding it: the declared count comes from a range request for the first bytes of
each .gam, and the delivered count from the local .trans.

The ratio here is a lower bound on the pipeline's own, which compares records *read*
(written plus dropped); rows dropped as unmappable are invisible from the output
side. For the completeness argument that is the more relevant quantity anyway --
what a consumer receives, not what the reader consumed.

Usage:
  python3 build_completeness_table.py [--min-stage 3] [--out reports/kurucz-incomplete]
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import gzip
import re
import subprocess
from pathlib import Path

from process_kurucz_atom import (
    BASE_URL,
    ELEMENTS,
    NUMBER_TO_ELEMENT,
    TRANSITION_COMPLETENESS_THRESHOLD,
    declared_line_count,
    discover_ion,
    roman,
)

TREES = ("Kurucz-Only-data", "Kurucz-Nist-Overlap-data")
HEADER_BYTES = 300


def fetch_header(code: str, gam: str, timeout: int) -> str:
    """Return the first record of a .gam without downloading the whole file."""
    url = f"{BASE_URL}/{code}/{gam}"
    if gam.endswith(("-gz", ".gz")):
        raw = subprocess.run(
            ["curl", "-sfL", url, "-m", str(timeout * 3)], capture_output=True
        ).stdout
        try:
            return gzip.decompress(raw).decode("latin-1", "replace")
        except (OSError, EOFError):
            return ""
    return subprocess.run(
        ["curl", "-sfL", "-r", f"0-{HEADER_BYTES}", url, "-m", str(timeout)],
        capture_output=True, text=True, errors="replace",
    ).stdout


def local_trans(ion: str) -> tuple[int, str]:
    """Line count and location of the delivered .trans for one ion."""
    for tree in TREES:
        path = Path(tree) / ion / f"{ion.replace('-', '_', 1)}__Kurucz.trans"
        if path.exists():
            with path.open("rb") as handle:
                return sum(chunk.count(b"\n") for chunk in iter(lambda: handle.read(1 << 22), b"")), str(path)
    return 0, ""


def survey(element: str, charge: int, timeout: int) -> dict | None:
    ion = f"{element}-{roman(charge + 1)}"
    delivered, path = local_trans(ion)
    if not path:
        return None
    discovery = discover_ion(element, charge, timeout, BASE_URL)
    declared: dict[str, int] = {}
    for group in discovery.groups:
        count = declared_line_count(fetch_header(discovery.code, group.gam, timeout))
        if count is not None:
            declared[group.gam] = count
    total = sum(declared.values())
    return {
        "ion": ion,
        "element": element,
        "stage": roman(charge + 1),
        "charge": charge,
        "code": discovery.code,
        "gam_files": ";".join(sorted(declared)),
        "declared_transitions": total,
        "delivered_transitions": delivered,
        "completeness": f"{delivered / total:.6f}" if total else "",
        "complete": "yes" if total and delivered / total >= TRANSITION_COMPLETENESS_THRESHOLD else "no",
        "trans_sources": ";".join(discovery.trans_sources),
        "path": path,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-stage", type=int, default=3, help="Lowest spectroscopic stage to include.")
    parser.add_argument("--out", default="reports/kurucz-incomplete")
    parser.add_argument("--timeout", type=int, default=45)
    parser.add_argument("--workers", type=int, default=10)
    args = parser.parse_args()

    targets = []
    for tree in TREES:
        for path in sorted(Path(tree).glob("*/")):
            name = path.name
            if "-" not in name:
                continue
            element, stage = name.split("-", 1)
            if element not in ELEMENTS:
                continue
            try:
                charge = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI"].index(stage)
            except ValueError:
                continue
            if charge + 1 >= args.min_stage:
                targets.append((element, charge))
    targets = sorted(set(targets), key=lambda t: (ELEMENTS[t[0]], t[1]))
    print(f"surveying {len(targets)} ions at stage {args.min_stage}+")

    rows = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(survey, e, c, args.timeout): (e, c) for e, c in targets}
        for done, future in enumerate(concurrent.futures.as_completed(futures), 1):
            try:
                row = future.result()
            except Exception as exc:  # noqa: BLE001 - one bad ion must not sink the survey
                element, charge = futures[future]
                print(f"  {element}-{roman(charge + 1)}: {type(exc).__name__}: {exc}")
                continue
            if row:
                rows.append(row)
            if done % 25 == 0:
                print(f"  [{done}/{len(targets)}]")

    rows.sort(key=lambda r: float(r["completeness"]) if r["completeness"] else -1.0)
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    target = out_dir / "transition-completeness.csv"
    with target.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    incomplete = [r for r in rows if r["complete"] == "no" and r["completeness"]]
    missing = sum(r["declared_transitions"] - r["delivered_transitions"] for r in incomplete)
    print(f"\nwrote {target} ({len(rows)} ions)")
    print(f"below the {TRANSITION_COMPLETENESS_THRESHOLD:.0%} threshold: {len(incomplete)} ions")
    print(f"transitions declared but not delivered by those ions: {missing:,}\n")
    for row in incomplete:
        print(f"  {row['ion']:<9} {float(row['completeness']):>9.4%}  "
              f"{row['delivered_transitions']:>12,} / {row['declared_transitions']:>12,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
