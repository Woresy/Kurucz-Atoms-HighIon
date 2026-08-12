#!/usr/bin/env python3
"""Build flat ExoMol-style folders for ions available only in Kurucz.

The ``Comparison`` sheet in ``reports/kurucz_vs_nist.xlsx`` marks these ions
with ``K``.  By default stages I and II are excluded, so processing starts at
stage III (charge 2).  Final files are written as::

    Kurucz-Only-data/Cr-III/Cr_III__Kurucz.{states,trans,pf}

The underlying processor briefly creates an ``exomol`` subdirectory; this
driver moves its files up one level after each successful run.  Existing
non-empty final files make a species resumable without downloading it again.
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xlsx", type=Path, default=Path("reports/kurucz_vs_nist.xlsx"))
    parser.add_argument("--processor", type=Path, default=Path("process_kurucz_atom.py"))
    parser.add_argument("--output-root", type=Path, default=Path("Kurucz-Only-data"))
    parser.add_argument("--min-charge", type=int, default=2,
                        help="Minimum charge (default 2, excluding stages I and II).")
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--start-at", help="Resume at an ion such as Cu-V (inclusive).")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def selected_ions(workbook: Path, min_charge: int = 2) -> list[tuple[str, str, int]]:
    table = pd.read_excel(workbook, sheet_name="Comparison")
    metadata = {"Z", "Element", "# ions"}
    stages = [column for column in table.columns if column not in metadata]
    ions: list[tuple[str, str, int]] = []
    for _, row in table.iterrows():
        element = str(row["Element"]).strip()
        if not element.isalpha() or element in {"D", "T"}:
            continue
        for charge, stage in enumerate(stages):
            if charge >= min_charge and str(row[stage]).strip().upper() == "K":
                ions.append((element, str(stage), charge))
    return ions


def flatten_output(output_root: Path, name: str) -> None:
    ion_dir = output_root / name
    exomol_dir = ion_dir / "exomol"
    if not exomol_dir.is_dir():
        return
    for source in exomol_dir.iterdir():
        if source.is_file():
            target = ion_dir / source.name
            if target.exists():
                target.unlink()
            shutil.move(str(source), str(target))
    exomol_dir.rmdir()


def complete(output_root: Path, element: str, stage: str) -> bool:
    """Return whether the stable outputs needed for resume are present.

    A transition file is not mandatory: a few published Kurucz ion directories
    (currently Pd-IV in this selection) contain levels and a partition function
    but no transition source at all.
    """
    stem = f"{element}_{stage}__Kurucz"
    ion_dir = output_root / f"{element}-{stage}"
    states = ion_dir / f"{stem}.states"
    pf = ion_dir / f"{stem}.pf"
    return states.is_file() and states.stat().st_size > 0 and pf.is_file() and pf.stat().st_size > 0


def main() -> int:
    args = parse_args()
    if not args.xlsx.is_file():
        raise SystemExit(f"Workbook not found: {args.xlsx}")
    if not args.processor.is_file():
        raise SystemExit(f"Processor not found: {args.processor}")

    ions = selected_ions(args.xlsx, args.min_charge)
    if args.start_at:
        names = [f"{element}-{stage}" for element, stage, _ in ions]
        try:
            ions = ions[names.index(args.start_at):]
        except ValueError:
            raise SystemExit(f"--start-at ion is not selected: {args.start_at}") from None

    print(f"Selected {len(ions)} Kurucz-only ions (charge >= {args.min_charge})", flush=True)
    if args.dry_run:
        for index, (element, stage, charge) in enumerate(ions, 1):
            print(f"{index:3d}\t{element}-{stage}\tcharge={charge}")
        return 0

    args.output_root.mkdir(parents=True, exist_ok=True)
    failures: list[str] = []
    skipped = 0
    for index, (element, stage, charge) in enumerate(ions, 1):
        name = f"{element}-{stage}"
        if not args.overwrite and complete(args.output_root, element, stage):
            skipped += 1
            print(f"[{index}/{len(ions)}] {name}: already complete", flush=True)
            continue
        command = [
            sys.executable, str(args.processor),
            "--element", element,
            "--charge", str(charge),
            "--data-root", str(args.output_root),
            "--timeout", str(args.timeout),
            "--no-save-raw",
            "--no-save-intermediate",
        ]
        if args.overwrite:
            command.append("--overwrite")
        print(f"\n=== Kurucz only [{index}/{len(ions)}] {name} ===", flush=True)
        result = subprocess.run(command, stdin=subprocess.DEVNULL)
        flatten_output(args.output_root, name)
        if result.returncode or not complete(args.output_root, element, stage):
            failures.append(name)
            print(f"FAILED or incomplete: {name}", flush=True)

    print(f"\nCompleted {len(ions) - len(failures)}/{len(ions)} ({skipped} already present)")
    if failures:
        print("Failed/incomplete ions: " + ", ".join(failures))
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
