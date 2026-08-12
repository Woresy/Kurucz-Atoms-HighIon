#!/usr/bin/env python3
"""Run process_kurucz_atom.py for ions present in both Kurucz and NIST.

The overlap workbook's ``Comparison`` sheet uses ``B`` for an ion present in
both databases.  Roman-numeral column I maps to charge 0, II to charge 1, etc.
"""

import argparse
import subprocess
import sys
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Process only ions marked B (both databases) in kurucz_vs_nist.xlsx."
    )
    parser.add_argument(
        "--xlsx",
        default="reports/kurucz_vs_nist.xlsx",
        help="Kurucz/NIST comparison workbook.",
    )
    parser.add_argument(
        "--processor",
        default="process_kurucz_atom.py",
        help="Kurucz processing script.",
    )
    parser.add_argument("--data-root", default="Kurucz-data")
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--states-only", action="store_true")
    parser.add_argument(
        "--start-at",
        help="Resume at this ion, for example Ti-II (the named ion is included).",
    )
    parser.add_argument("--dry-run", action="store_true", help="Print the selection without processing.")
    return parser.parse_args()


def overlap_ions(workbook: Path) -> list[tuple[str, str, int]]:
    table = pd.read_excel(workbook, sheet_name="Comparison")
    metadata = {"Z", "Element", "# ions"}
    stage_columns = [column for column in table.columns if column not in metadata]
    ions: list[tuple[str, str, int]] = []
    for _, row in table.iterrows():
        element = str(row["Element"]).strip()
        # D, T, 3He and any other isotope-only rows cannot be passed to the Kurucz
        # element-symbol processor.  They are never marked B, but keep the guard
        # explicit so a future workbook cannot launch an invalid request.
        if not element.isalpha() or element in {"D", "T"}:
            continue
        for charge, stage in enumerate(stage_columns):
            if str(row[stage]).strip().upper() == "B":
                ions.append((element, str(stage), charge))
    return ions


def main() -> int:
    args = parse_args()
    workbook = Path(args.xlsx)
    processor = Path(args.processor)
    if not workbook.is_file():
        raise SystemExit(f"Workbook not found: {workbook}")
    if not processor.is_file():
        raise SystemExit(f"Processor not found: {processor}")

    ions = overlap_ions(workbook)
    if args.start_at:
        names = [f"{element}-{stage}" for element, stage, _ in ions]
        try:
            ions = ions[names.index(args.start_at):]
        except ValueError:
            raise SystemExit(f"--start-at ion is not in the overlap: {args.start_at}") from None

    print(f"Selected {len(ions)} Kurucz/NIST overlap ions", flush=True)
    if args.dry_run:
        for index, (element, stage, charge) in enumerate(ions, 1):
            print(f"{index:3d}\t{element}-{stage}\tcharge={charge}")
        return 0

    failures: list[str] = []
    for index, (element, stage, charge) in enumerate(ions, 1):
        name = f"{element}-{stage}"
        command = [
            sys.executable,
            str(processor),
            "--element", element,
            "--charge", str(charge),
            "--data-root", args.data_root,
            "--timeout", str(args.timeout),
        ]
        if args.overwrite:
            command.append("--overwrite")
        if args.states_only:
            command.append("--states-only")
        print(f"\n=== overlap [{index}/{len(ions)}] {name} ===", flush=True)
        result = subprocess.run(command, stdin=subprocess.DEVNULL)
        if result.returncode:
            failures.append(name)
            print(f"FAILED {name}: exit code {result.returncode}", flush=True)

    print(f"\nCompleted: {len(ions) - len(failures)}/{len(ions)}", flush=True)
    if failures:
        print("Failed ions: " + ", ".join(failures), flush=True)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
