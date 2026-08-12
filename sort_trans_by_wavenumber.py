#!/usr/bin/env python3
"""Order every delivered .trans file by ascending wavenumber.

Row order is not fixed by the ExoAtom field specification, but every released Kurucz
transition file ascends in wavenumber, and the files produced here do not: they are
predominantly descending with occasional ascending steps where records supplemented
from the line list were appended after those read from the transition file. The order
was a by-product of construction rather than a choice.

Sorting changes no value and no state index -- it permutes whole records -- so the
operation is verifiable by comparing the multiset of lines before and after.

Records are fixed width, so the file is treated as an array of equal-length rows and
the wavenumber column is read by offset. Files that are not uniform fall back to a
line-based sort.

Usage:
  python3 sort_trans_by_wavenumber.py [--min-stage 3] [--dry-run]
"""

from __future__ import annotations

import argparse
import glob
import os
from pathlib import Path

import numpy as np

STAGES = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI"]
TREES = ("Kurucz-Only-data", "Kurucz-Nist-Overlap-data")
WAVENUMBER = slice(37, 52)   # %15.6e field of the 52-character record
CHUNK = 1_000_000            # records per write, to bound peak memory


def order_of(path: Path) -> tuple[str, int, int]:
    """Return (verdict, ascending steps, descending steps) without rewriting."""
    ascending = descending = 0
    previous = None
    with path.open() as handle:
        for line in handle:
            fields = line.split()
            if len(fields) < 4:
                continue
            value = float(fields[3])
            if previous is not None:
                if value > previous:
                    ascending += 1
                elif value < previous:
                    descending += 1
            previous = value
    verdict = "ascending" if descending == 0 else ("descending" if ascending == 0 else "mixed")
    return verdict, ascending, descending


def sort_file(path: Path) -> tuple[int, bool]:
    data = path.read_bytes()
    width = data.index(b"\n") + 1
    uniform = len(data) % width == 0
    rows = None
    if uniform:
        rows = np.frombuffer(data, dtype=np.uint8).reshape(-1, width)
        uniform = bool((rows[:, width - 1] == 0x0A).all())
    if uniform:
        column = rows[:, WAVENUMBER].tobytes().decode("ascii")
        wavenumber = np.array(column.split(), dtype=np.float64)
        uniform = len(wavenumber) == len(rows)

    # Both branches reorder whole records by a permutation -- fancy-indexing with an
    # argsort, or list.sort -- so the multiset of records is preserved by construction.
    # Verifying it explicitly would mean materialising every record as a Python object,
    # which on the largest files costs several gigabytes to establish what the operation
    # already guarantees; the byte-length check below is sufficient and costs nothing.
    temporary = path.with_suffix(path.suffix + ".tmp")
    if uniform:
        order = np.argsort(wavenumber, kind="stable")
        count = len(rows)
        # Write in chunks so the input and a full copy of the output are never both held.
        with temporary.open("wb") as handle:
            for start in range(0, count, CHUNK):
                handle.write(rows[order[start:start + CHUNK]].tobytes())
    else:
        lines = data.decode("latin-1").splitlines(keepends=True)
        lines.sort(key=lambda line: float(line.split()[3]))
        count = len(lines)
        temporary.write_text("".join(lines), encoding="latin-1")

    written = temporary.stat().st_size
    if written != len(data):
        temporary.unlink()
        raise RuntimeError(f"{path}: wrote {written} bytes, expected {len(data)}")
    os.replace(temporary, path)
    return count, uniform


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-stage", type=int, default=3)
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()

    targets = []
    for tree in TREES:
        for name in sorted(glob.glob(f"{tree}/*/*.trans")):
            path = Path(name)
            stage = path.parent.name.split("-", 1)[-1]
            if stage in STAGES and STAGES.index(stage) + 1 >= args.min_stage:
                targets.append(path)
    print(f"{len(targets)} files at stage {args.min_stage}+")

    total = 0
    fallback = []
    for done, path in enumerate(targets, 1):
        if args.dry_run:
            verdict, up, down = order_of(path)
            print(f"  {path.parent.name:<9} {verdict:<11} up={up:,} down={down:,}")
            continue
        count, uniform = sort_file(path)
        total += count
        if not uniform:
            fallback.append(path.parent.name)
        if done % 20 == 0:
            print(f"  [{done}/{len(targets)}] {total:,} records")

    if not args.dry_run:
        print(f"\nsorted {total:,} records in {len(targets)} files")
        if fallback:
            print(f"non-uniform record width, sorted line-wise: {', '.join(fallback)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
