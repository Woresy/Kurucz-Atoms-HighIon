#!/usr/bin/env python3
"""Reorder existing ExoMol .trans files so the upper state comes first.

process_kurucz_atom.py used to write transition endpoints in raw Kurucz-file
order, which is upper-first in only part of the rows; the ExoMol format
requires (upper, lower, A, wavenumber). This walks <root>/*/exomol/*.trans,
and on every row where E[id1] < E[id2] swaps the two 12-character ID fields
by byte slicing, leaving the A and wavenumber columns untouched. Files are
rewritten atomically and left alone when already fully ordered. A sibling
.zip archive, if present, is rebuilt from the current exomol files.

Usage:
  python3 fix_trans_order.py                 # fix Kurucz-data
  python3 fix_trans_order.py Kurucz-data test-out
  python3 fix_trans_order.py --only Fe-I     # single species
"""

from __future__ import annotations

import argparse
import time
import zipfile
from pathlib import Path

import numpy as np

# .trans row layout: [0:12] id1, [13:25] id2, [26:] A + wavenumber
ID1 = slice(0, 12)
ID2 = slice(13, 25)


def load_energies(states_path: Path) -> np.ndarray:
    ids, energies = [], []
    with states_path.open() as handle:
        for line in handle:
            tokens = line.split(None, 2)
            if len(tokens) < 2:
                continue
            ids.append(int(tokens[0]))
            energies.append(float(tokens[1]))
    e_by_id = np.full(max(ids) + 1, np.nan)
    e_by_id[ids] = energies
    return e_by_id


def fix_trans_file(trans_path: Path, e_by_id: np.ndarray) -> tuple[int, int]:
    """Rewrite trans_path with upper state first; return (total, swapped) rows."""
    temporary_path = trans_path.with_name(f".{trans_path.name}.fixorder.tmp")
    total = swapped = 0
    try:
        with trans_path.open() as src, temporary_path.open("w") as dst:
            for line in src:
                total += 1
                id1 = int(line[ID1])
                id2 = int(line[ID2])
                if e_by_id[id1] < e_by_id[id2]:
                    swapped += 1
                    line = line[ID2] + " " + line[ID1] + " " + line[26:]
                dst.write(line)
        if swapped:
            temporary_path.replace(trans_path)
        else:
            temporary_path.unlink()
    except BaseException:
        temporary_path.unlink(missing_ok=True)
        raise
    return total, swapped


def refresh_zip(exomol_dir: Path) -> Path | None:
    archives = list(exomol_dir.glob("*.zip"))
    if not archives:
        return None
    archive = archives[0]
    members = sorted(p for p in exomol_dir.glob("*") if p.suffix in {".states", ".trans", ".pf"})
    with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as bundle:
        for member in members:
            bundle.write(member, member.name)
    return archive


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("roots", nargs="*", default=["Kurucz-data"], type=Path)
    parser.add_argument("--only", help="process just this species directory name, e.g. Fe-I")
    parser.add_argument("--no-zip", action="store_true", help="skip rebuilding sibling .zip archives")
    args = parser.parse_args()

    species_dirs = []
    for root in args.roots:
        species_dirs += sorted(
            d for d in Path(root).iterdir()
            if d.is_dir() and (d / "exomol").is_dir() and (args.only is None or d.name == args.only)
        )

    grand_total = grand_swapped = files_changed = 0
    start = time.time()
    for species_dir in species_dirs:
        exomol_dir = species_dir / "exomol"
        for trans_path in sorted(exomol_dir.glob("*.trans")):
            states_candidates = list(exomol_dir.glob("*.states"))
            if not states_candidates:
                print(f"{species_dir.name}: no .states file, skipped")
                continue
            e_by_id = load_energies(states_candidates[0])
            total, swapped = fix_trans_file(trans_path, e_by_id)
            grand_total += total
            grand_swapped += swapped
            note = ""
            if swapped:
                files_changed += 1
                if not args.no_zip:
                    archive = refresh_zip(exomol_dir)
                    if archive:
                        note = f", rebuilt {archive.name}"
            print(f"{species_dir.name}: {total} rows, swapped {swapped} ({swapped / max(total, 1) * 100:.1f}%){note}", flush=True)

    elapsed = time.time() - start
    print(f"\ndone: {len(species_dirs)} species, {grand_total} rows scanned, "
          f"{grand_swapped} swapped, {files_changed} files rewritten, {elapsed:.0f}s")


if __name__ == "__main__":
    main()
