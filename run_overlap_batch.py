#!/usr/bin/env python3
"""Batch-run compare_nist.py over the non-I/II Kurucz-NIST overlap ions.

Runs one species at a time (each reads a large .trans file), logs failures
without stopping, and writes an aggregate table to
reports/nist-compare/SUMMARY_highions.csv + .md.
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OVERLAP = ROOT / "Kurucz-Nist-Overlap-data"
NIST = ROOT / "Nist-temp-data"
OUT = ROOT / "reports" / "nist-compare"


def has_all_files(dir_name: str) -> bool:
    u = dir_name.replace("-", "_")
    kur = OVERLAP / dir_name
    for ext in ("states", "trans", "pf"):
        p = kur / f"{u}__Kurucz.{ext}"
        if not (p.exists() and p.stat().st_size > 0):
            return False
        hits = [h for h in NIST.rglob(f"{u}__NIST.{ext}") if h.stat().st_size > 0]
        if not hits:
            return False
    return True


def parse_summary(text: str) -> dict:
    """Pull the headline numbers out of a per-species summary.txt."""
    def grab(pattern, default=""):
        m = re.search(pattern, text)
        return m.group(1) if m else default

    d = {
        "nist_levels": grab(r"NIST (\d+)\n?.*?== Energy", ""),
        "level_match_pct": grab(r"matched \d+/\d+ NIST levels \(([\d.]+)%\)"),
        "energy_exact_pct": grab(r"exact energy agreement \(\|dE\| <= 0\.001\): ([\d.]+)%"),
        "nist_lines": grab(r"NIST lines: (\d+);"),
        "lines_found_pct": grab(r"found in Kurucz trans: \d+ \(([\d.]+)% of all NIST lines\)"),
        "dlogA_median": grab(r"dlogA = log10\(A_Kurucz/A_NIST\): median ([+\-][\d.]+)"),
        "A_factor2_pct": grab(r"<= 0\.30 \(factor 2\): ([\d.]+)%"),
        "A_factor10_pct": grab(r"<= 1\.00 \(factor 10\): ([\d.]+)%"),
        "warning": "yes" if "WARNING" in text else "",
    }
    # Q max deviation from the "largest deviation: <ratio> at" line
    m = re.search(r"largest deviation: ([\d.]+) at", text)
    if m:
        d["Q_maxdev_pct"] = f"{abs(float(m.group(1)) - 1) * 100:.3f}"
    else:
        d["Q_maxdev_pct"] = ""
    # NIST-level count sits on the "Levels: Kurucz N, NIST M" line
    m = re.search(r"Levels: Kurucz \d+, NIST (\d+)", text)
    if m:
        d["nist_levels"] = m.group(1)
    return d


def main() -> None:
    ions = sorted(
        d.name for d in OVERLAP.iterdir()
        if d.is_dir() and not re.search(r"-(I|II)$", d.name)
    )
    todo = [i for i in ions if has_all_files(i)]
    skipped = [i for i in ions if i not in todo]
    print(f"{len(todo)} ions to run, {len(skipped)} skipped for missing files: {skipped}")

    results, failed = [], []
    for n, dir_name in enumerate(todo, 1):
        species = dir_name.replace("-", "_")
        print(f"[{n}/{len(todo)}] {species} ...", flush=True)
        log_path = OUT / f"{species}_run.log"
        proc = subprocess.run(
            [sys.executable, "compare_nist.py", "--species", species,
             "--nist-dir", "Nist-temp-data",
             "--kurucz-dir", f"Kurucz-Nist-Overlap-data/{dir_name}"],
            cwd=ROOT, capture_output=True, text=True,
        )
        log_path.write_text(proc.stdout + "\n----STDERR----\n" + proc.stderr)
        summary_file = OUT / species / "summary.txt"
        if proc.returncode != 0 or not summary_file.exists():
            failed.append(species)
            print(f"    FAILED (exit {proc.returncode}); see {log_path}")
            continue
        row = {"species": species}
        row.update(parse_summary(summary_file.read_text()))
        results.append(row)

    # ---- aggregate CSV ----
    cols = ["species", "nist_levels", "level_match_pct", "energy_exact_pct",
            "nist_lines", "lines_found_pct", "dlogA_median",
            "A_factor2_pct", "A_factor10_pct", "Q_maxdev_pct", "warning"]
    csv_path = OUT / "SUMMARY_highions.csv"
    with csv_path.open("w") as f:
        f.write(",".join(cols) + "\n")
        for r in results:
            f.write(",".join(str(r.get(c, "")) for c in cols) + "\n")

    # ---- aggregate Markdown ----
    md = [
        "# Kurucz vs NIST 对比：三价及以上离子（110 个）",
        "",
        f"运行 {len(results)} 个成功，{len(failed)} 个失败，{len(skipped)} 个因缺文件跳过。",
        "详细结果见各物种目录 `reports/nist-compare/<物种>/`。",
        "",
        "| 物种 | NIST能级 | 能级匹配% | 能量精确% | NIST谱线 | 谱线找到% | dlogA中位 | A因子2% | A因子10% | Q最大偏差% | 备注 |",
        "|------|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|-----:|:----|",
    ]
    for r in results:
        md.append(
            "| {species} | {nist_levels} | {level_match_pct} | {energy_exact_pct} "
            "| {nist_lines} | {lines_found_pct} | {dlogA_median} | {A_factor2_pct} "
            "| {A_factor10_pct} | {Q_maxdev_pct} | {warning} |".format(**{k: r.get(k, "") for k in cols})
        )
    if failed:
        md += ["", "## 失败的物种", ", ".join(failed)]
    if skipped:
        md += ["", "## 因缺文件跳过", ", ".join(skipped)]
    (OUT / "SUMMARY_highions.md").write_text("\n".join(md) + "\n")

    print(f"\nDone. {len(results)} ok, {len(failed)} failed, {len(skipped)} skipped.")
    print(f"Aggregate: {csv_path}")
    print(f"Aggregate: {OUT / 'SUMMARY_highions.md'}")


if __name__ == "__main__":
    main()
