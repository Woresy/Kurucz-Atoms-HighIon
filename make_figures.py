#!/usr/bin/env python3
"""Generate the dissertation figures that aggregate across species.

The per-species plots under reports/nist-compare/ are diagnostics; these are the
figures that carry an argument, so they are built here from the same underlying
CSVs rather than assembled by hand. Output goes to essay-overleaf/figures/ as PDF,
which is what \\includegraphics should prefer.

Usage:
  python3 make_figures.py [--out essay-overleaf/figures]
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import statistics
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

STAGES = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X"]
# Incidence of the energy-field overflow of Section 5.6, measured 2026-08-10 over
# the delivered dataset before the fix. See essay/evidence-map.md.
OVERFLOW_BY_STAGE = {"III": 0.0, "IV": 1.1, "V": 6.3, "VI": 38.8,
                     "VII": 56.4, "VIII": 60.2, "IX": 61.5, "X": 58.8}

plt.rcParams.update({
    "font.size": 9,
    "axes.grid": True,
    "grid.alpha": 0.3,
    "grid.linewidth": 0.5,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "figure.dpi": 150,
})
INK = "#1f3b57"
ACCENT = "#a8322d"


def species_dirs(min_stage: int = 3):
    for path in sorted(glob.glob("reports/nist-compare/*/")):
        name = os.path.basename(path.rstrip("/"))
        if "_" not in name:
            continue
        element, stage = name.split("_", 1)
        if stage in STAGES and STAGES.index(stage) + 1 >= min_stage:
            yield name, element, stage, path


def save(fig, out: Path, name: str) -> None:
    for suffix in ("pdf", "png"):
        fig.savefig(out / f"{name}.{suffix}", bbox_inches="tight")
    plt.close(fig)
    print(f"  {name}.pdf")


def fig_overflow_by_stage(out: Path) -> None:
    """Section 5.6: the defect's incidence against ionization stage."""
    stages = [s for s in STAGES if s in OVERFLOW_BY_STAGE]
    values = [OVERFLOW_BY_STAGE[s] for s in stages]
    fig, ax = plt.subplots(figsize=(4.4, 2.9))
    ax.plot(range(len(stages)), values, "o-", color=ACCENT, lw=1.6, ms=5)
    ax.set_xticks(range(len(stages)))
    ax.set_xticklabels(stages)
    ax.set_xlabel("Ionization stage")
    ax.set_ylabel("Lines exceeding the field width (%)")
    ax.set_ylim(-3, 70)
    ax.annotate("neutral and singly ionized\nspecies cannot reach $10^6$ cm$^{-1}$",
                xy=(0, 0), xytext=(0.6, 22), fontsize=7.5, color="0.35",
                arrowprops=dict(arrowstyle="->", color="0.5", lw=0.8))
    save(fig, out, "overflow_by_stage")


def fig_reference_levels_by_stage(out: Path) -> None:
    """Sections 6.4.3 and 6.6.3: how the reference thins with ionization stage."""
    per = {}
    for name, element, stage, path in species_dirs():
        summary = Path(path) / "summary.txt"
        if not summary.exists():
            continue
        import re
        match = re.search(r"Levels: Kurucz (\d+), NIST (\d+)", summary.read_text(errors="replace"))
        if match:
            per.setdefault(stage, []).append((int(match.group(1)), int(match.group(2))))
    stages = [s for s in STAGES if s in per]
    kurucz = [statistics.median(k for k, _ in per[s]) for s in stages]
    reference = [statistics.median(n for _, n in per[s]) for s in stages]
    counts = [len(per[s]) for s in stages]

    fig, ax = plt.subplots(figsize=(4.8, 3.0))
    x = np.arange(len(stages))
    ax.semilogy(x, kurucz, "o-", color=INK, lw=1.6, ms=5, label="This work")
    ax.semilogy(x, reference, "s-", color=ACCENT, lw=1.6, ms=5, label="Reference")
    for i, n in enumerate(counts):
        ax.annotate(f"$n$={n}", (i, reference[i]), textcoords="offset points",
                    xytext=(0, -14), ha="center", fontsize=6.5, color="0.4")
    ax.set_xticks(x)
    ax.set_xticklabels(stages)
    ax.set_xlabel("Ionization stage")
    ax.set_ylabel("Median levels per species")
    ax.legend(frameon=False, fontsize=8)
    save(fig, out, "reference_levels_by_stage")


def fig_dlogA_all_species(out: Path) -> None:
    """Section 6.5: transition-probability agreement pooled over the sample."""
    values = []
    for name, element, stage, path in species_dirs():
        matched = Path(path) / "trans_matched.csv"
        if not matched.exists():
            continue
        with matched.open() as handle:
            for row in csv.DictReader(handle):
                raw = row.get("dlogA", "")
                if raw in ("", "nan"):
                    continue
                try:
                    values.append(float(raw))
                except ValueError:
                    continue
    values = np.array(values)
    fig, ax = plt.subplots(figsize=(4.8, 3.0))
    ax.hist(np.clip(values, -2, 2), bins=80, color=INK, alpha=0.85)
    for edge, style in ((0.30, "--"), (1.00, ":")):
        for sign in (-1, 1):
            ax.axvline(sign * edge, color=ACCENT, ls=style, lw=1.0)
    inside2 = np.mean(np.abs(values) <= 0.30) * 100
    inside10 = np.mean(np.abs(values) <= 1.00) * 100
    ax.set_xlabel(r"$\Delta\log A = \log_{10}(A_\mathrm{this\,work}/A_\mathrm{ref})$")
    ax.set_ylabel("Transitions")
    ax.set_title(f"{len(values):,} matched transitions; "
                 f"{inside2:.0f}% within a factor of 2, {inside10:.0f}% within 10",
                 fontsize=8)
    save(fig, out, "dlogA_all_species")
    print(f"    pooled: n={len(values):,}, |dlogA|<=0.30: {inside2:.1f}%, <=1.00: {inside10:.1f}%")


def fig_pf_ratio_ti_ix(out: Path) -> None:
    """Section 6.6.2: the partition-function departure for Ti IX."""
    path = Path("reports/nist-compare/Ti_IX/pf_compare.csv")
    if not path.exists():
        print("  Ti_IX/pf_compare.csv missing; skipped")
        return
    temperature, ratio = [], []
    with path.open() as handle:
        for row in csv.DictReader(handle):
            try:
                temperature.append(float(row["T"]))
                ratio.append(float(row["ratio"]))
            except (ValueError, KeyError):
                continue
    fig, ax = plt.subplots(figsize=(4.6, 2.9))
    ax.semilogx(temperature, ratio, "-", color=ACCENT, lw=1.6)
    ax.axhline(1.0, color="0.5", lw=0.8, ls="--")
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel(r"$Q_\mathrm{this\,work}\,/\,Q_\mathrm{ref}$")
    ax.set_title("Ti IX", fontsize=9)
    save(fig, out, "pf_ratio_ti_ix")
    floor = min(ratio)
    print(f"    Ti IX: ratio falls to {floor:.3f} (max deviation {abs(1-floor)*100:.1f}%)")


def fig_coverage_matrix(out: Path) -> None:
    """Section 3.6: the element x stage coverage matrix."""
    source = Path("reports/ion_coverage.md")
    if not source.exists():
        print("  ion_coverage.md missing; skipped")
        return
    lines = source.read_text().splitlines()
    header = [c.strip() for c in lines[0].split("|")]
    columns = {c: i for i, c in enumerate(header) if c in STAGES}
    order = [s for s in STAGES if s in columns]
    elements, grid = [], []
    for line in lines[2:]:
        cells = [c.strip() for c in line.split("|")]
        if len(cells) <= max(columns.values()) or not cells[1].isdigit():
            continue
        elements.append(cells[2])
        grid.append([1 if cells[columns[s]] == "✓" else 0 for s in order])
    data = np.array(grid, dtype=float)
    # Stages I and II are outside this work's scope; shade them differently.
    scope = np.array([[0.45 if order[j] in ("I", "II") else 1.0
                       for j in range(len(order))] for _ in elements])
    fig, ax = plt.subplots(figsize=(5.0, 0.20 * len(elements) + 1.0))
    ax.imshow(data * scope, cmap="Blues", vmin=0, vmax=1.15, aspect="auto")
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order)
    ax.set_yticks(range(len(elements)))
    ax.set_yticklabels(elements, fontsize=6)
    ax.set_xlabel("Ionization stage")
    ax.grid(False)
    ax.set_xticks(np.arange(-0.5, len(order), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(elements), 1), minor=True)
    ax.grid(which="minor", color="white", lw=0.6)
    ax.tick_params(which="minor", length=0)
    # The two shades carry the argument of the figure, so they need naming.
    from matplotlib.patches import Patch
    cmap = plt.get_cmap("Blues")
    ax.legend(handles=[Patch(facecolor=cmap(1.0 / 1.15), label="This work (stage III+)"),
                       Patch(facecolor=cmap(0.45 / 1.15), label="Existing scope (I, II)")],
              loc="upper center", bbox_to_anchor=(0.5, -0.045), ncol=2,
              frameon=False, fontsize=7.5, handlelength=1.2, handleheight=1.0)
    save(fig, out, "coverage_matrix")
    print(f"    {int(data.sum())} ions, of which "
          f"{int(data[:, [order.index(s) for s in order if s not in ('I','II')]].sum())} at stage III+")



def fig_pf_lowering_sensitivity(out: Path) -> None:
    """Section 4.7: how much the choice of potential-lowering column matters."""
    import urllib.request
    species = [("2402", "Cr III"), ("2604", "Fe V"), ("2607", "Fe VIII"),
               ("2808", "Ni IX"), ("3008", "Zn IX")]
    fig, ax = plt.subplots(figsize=(4.8, 3.0))
    colours = plt.cm.viridis(np.linspace(0.05, 0.85, len(species)))
    for (code, name), colour in zip(species, colours):
        path = Path(f"/tmp/kz/pf{code}.dat")
        if not path.exists():
            try:
                urllib.request.urlretrieve(
                    f"http://kurucz.harvard.edu/atoms/{code}/partfn{code}.dat", path)
            except Exception:
                continue
        rows = []
        for line in path.read_text(errors="replace").splitlines():
            parts = line.split()
            if len(parts) >= 10:
                try:
                    rows.append([float(x) for x in parts[2:10]])
                except ValueError:
                    continue
        table = np.array(rows)
        keep = table[:, 1] > 0
        # Column 6 is the most severe lowering that still carries a value; the
        # seventh is tabulated as zero throughout and is not a physical result.
        ax.semilogx(table[keep, 0], table[keep, 6] / table[keep, 1],
                    lw=1.6, color=colour, label=name)
    ax.axhline(1.0, color="0.5", lw=0.8, ls="--")
    ax.set_xlabel("Temperature (K)")
    ax.set_ylabel(r"$Q(-16000)\,/\,Q(-500)$")
    ax.set_ylim(0, 1.08)
    ax.legend(frameon=False, fontsize=7.5, loc="lower left")
    save(fig, out, "pf_lowering_sensitivity")



def fig_validation_axes(out: Path) -> None:
    """Section 6.1: what each validation axis is capable of rejecting.

    Organised by capability rather than by procedure, because the chapter's argument
    is that a validation result is a statement about the property tested. The final
    row is the point: a defect none of the three axes examines is not caught by
    running them more often.
    """
    axes = ["Axis 1\ninternal", "Axis 2\nregression", "Axis 3\nexternal"]
    defects = [
        "Malformed\noutput",
        "Endpoints vs\nwavenumber",
        "Incomplete\nsource",
        "Method alters\nvalues",
        "Error in\nthe source",
        "Property not\nexamined",
    ]
    # 2 = rejects it, 1 = partially, 0 = structurally cannot
    grid = np.array([
        [2, 2, 2, 0, 0, 0],
        [0, 0, 0, 2, 0, 0],
        [1, 1, 1, 0, 2, 0],
    ])
    fig, ax = plt.subplots(figsize=(6.6, 2.4))
    ax.imshow(grid, cmap="Blues", vmin=0, vmax=2.6, aspect="auto")
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            mark = {2: "\u2713", 1: "\u25cb", 0: ""}[grid[i, j]]
            ax.text(j, i, mark, ha="center", va="center", fontsize=13,
                    color="white" if grid[i, j] == 2 else "0.25")
    ax.set_xticks(range(len(defects)))
    ax.set_xticklabels(defects, fontsize=7.5)
    ax.set_yticks(range(len(axes)))
    ax.set_yticklabels(axes, fontsize=8)
    ax.set_xticks(np.arange(-0.5, len(defects), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(axes), 1), minor=True)
    ax.grid(which="minor", color="white", lw=1.4)
    ax.grid(False)
    ax.tick_params(which="minor", length=0)
    ax.tick_params(axis="x", length=0, pad=4)
    ax.annotate("caught by neither:\nfound by separate audit", xy=(5, 1), xytext=(5, -1.15),
                ha="center", va="center", fontsize=7, color=ACCENT, style="italic",
                arrowprops=dict(arrowstyle="->", color=ACCENT, lw=0.9))
    ax.set_ylim(2.5, -1.7)
    save(fig, out, "validation_axes")



def fig_isoelectronic_spectra(out: Path) -> None:
    """Section 7.4: LTE spectra of a 19-electron sequence over the full range.

    Computed directly from the delivered files rather than by the external consumer,
    so that the wavelength range is not fixed in advance -- the point of the figure is
    where the strong lines fall, which a pre-chosen window cannot show.
    """
    import glob, math, os
    C2, T = 1.4388, 20000.0
    fig, axes = plt.subplots(4, 1, figsize=(6.2, 5.4), sharex=True)
    ions = ["Ti-IV", "V-V", "Cr-VI", "Mn-VII"]
    colours = plt.cm.viridis(np.linspace(0.05, 0.8, len(ions)))
    for ax, ion, colour in zip(axes, ions, colours):
        found = [x for x in glob.glob("*-data/" + ion) if os.path.isdir(x)]
        if not found:
            continue
        directory = found[0]
        energy, weight = {}, {}
        for line in open(glob.glob(directory + "/*.states")[0]):
            parts = line.split()
            if len(parts) >= 4:
                energy[int(parts[0])] = float(parts[1])
                weight[int(parts[0])] = int(parts[2])
        grid = [(float(p[0]), float(p[1]))
                for p in (l.split() for l in open(glob.glob(directory + "/*.pf")[0])) if len(p) >= 2]
        partition = min(grid, key=lambda x: abs(x[0] - T))[1]
        lam, intensity = [], []
        for line in open(glob.glob(directory + "/*.trans")[0]):
            parts = line.split()
            if len(parts) < 4:
                continue
            upper, lower, a_value, wavenumber = int(parts[0]), int(parts[1]), float(parts[2]), float(parts[3])
            if wavenumber <= 0 or upper not in weight or lower not in energy:
                continue
            value = (a_value * weight[upper] / (8 * math.pi * 2.99792458e10 * wavenumber ** 2)
                     * math.exp(-C2 * energy[lower] / T) / partition
                     * (1 - math.exp(-C2 * wavenumber / T)))
            lam.append(1e7 / wavenumber)
            intensity.append(value)
        lam, intensity = np.array(lam), np.array(intensity)
        peak = intensity.max()
        ax.vlines(lam, 0, intensity / peak, lw=0.6, color=colour)
        ax.axvspan(275, 330, color="0.85", zorder=0)
        ax.set_xscale("log")
        ax.set_xlim(8, 400)
        ax.set_ylim(0, 1.15)
        ax.set_yticks([])
        ax.text(0.012, 0.72, ion.replace("-", " "), transform=ax.transAxes, fontsize=8.5)
        strongest = lam[int(np.argmax(intensity))]
        ax.annotate(f"{strongest:.1f} nm", xy=(strongest, 1.0), xytext=(strongest * 1.5, 0.85),
                    fontsize=7, color=ACCENT,
                    arrowprops=dict(arrowstyle="-", color=ACCENT, lw=0.7))
    axes[0].text(300, 1.28, "window of the\nearlier figure",
                 ha="center", fontsize=6.5, color="0.4")
    axes[-1].set_xlabel("Vacuum wavelength (nm)")
    fig.supylabel("Relative intensity at 20 000 K", fontsize=8)
    fig.subplots_adjust(hspace=0.12)
    save(fig, out, "isoelectronic_spectra")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="essay-overleaf/figures")
    args = parser.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    print(f"writing to {out}/")
    fig_coverage_matrix(out)
    fig_overflow_by_stage(out)
    fig_reference_levels_by_stage(out)
    fig_dlogA_all_species(out)
    fig_pf_ratio_ti_ix(out)
    fig_pf_lowering_sensitivity(out)
    fig_validation_axes(out)
    fig_isoelectronic_spectra(out)
    fig_sc_iii_worked_example(out)
    return 0



def fig_sc_iii_worked_example(out: Path) -> None:
    """The two panels of the Section 6 worked example, both read back from the delivered files.

    Left: the stick spectrum the external consumer produces from them. Right: radiative
    lifetimes it recomputes from the transition list, against those the source supplies.
    The right panel is the one that carries an argument -- agreement there is evidence the
    transition probabilities are internally consistent with the lifetimes they imply.
    """
    sticks = list(csv.DictReader(open("reports/Sc-III/stick_spectrum.csv")))
    wavenumber = np.array([float(row["wavenumber_cm-1"]) for row in sticks])
    intensity = np.array([float(row["intensity"]) for row in sticks])
    intensity = intensity / intensity.max()

    lifetimes = [row for row in csv.DictReader(open("reports/Sc-III/lifetime_comparison.csv"))
                 if row["relative_error"]]
    supplied = np.array([float(row["lifetime"]) for row in lifetimes])
    recomputed = np.array([float(row["calculated_lifetime_s"]) for row in lifetimes])
    error = np.abs(np.array([float(row["relative_error"]) for row in lifetimes]))

    fig, (left, right) = plt.subplots(1, 2, figsize=(7.4, 3.0))

    floor = 1e-8
    visible = intensity > floor
    left.vlines(wavenumber[visible] / 1e3, floor, intensity[visible],
                color=INK, linewidth=0.5, alpha=0.85)
    left.set_yscale("log")
    left.set_ylim(floor, 3.0)
    left.set_xlabel(r"Wavenumber ($10^3$ cm$^{-1}$)")
    left.set_ylabel("Relative intensity at 1000 K")
    left.set_title(f"(a) Stick spectrum, {visible.sum():,} of {len(sticks):,} lines shown",
                   fontsize=8, loc="left")

    span = [0.8 * min(supplied.min(), recomputed.min()), 1.25 * max(supplied.max(), recomputed.max())]
    right.plot(span, span, color="0.6", linewidth=0.8, zorder=1)
    right.scatter(supplied, recomputed, s=9, color=ACCENT, alpha=0.7,
                  edgecolors="none", zorder=2)
    right.set_xscale("log")
    right.set_yscale("log")
    right.set_xlim(span)
    right.set_ylim(span)
    right.set_xlabel("Lifetime supplied by the source (s)")
    right.set_ylabel("Lifetime recomputed (s)")
    right.set_title(f"(b) Radiative lifetimes, {len(lifetimes)} levels", fontsize=8, loc="left")
    right.text(0.04, 0.94,
               f"median |error| {statistics.median(error) * 100:.2f}%\nmax {error.max() * 100:.2f}%",
               transform=right.transAxes, fontsize=7.5, va="top", color=INK)

    fig.tight_layout()
    save(fig, out, "sc_iii_worked_example")

if __name__ == "__main__":
    raise SystemExit(main())
