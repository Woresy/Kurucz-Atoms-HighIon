"""Publication-ready plots for generic Kurucz/NIST validation outputs."""

from __future__ import annotations

import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .reporting import TRANSITION_CLASSES, log_intensity_bounds


LOGGER = logging.getLogger(__name__)

COLORS = {"Kurucz": "#2166ac", "NIST": "#d73027"}


def _save(figure: plt.Figure, output: Path, stem: str) -> list[str]:
    png = output / f"{stem}.png"
    pdf = output / f"{stem}.pdf"
    figure.savefig(png, dpi=220, bbox_inches="tight")
    figure.savefig(pdf, bbox_inches="tight")
    plt.close(figure)
    return [png.name, pdf.name]


def plot_state_comparison(output: Path, ion: str) -> list[str]:
    frame = pd.read_csv(output / "kurucz_vs_nist_states.csv")
    matched = frame[frame["classification"] != "UNMATCHED"].copy()
    if matched.empty:
        raise ValueError(f"no matched Kurucz/NIST states for {ion}")
    figure, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].scatter(
        matched["energy_nist_cm-1"], matched["energy_kurucz_cm-1"],
        s=15, alpha=0.65, color=COLORS["Kurucz"],
    )
    limits = [
        min(matched["energy_nist_cm-1"].min(), matched["energy_kurucz_cm-1"].min()),
        max(matched["energy_nist_cm-1"].max(), matched["energy_kurucz_cm-1"].max()),
    ]
    axes[0].plot(limits, limits, "k--", linewidth=1)
    axes[0].set(
        xlabel=r"NIST energy (cm$^{-1}$)",
        ylabel=r"Kurucz energy (cm$^{-1}$)",
        title=f"{ion}: matched energy levels",
    )
    axes[1].hist(
        matched["delta_energy_cm-1"].dropna(),
        bins=min(60, max(10, int(np.sqrt(len(matched))))),
        color=COLORS["Kurucz"], alpha=0.82,
    )
    axes[1].axvline(0, color="black", linestyle="--", linewidth=1)
    axes[1].set(
        xlabel=r"$E_{\rm Kurucz}-E_{\rm NIST}$ (cm$^{-1}$)",
        ylabel="Matched states",
        title="Energy residual distribution",
    )
    for axis in axes:
        axis.grid(alpha=0.22)
    figure.tight_layout()
    return _save(figure, output, "state_energy_comparison")


def plot_transition_comparison(output: Path, ion: str) -> list[str]:
    frame = pd.read_csv(output / "kurucz_vs_nist_transitions.csv")
    valid = frame[
        (frame["A_nist_s-1"] > 0) & (frame["A_kurucz_s-1"] > 0)
    ]
    if valid.empty:
        raise ValueError(f"no transitions with positive A on both sides for {ion}")
    figure, axes = plt.subplots(1, 2, figsize=(12, 5))
    x = np.log10(valid["A_nist_s-1"])
    y = np.log10(valid["A_kurucz_s-1"])
    axes[0].scatter(x, y, s=11, alpha=0.4, color=COLORS["Kurucz"])
    limits = [min(x.min(), y.min()), max(x.max(), y.max())]
    axes[0].plot(limits, limits, "k--", linewidth=1)
    axes[0].set(
        xlabel=r"log$_{10}$ A(NIST) (s$^{-1}$)",
        ylabel=r"log$_{10}$ A(Kurucz) (s$^{-1}$)",
        title=f"{ion}: Einstein A comparison",
    )
    axes[1].hist(
        valid["delta_logA"].clip(-2, 2),
        bins=min(70, max(10, int(np.sqrt(len(valid))))),
        color=COLORS["NIST"], alpha=0.82,
    )
    axes[1].axvline(0, color="black", linestyle="--", linewidth=1)
    axes[1].set(
        xlabel=r"log$_{10}$(A$_{\rm Kurucz}$/A$_{\rm NIST}$), clipped to ±2",
        ylabel="Matched transitions",
        title="A-value residual distribution",
    )
    for axis in axes:
        axis.grid(alpha=0.22)
    figure.tight_layout()
    return _save(figure, output, "transition_A_comparison")


def plot_partition_comparison(output: Path, ion: str) -> list[str]:
    frame = pd.read_csv(output / "kurucz_vs_nist_pf.csv")
    figure, axes = plt.subplots(
        2, 1, figsize=(8, 7), sharex=True,
        gridspec_kw={"height_ratios": [2, 1]},
    )
    axes[0].plot(frame["temperature"], frame["Q_kurucz"],
                 label="Kurucz", color=COLORS["Kurucz"], linewidth=2)
    axes[0].plot(frame["temperature"], frame["Q_nist"],
                 label="NIST", color=COLORS["NIST"], linewidth=2, linestyle="--")
    axes[0].set(ylabel="Partition function Q(T)",
                title=f"{ion}: Kurucz vs NIST partition functions")
    axes[0].legend()
    axes[1].plot(frame["temperature"], frame["relative_difference"],
                 color="#4d4d4d")
    axes[1].axhline(0, color="black", linestyle="--", linewidth=1)
    axes[1].set(
        xlabel="Temperature (K)",
        ylabel="(QK − QN) / QN",
    )
    for axis in axes:
        axis.grid(alpha=0.22)
    figure.tight_layout()
    return _save(figure, output, "partition_function_cross_source")


def _chunk_temperatures(chunk: pd.DataFrame) -> pd.Series:
    """Temperature per row, from the numeric column when the writer supplied it."""
    if "temperature_K" in chunk.columns:
        return pd.to_numeric(chunk["temperature_K"], errors="coerce")
    return pd.to_numeric(
        chunk["temperature_source"].astype(str).str.extract(r"T([0-9.]+)K")[0],
        errors="coerce",
    )


def _read_temperature_spectra(
    path: Path, top_n: int = 8_000
) -> tuple[dict[float, pd.DataFrame], dict[float, int]]:
    """Stream a possibly multi-gigabyte spectrum and retain strong plot lines.

    Returns the retained lines per temperature plus the full line count per
    temperature, so a figure can state how much of the spectrum it is showing.
    """
    strongest: dict[float, pd.DataFrame] = {}
    totals: dict[float, int] = {}
    for chunk in pd.read_csv(path, chunksize=400_000):
        temperatures = _chunk_temperatures(chunk)
        for temperature in temperatures.dropna().unique():
            selected = chunk[np.isclose(temperatures, temperature)]
            key = float(temperature)
            totals[key] = totals.get(key, 0) + len(selected)
            previous = strongest.get(key)
            if previous is not None:
                selected = pd.concat([previous, selected], ignore_index=True)
            strongest[key] = selected.nlargest(min(top_n, len(selected)), "intensity")
    if not strongest:
        raise ValueError(f"no temperature-tagged spectrum rows in {path}")
    return strongest, totals


def _decade_ticks(
    floor_exponent: float, ceiling_exponent: float, max_ticks: int = 8
) -> tuple[list[float], list[str]]:
    """Tick offsets (in decades above the floor) labelled with real intensities."""
    first = int(np.ceil(floor_exponent))
    last = int(np.floor(ceiling_exponent))
    step = max(1, int(np.ceil((last - first + 1) / max_ticks)))
    exponents = list(range(first, last + 1, step))
    return (
        [exponent - floor_exponent for exponent in exponents],
        [rf"$10^{{{exponent}}}$" for exponent in exponents],
    )


def plot_spectrum_comparison(output: Path, ion: str) -> list[str]:
    """Compare both stick spectra in absolute intensity on a shared log scale.

    Neither source is normalised: an order-of-magnitude offset between Kurucz
    and NIST is a real result of the comparison and has to stay visible.
    """
    kurucz_by_temperature, kurucz_totals = _read_temperature_spectra(
        output / "pyexocross_kurucz_spectrum.csv"
    )
    nist_by_temperature, nist_totals = _read_temperature_spectra(
        output / "pyexocross_nist_spectrum.csv"
    )
    temperatures = sorted(
        set(kurucz_by_temperature) & set(nist_by_temperature)
    )
    if not temperatures:
        raise ValueError(f"no common Kurucz/NIST spectrum temperatures for {ion}")

    outputs: list[str] = []
    for temperature in temperatures:
        temperature_label = f"{temperature:g} K"
        temperature_token = f"{temperature:g}K"
        frames = {
            "Kurucz": kurucz_by_temperature[temperature],
            "NIST": nist_by_temperature[temperature],
        }
        totals = {
            "Kurucz": kurucz_totals[temperature],
            "NIST": nist_totals[temperature],
        }
        # Anchoring the floor to the weakest source *and class* guarantees each
        # stays on the axis; the stronger ones then sit higher by their true
        # margin. Forbidden lines are usually far weaker than allowed ones, and
        # a floor set from the allowed lines alone would drop them off the plot.
        floor, ceiling = log_intensity_bounds([
            float(subset["intensity"].max())
            for frame in frames.values()
            for _, subset in frame.groupby("transition_type", sort=True)
            if not subset.empty
        ])
        floor_exponent = float(np.log10(floor))
        ceiling_exponent = float(np.log10(ceiling))
        span = ceiling_exponent - floor_exponent
        tick_offsets, tick_labels = _decade_ticks(floor_exponent, ceiling_exponent)

        # Mirrored view: bar length is the number of decades above the shared
        # floor, so the two halves use one common absolute intensity scale.
        figure, axis = plt.subplots(figsize=(13, 6.5))
        for source, frame in frames.items():
            visible = frame[frame["intensity"] >= floor]
            direction = 1.0 if source == "Kurucz" else -1.0
            for name_class, (class_label, class_color, width) in (
                TRANSITION_CLASSES.items()
            ):
                subset = visible[visible["transition_type"] == name_class]
                if subset.empty:
                    continue
                axis.vlines(
                    subset["vacuum_wavelength_nm"], 0,
                    direction * (np.log10(subset["intensity"]) - floor_exponent),
                    color=class_color if class_color else COLORS[source],
                    alpha=0.62 if class_color is None else 0.95,
                    linewidth=0.8 * width,
                    label=f"{source} — {class_label}: {len(subset)}",
                )
            axis.plot([], [], " ", label=f"    ({totals[source]} lines in range)")
        axis.axhline(0, color="black", linewidth=0.8)
        axis.set_ylim(-span * 1.05, span * 1.05)
        # A tick exactly on the mirror axis must not be emitted twice.
        mirrored = {
            -offset: label for offset, label in zip(tick_offsets, tick_labels)
        }
        mirrored.update(zip(tick_offsets, tick_labels))
        axis.set_yticks(sorted(mirrored))
        axis.set_yticklabels([mirrored[offset] for offset in sorted(mirrored)])
        axis.set(
            xlabel="Vacuum wavelength (nm)",
            ylabel=(
                "Absolute LTE line intensity, log scale\n"
                "(Kurucz ↑ / NIST ↓)"
            ),
            title=(
                f"{ion}: mirrored PyExoCross spectra at {temperature_label} — "
                f"shared absolute scale, floor {floor:.2e}"
            ),
        )
        axis.grid(alpha=0.22)
        axis.legend(loc="upper right", fontsize="small")
        figure.tight_layout()
        outputs += _save(
            figure, output,
            f"pyexocross_spectrum_{temperature_token}_absolute_comparison",
        )

        # Stacked view: same absolute limits on both panels, so panel-to-panel
        # heights remain directly comparable.
        figure, axes = plt.subplots(2, 1, figsize=(13, 8), sharex=True, sharey=True)
        for axis, (source, frame) in zip(axes, frames.items()):
            visible = frame[frame["intensity"] >= floor]
            for name_class, (class_label, class_color, width) in (
                TRANSITION_CLASSES.items()
            ):
                subset = visible[visible["transition_type"] == name_class]
                if subset.empty:
                    continue
                axis.vlines(
                    subset["vacuum_wavelength_nm"], floor, subset["intensity"],
                    color=class_color if class_color else COLORS[source],
                    alpha=0.72 if class_color is None else 0.95,
                    linewidth=0.75 * width,
                    label=f"{class_label} ({len(subset)})",
                )
            axis.legend(loc="upper left", fontsize="x-small")
            axis.set_yscale("log")
            axis.set_ylim(floor, ceiling * 10.0)
            axis.set_ylabel(f"{source}\nabsolute intensity")
            axis.grid(alpha=0.22)
            axis.text(
                0.995, 0.93,
                f"{len(visible)} lines above floor (of {totals[source]} in range)",
                transform=axis.transAxes, ha="right", va="top", fontsize="small",
            )
        axes[0].set_title(
            f"{ion}: PyExoCross spectra at {temperature_label} "
            "(shared absolute-intensity axis)"
        )
        axes[-1].set_xlabel("Vacuum wavelength (nm)")
        figure.tight_layout()
        outputs += _save(
            figure, output,
            f"pyexocross_spectrum_{temperature_token}_stacked_absolute",
        )
    return outputs


def _has_data_rows(path: Path) -> bool:
    """True when the CSV exists and carries at least one row below its header."""
    if not path.exists():
        return False
    with path.open(encoding="utf-8", errors="ignore") as handle:
        return sum(1 for _ in zip(range(2), handle)) >= 2


def generate_cross_source_plots(output: Path, ion: str) -> list[str]:
    """Generate every plot whose source CSVs are available and non-empty.

    A plot that cannot be produced is skipped with a warning: one unusable
    comparison must not cost the caller the figures that did work.
    """
    plotters = [
        (("kurucz_vs_nist_states.csv",), plot_state_comparison),
        (("kurucz_vs_nist_transitions.csv",), plot_transition_comparison),
        (("kurucz_vs_nist_pf.csv",), plot_partition_comparison),
        (
            ("pyexocross_kurucz_spectrum.csv", "pyexocross_nist_spectrum.csv"),
            plot_spectrum_comparison,
        ),
    ]
    generated: list[str] = []
    for required, plotter in plotters:
        if not all(_has_data_rows(output / name) for name in required):
            continue
        try:
            generated.extend(plotter(output, ion))
        except Exception as exc:
            LOGGER.warning(
                "%s: skipping %s (%s: %s)",
                ion, plotter.__name__, type(exc).__name__, exc,
            )
    return generated
