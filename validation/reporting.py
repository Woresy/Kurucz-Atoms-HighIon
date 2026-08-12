"""Plots and Markdown/JSON reporting for a validation run."""

from __future__ import annotations

import json
from pathlib import Path
import re
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Stick intensities span well over a hundred decades, so every spectrum axis is
# logarithmic and bounded from below. The floor is placed this many decades
# under the *weakest* series' strongest line, which keeps every series visible
# while leaving their absolute separation on the axis instead of hiding it.
SPECTRUM_FLOOR_DECADES = 6.0

TEMPERATURE_PATTERN = re.compile(r"T([0-9.]+)K")

# PyExoCross stick columns for an LTE run. Only the first two are guaranteed;
# the state descriptors are used when present and verified before being trusted.
STICK_COLUMNS = [
    "wavenumber_cm-1", "intensity",
    "J_upper", "energy_upper_cm-1", "J_lower", "energy_lower_cm-1",
]

# Colour, width and label per transition class. "E1_candidate" is not a claim
# that a line is allowed: it means the states file carried no parity, so only
# the angular-momentum half of the selection rule could be tested.
# Forbidden lines are few and are the whole point of the classification, so
# they get a hue far from both source colours plus extra width; a colour merely
# adjacent to the source red was not separable at stick widths.
TRANSITION_CLASSES = {
    "E1_allowed": ("E1 allowed", None, 1.0),
    "E1_candidate": ("E1 (ΔJ only, no parity in states)", None, 1.0),
    "forbidden": ("forbidden", "#00A03C", 2.4),
}


def read_state_parities(path: Path) -> pd.DataFrame | None:
    """Energy/J to parity, or None when the states file carries no parity column.

    Kurucz states files record no parity at all, so this is the honest signal
    that only the angular-momentum half of the E1 rule can be evaluated.
    """
    rows: list[tuple[float, float, str]] = []
    with path.open(encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            tokens = line.split()
            if len(tokens) < 4:
                continue
            if tokens[-1] not in {"+", "-"}:
                return None
            try:
                rows.append((round(float(tokens[1]), 4), float(tokens[3]), tokens[-1]))
            except ValueError:
                return None
    if not rows:
        return None
    return pd.DataFrame(
        rows, columns=["energy_key", "J", "parity"]
    ).drop_duplicates(["energy_key", "J"])


def classify_transitions(
    frame: pd.DataFrame, parities: pd.DataFrame | None
) -> pd.Series:
    """Label each line by the electric-dipole selection rule.

    E1 requires a parity change, ``ΔJ`` in {0, ±1}, and excludes ``J=0 -> J=0``.
    A ``ΔJ`` violation rules out E1 on its own, so lines can be identified as
    forbidden even for a source that publishes no parity.
    """
    if not {"J_upper", "J_lower"}.issubset(frame.columns):
        return pd.Series("unknown", index=frame.index, dtype=object)
    delta_j = (frame["J_upper"] - frame["J_lower"]).abs()
    angular_ok = (delta_j <= 1) & ~(
        (frame["J_upper"] == 0) & (frame["J_lower"] == 0)
    )
    if parities is None:
        return pd.Series(
            np.where(angular_ok, "E1_candidate", "forbidden"),
            index=frame.index, dtype=object,
        )
    lookup = parities.set_index(["energy_key", "J"])["parity"]
    upper = lookup.reindex(
        pd.MultiIndex.from_arrays([
            frame["energy_upper_cm-1"].round(4), frame["J_upper"],
        ])
    ).to_numpy()
    lower = lookup.reindex(
        pd.MultiIndex.from_arrays([
            frame["energy_lower_cm-1"].round(4), frame["J_lower"],
        ])
    ).to_numpy()
    known = pd.notna(upper) & pd.notna(lower)
    parity_changes = known & (upper != lower)
    return pd.Series(
        np.select(
            [angular_ok & parity_changes, angular_ok & ~known],
            ["E1_allowed", "E1_candidate"],
            default="forbidden",
        ),
        index=frame.index, dtype=object,
    )


def log_intensity_bounds(
    maxima: list[float], decades: float = SPECTRUM_FLOOR_DECADES
) -> tuple[float, float]:
    """Shared (floor, ceiling) bounds for a logarithmic absolute intensity axis."""
    positive = [
        float(value) for value in maxima
        if np.isfinite(value) and float(value) > 0
    ]
    if not positive:
        raise ValueError("no positive stick intensities to plot")
    return min(positive) * 10.0 ** -decades, max(positive)


def save_partition_plots(comparison: pd.DataFrame, output_dir: Path) -> None:
    """Write linear, logarithmic, and relative-error partition plots."""
    fig, axis = plt.subplots()
    for column, label in [
        ("Q_file", "file"), ("Q_pyexocross", "PyExoCross"),
        ("Q_reference", "reference"),
    ]:
        axis.plot(comparison["temperature_K"], comparison[column], label=label)
    axis.set(xlabel="Temperature (K)", ylabel="Q(T) [dimensionless]")
    axis.legend()
    fig.tight_layout()
    fig.savefig(output_dir / "partition_function.png", dpi=160)
    axis.set_yscale("log")
    fig.tight_layout()
    fig.savefig(output_dir / "partition_function_log.png", dpi=160)
    plt.close(fig)

    fig, axis = plt.subplots()
    axis.plot(comparison["temperature_K"], comparison["relative_difference"])
    axis.axhline(0, color="black", linewidth=0.7)
    axis.set(xlabel="Temperature (K)", ylabel="(Q_PyExoCross - Q_file) / Q_file")
    fig.tight_layout()
    fig.savefig(output_dir / "partition_function_relative_error.png", dpi=160)
    plt.close(fig)


def normalize_stick_output(
    raw_paths: list[Path],
    output_dir: Path,
    states_path: Path | None = None,
    top_n: int = 20,
) -> dict[str, Any]:
    """Normalize PyExoCross stick output and produce vacuum-axis plots."""
    candidates = [
        path for path in raw_paths
        if path.suffix.lower() in {".stick", ".txt", ".csv", ".dat"}
    ]
    frames: list[pd.DataFrame] = []
    for path in candidates:
        try:
            frame = pd.read_csv(path, sep=r"\s+", header=None, comment="#")
        except (OSError, pd.errors.ParserError):
            continue
        if frame.shape[1] >= 2:
            width = 6 if frame.shape[1] >= 6 else 2
            frame = frame.iloc[:, :width].copy()
            frame.columns = STICK_COLUMNS[:width]
            if width == 6 and not np.allclose(
                frame["energy_upper_cm-1"] - frame["energy_lower_cm-1"],
                frame["wavenumber_cm-1"], atol=1e-3, rtol=1e-6,
            ):
                # The state descriptors are not where this PyExoCross build put
                # them; drop them rather than classify transitions from them.
                frame = frame.iloc[:, :2]
            frame["temperature_source"] = path.name
            # Recorded numerically here so downstream plotting never has to
            # re-parse the PyExoCross file-naming convention.
            match = TEMPERATURE_PATTERN.search(path.name)
            frame["temperature_K"] = float(match.group(1)) if match else np.nan
            frames.append(frame)
    if not frames:
        return {"generated": False, "raw_files": [str(path) for path in raw_paths]}
    combined = pd.concat(frames, ignore_index=True)
    combined = combined[
        np.isfinite(combined["wavenumber_cm-1"])
        & (combined["wavenumber_cm-1"] > 0)
        & np.isfinite(combined["intensity"])
    ].copy()
    combined["vacuum_wavelength_nm"] = 1e7 / combined["wavenumber_cm-1"]
    parities = read_state_parities(states_path) if states_path else None
    combined["transition_type"] = classify_transitions(combined, parities)
    combined.to_csv(output_dir / "stick_spectrum.csv", index=False)
    strongest = combined.nlargest(top_n, "intensity")
    strongest.to_csv(output_dir / "strongest_lines.csv", index=False)
    # One panel per temperature, each on its own logarithmic scale. Temperatures
    # in one window can differ by tens of decades, so a shared axis would either
    # bury the coldest spectrum or turn the hottest into a solid block. Each
    # panel is also cut at its own floor: rendering every faint line of a
    # full-spectrum ion such as Fe III can exhaust memory, and those lines carry
    # no visible information. The CSV remains complete.
    panels = []
    for source, group in combined.groupby("temperature_source", sort=True):
        classes = list(group.groupby("transition_type", sort=True))
        try:
            # Anchoring on the weakest class keeps forbidden lines on the axis
            # even when they are orders of magnitude below the allowed ones.
            floor, ceiling = log_intensity_bounds(
                [subset["intensity"].max() for _, subset in classes]
            )
        except ValueError:
            continue
        visible = group[group["intensity"] >= floor].nlargest(20_000, "intensity")
        temperature = group["temperature_K"].dropna()
        label = f"{temperature.iloc[0]:g} K" if len(temperature) else str(source)
        order = float(temperature.iloc[0]) if len(temperature) else np.inf
        panels.append((order, label, visible, len(group), floor, ceiling))
    if not panels:
        panels = [(0.0, "no positive intensities", combined.iloc[:0], 0, 1e-1, 1.0)]
    panels.sort(key=lambda panel: (panel[0], panel[1]))
    colors = plt.get_cmap("viridis")(np.linspace(0.05, 0.75, len(panels)))
    for x, xlabel, name in [
        ("wavenumber_cm-1", "Wavenumber (cm$^{-1}$)", "stick_spectrum_wavenumber.png"),
        ("vacuum_wavelength_nm", "Vacuum wavelength (nm)", "stick_spectrum_wavelength.png"),
    ]:
        fig, axes = plt.subplots(
            len(panels), 1, sharex=True, squeeze=False,
            figsize=(9, 2.6 * len(panels) + 0.8),
        )
        for axis, (_, label, visible, total, floor, ceiling), color in zip(
            axes[:, 0], panels, colors
        ):
            # Forbidden lines are drawn last and in their own colour so they
            # stay legible where allowed lines are dense.
            for name_class, (class_label, class_color, width) in (
                TRANSITION_CLASSES.items()
            ):
                subset = visible[visible["transition_type"] == name_class]
                if subset.empty:
                    continue
                axis.vlines(
                    subset[x], floor, subset["intensity"], linewidth=0.5 * width,
                    color=class_color if class_color else color,
                    label=f"{class_label} ({len(subset)})",
                )
            axis.set_yscale("log")
            axis.set_ylim(floor, ceiling * 10.0)
            axis.set_ylabel(f"{label}\nline intensity")
            axis.grid(alpha=0.22)
            axis.legend(loc="upper right", fontsize="x-small")
            axis.text(
                0.995, 0.02, f"{len(visible)} lines above {floor:.1e} (of {total})",
                transform=axis.transAxes, ha="right", va="bottom", fontsize="small",
            )
        axes[0, 0].set_title("LTE stick spectrum (PyExoCross output units)")
        axes[-1, 0].set_xlabel(xlabel)
        fig.tight_layout()
        fig.savefig(output_dir / name, dpi=160)
        plt.close(fig)
    return {
        "generated": True,
        "line_count": len(combined),
        "raw_files": [str(path) for path in raw_paths],
        "wavelength_medium": "vacuum",
    }


def write_report(
    output_dir: Path,
    summary: dict[str, Any],
    sections: dict[str, Any],
) -> None:
    """Write the required structured summary and human-readable report."""
    (output_dir / "validation_summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=True), encoding="utf-8"
    )
    metadata = summary["metadata"]
    lines = [
        f"# PyExoCross validation: {metadata['element']} {metadata['spectroscopic_label']}",
        "",
        f"Overall status: **{summary['overall_status']}**",
        "",
        "## Inputs and versions",
        "",
        f"- Source: {metadata['source_database']}",
        f"- States: `{summary['inputs']['states']}`",
        f"- Transitions: `{summary['inputs']['trans']}`",
        f"- Partition function: `{summary['inputs']['pf']}`",
        f"- PyExoCross: {summary['pyexocross']['version']}",
        f"- Code version: {metadata['code_version']}",
        "",
        "## Charge-aware identity",
        "",
        f"- Z = {metadata['atomic_number']}",
        f"- Spectroscopic stage = {metadata['spectroscopic_stage']}",
        f"- Charge = +{metadata['charge']}",
        f"- Bound electrons = {metadata['electron_count']}",
        "",
        "## Validation categories",
        "",
        "- Format validation: independent parsers checked all three input schemas.",
        "- Internal consistency: state references and transition energy differences were checked.",
        "- Independent physical validation: only comparisons with genuinely external source paths qualify.",
        "- The Kurucz `.pf` and `.states` may share the same source/model; their agreement is a circular consistency check, not independent proof.",
        "",
    ]
    for heading, body in sections.items():
        lines.extend([f"## {heading}", "", str(body), ""])
    lines.extend([
        "## Scientific limitations",
        "",
        "- A finite states list can underestimate Q(T), especially at high temperature.",
        "- No ionization-limit truncation, continuum lowering, or occupation-probability model was supplied or invented.",
        "- Missing radiative channels make summed A too small and calculated lifetimes too long; metastable states are not automatically errors.",
        "- Cross-sections, when requested without measured broadening data, are exploratory format-validation products, not precision predictions.",
        "- Successful PyExoCross execution demonstrates compatibility and partial consistency, not completeness or universal physical correctness.",
        "",
    ])
    (output_dir / "validation_report.md").write_text("\n".join(lines), encoding="utf-8")
