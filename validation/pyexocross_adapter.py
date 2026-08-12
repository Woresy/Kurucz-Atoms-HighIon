"""Version-isolated adapter around the installed PyExoCross public API."""

from __future__ import annotations

from dataclasses import dataclass
import importlib
import json
import logging
import os
from pathlib import Path
import shutil
import time
from typing import Any

import numpy as np
import pandas as pd

from .calculations import reference_partition_function
from .exomol_validator import read_states
from .metadata import IonMetadata

LOGGER = logging.getLogger(__name__)


class PyExoCrossUnavailable(RuntimeError):
    """Raised when the optional runtime dependency cannot be imported."""


class PyExoCrossExecutionError(RuntimeError):
    """Raised with operation and input context when PyExoCross fails."""


@dataclass
class StagedDataset:
    """Canonical ExoAtom dataset layout expected by PyExoCross."""

    read_root: Path
    dataset_dir: Path
    atom: str
    dataset: str
    states_path: Path
    trans_path: Path
    pf_path: Path
    definition_path: Path


def pyexocross_info() -> dict[str, str]:
    """Return actual installed version and module path."""
    try:
        module = importlib.import_module("pyexocross")
    except ImportError as exc:
        raise PyExoCrossUnavailable(
            "PyExoCross is unavailable; install pyexocross before running this operation"
        ) from exc
    return {
        "version": str(getattr(module, "__version__", "unknown")),
        "module_path": str(Path(module.__file__).resolve()),
    }


# Field names and cfmt strings PyExoCross reads to lay out the .states columns.
# Names are the project's own; ExoAtom's Table 11 calls the first four ID, E,
# gtot and J, but renaming them here would change how PyExoCross labels its
# output columns, so the published names are carried in "desc" instead.
_STATES_FILE_FIELDS = [
    {"name": "i", "desc": "Unique integer identifier for the energy level (ExoAtom: ID)",
     "ffmt": "I12", "cfmt": "%12d"},
    {"name": "E", "desc": "State energy in cm^-1", "ffmt": "F12.6", "cfmt": "%12.6f"},
    {"name": "g", "desc": "State degeneracy (ExoAtom: gtot)", "ffmt": "I6", "cfmt": "%6d"},
    {"name": "J", "desc": "Total angular momentum quantum number", "ffmt": "F7.1", "cfmt": "%7.1f"},
    {"name": "unc", "desc": "Uncertainty in the state energy in cm^-1", "ffmt": "F12.6", "cfmt": "%12.6f"},
    {"name": "tau", "desc": "Radiative lifetime in s", "ffmt": "ES12.4", "cfmt": "%12.4E"},
    {"name": "gfac", "desc": "Lande g-factor (ExoAtom: gfactor)", "ffmt": "F10.6", "cfmt": "%10.6f"},
    {"name": "configuration", "desc": "Electronic configuration (ExoAtom: qn:configuration)",
     "ffmt": "A12", "cfmt": "%-12s"},
    {"name": "term", "desc": "Term symbol (ExoAtom: qn:LSCoupling)", "ffmt": "A7", "cfmt": "%-7s"},
    {"name": "source", "desc": "NI (measured, NIST) or CA (calculated)", "ffmt": "A2", "cfmt": "%2s"},
]

_TRANSITIONS_FILE_FIELDS = [
    {"name": "i", "desc": "Upper state ID", "ffmt": "I12", "cfmt": "%12d"},
    {"name": "f", "desc": "Lower state ID", "ffmt": "I12", "cfmt": "%12d"},
    {"name": "A", "desc": "Einstein A coefficient in s^-1", "ffmt": "ES10.4", "cfmt": "%10.4e"},
    {"name": "Wavenumber", "desc": "Transition wavenumber in cm^-1", "ffmt": "ES15.6", "cfmt": "%15.6e"},
]


def _definition(
    metadata: IonMetadata,
    states: pd.DataFrame | None = None,
    transition_count: int | None = None,
    max_wavenumber: float | None = None,
    max_temperature: float | None = None,
) -> dict[str, Any]:
    """Build the ExoAtom .adef.json for the staged dataset.

    The schema follows Ni et al. (2026), ExoAtom §3.5.1 / Table 11. PyExoCross
    1.1.9 reads only ``species.mass_in_Da``, ``dataset.predis`` and the
    ``dataset.states`` block, and parses the file with ``pd.read_json(orient=
    "columns")``, so every top-level key must map to an object. The remaining
    published fields are additive and describe the dataset for a human reader.
    """
    return {
        "species": {
            "atom": metadata.element,
            "ordinary_formula": metadata.element,
            "spectroscopic_notation": f"{metadata.element} {metadata.spectroscopic_label}",
            "charge": metadata.charge,
            "name": f"{metadata.element} {metadata.spectroscopic_label}",
            "mass_in_Da": metadata.atomic_mass_da,
        },
        "dataset": {
            "name": metadata.source_database,
            "doi": "",
            "max_temperature": max_temperature,
            "num_pressure_broadeners": 0,
            "nxsec_files": 0,
            "nkcoeff_files": 0,
            # Kurucz publishes gf values, not dipole matrix elements, and this
            # project derives neither cooling functions nor specific heats.
            "dipole_available": False,
            "cooling_function_available": False,
            "specific_heat_available": False,
            "Ionisation": None,
            # PyExoCross 1.1.9 keys, kept exactly as that reader expects them.
            "predis": False,
            "states": {
                "number_of_states": None if states is None else int(len(states)),
                "max_energy": None if states is None or states.empty
                else float(states["energy"].max()),
                "uncertainty_available": True,
                "uncertainties_available": True,
                "lifetime_available": True,
                "lande_g_available": True,
                "num_quanta": 2,
                "states_file_fields": _STATES_FILE_FIELDS,
            },
            "transitions": {
                "number_of_transitions": transition_count,
                "number_of_transition_files": 1,
                "max_wavenumber": max_wavenumber,
                "transitions_file_fields": _TRANSITIONS_FILE_FIELDS,
            },
        },
    }


def stage_local_files(
    states: Path,
    trans: Path,
    pf: Path,
    metadata: IonMetadata,
    root: Path,
    dataset: str = "Kurucz",
    states_frame: pd.DataFrame | None = None,
    transition_count: int | None = None,
    max_wavenumber: float | None = None,
    max_temperature: float | None = None,
) -> StagedDataset:
    """Copy local files into a reproducible, non-destructive ExoAtom layout."""
    atom = f"{metadata.element}_{metadata.spectroscopic_label}"
    dataset_dir = root / atom / dataset
    dataset_dir.mkdir(parents=True, exist_ok=True)
    prefix = f"{atom}__{dataset}"
    staged_states = dataset_dir / f"{prefix}.states"
    staged_trans = dataset_dir / f"{prefix}.trans"
    staged_pf = dataset_dir / f"{prefix}.pf"
    definition = dataset_dir / f"{prefix}.adef.json"
    # The staged .states and .trans are only ever read, so a hard link avoids
    # duplicating multi-GB files. The staged .pf is rewritten in place by
    # _ensure_pf_temperatures, so it must be an independent copy: a hard link
    # would make that write reach through to the source dataset.
    for source, target, may_link in [
        (states, staged_states, True),
        (trans, staged_trans, True),
        (pf, staged_pf, False),
    ]:
        if source.resolve() == target.resolve():
            continue
        if target.exists():
            target.unlink()
        if may_link:
            try:
                os.link(source, target)
                continue
            except OSError:
                pass
        shutil.copy2(source, target)
    definition.write_text(
        json.dumps(
            _definition(
                metadata, states_frame, transition_count, max_wavenumber, max_temperature
            ),
            indent=2,
        ),
        encoding="utf-8",
    )
    return StagedDataset(
        read_root=root,
        dataset_dir=dataset_dir,
        atom=atom,
        dataset=dataset,
        states_path=staged_states,
        trans_path=staged_trans,
        pf_path=staged_pf,
        definition_path=definition,
    )


class PyExoCrossAdapter:
    """Small public-API facade that records calls and isolates v1.1.x layout."""

    def __init__(
        self,
        staged: StagedDataset,
        output_dir: Path,
        log_path: Path,
    ) -> None:
        info = pyexocross_info()
        self.version = info["version"]
        self.staged = staged
        self.output_dir = output_dir.resolve()
        self.log_path = log_path.resolve()
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.calls: list[dict[str, object]] = []

    def _base_kwargs(self) -> dict[str, object]:
        return {
            "database": "ExoAtom",
            "atom": self.staged.atom,
            "dataset": self.staged.dataset,
            "read_path": str(self.staged.read_root.resolve()) + "/",
            "save_path": str(self.output_dir) + "/",
            "logs_path": str(self.log_path),
            "ncputrans": 1,
            "ncpufiles": 1,
            "chunk_size": 100_000,
            "cache": "none",
        }

    def _invoke(self, operation: str, **parameters: object) -> None:
        try:
            px = importlib.import_module("pyexocross")
            function = getattr(px, operation)
            kwargs = self._base_kwargs()
            kwargs.update(parameters)
            LOGGER.info("PyExoCross %s %s", operation, kwargs)
            function(**kwargs)
            self.calls.append({
                "operation": operation,
                "status": "success",
                "parameters": {key: str(value) for key, value in kwargs.items()},
            })
        except Exception as exc:
            self.calls.append({
                "operation": operation,
                "status": "failed",
                "error": f"{type(exc).__name__}: {exc}",
            })
            raise PyExoCrossExecutionError(
                f"PyExoCross {operation} failed for {self.staged.dataset_dir}: {exc}"
            ) from exc

    def calculate_partition_function(
        self, temperatures_k: np.ndarray
    ) -> tuple[np.ndarray, list[Path]]:
        """Use the public high-level API and return its saved Q values."""
        temperatures = np.asarray(temperatures_k, dtype=float)
        rounded = np.round(temperatures).astype(int)
        if not np.allclose(temperatures, rounded):
            raise ValueError("PyExoCross 1.1.9 PF API requires integer temperature grid")
        if np.any(rounded < 1):
            raise ValueError("PyExoCross 1.1.9 PF API requires temperatures of 1 K or more")
        # The API only exposes the uniform grid Ntemp, 2*Ntemp, ... Tmax, so
        # ntemp == tmax == T makes that grid exactly {T}. One call per requested
        # temperature therefore costs one state sum each. Deriving a single
        # shared step from the requested values instead makes the cost Tmax/step,
        # which explodes to Tmax points whenever they are mutually prime.
        values: list[float] = []
        archived: list[Path] = []
        for temperature in rounded:
            for stale in self.output_dir.rglob("*.pf"):
                if stale not in archived:
                    stale.unlink()
            self._invoke(
                "partition_functions",
                ntemp=int(temperature),
                tmax=int(temperature),
            )
            produced = [
                path for path in self.output_dir.rglob("*.pf") if path not in archived
            ]
            if len(produced) != 1:
                raise PyExoCrossExecutionError(
                    f"expected one .pf file for {temperature} K, found {len(produced)}"
                )
            raw = pd.read_csv(produced[0], sep=r"\s+", header=None, names=["T", "Q"])
            match = raw.loc[raw["T"].round().astype(int) == int(temperature), "Q"]
            if match.empty:
                raise PyExoCrossExecutionError(
                    f"{temperature} K absent from {produced[0]}"
                )
            values.append(float(match.iloc[-1]))
            # Keep every temperature's raw output; each call overwrites the same
            # PyExoCross filename otherwise.
            target = produced[0].with_name(
                f"{produced[0].stem}_T{int(temperature)}K.pf"
            )
            produced[0].replace(target)
            archived.append(target)
        return np.array(values, dtype=float), archived

    def calculate_lifetimes(self) -> list[Path]:
        """Run the public lifetime workflow and return new raw outputs."""
        before = set(self.output_dir.rglob("*.states"))
        self._invoke("lifetimes", compress=False)
        generated = sorted((self.output_dir / "lifetime").rglob("*.states"))
        if generated:
            return generated
        return sorted(set(self.output_dir.rglob("*.states")) - before)

    def calculate_stick_spectrum(
        self,
        temperatures_k: list[float],
        min_wavenumber: float,
        max_wavenumber: float,
        absorption_emission: str = "Ab",
    ) -> list[Path]:
        """Run LTE stick spectra in vacuum wavenumber space."""
        self._ensure_pf_temperatures(temperatures_k)
        before = set(self.output_dir.rglob("*"))
        # PyExoCross encodes the wavenumber range in the filename, so a rerun
        # over a different range leaves the previous range's files in place.
        # Anything not written by this call must be excluded or the two ranges
        # silently merge into one spectrum.
        started = time.time()
        # PyExoCross 1.1.9 can construct a Q array shorter than T_list for some
        # sparse/non-uniform .pf grids, then raises IndexError on the third
        # temperature. One public call per temperature avoids that version bug
        # and produces the same documented LTE result without patching the package.
        for temperature in temperatures_k:
            self._invoke(
                "stick_spectra",
                temperatures=[temperature],
                min_range=min_wavenumber,
                max_range=max_wavenumber,
                wn_wl="WN",
                wn_wl_unit="cm-1",
                abs_emi=absorption_emission,
                nlte_method="L",
                plot=False,
            )
        generated = sorted(
            path for path in self.output_dir.rglob("*.stick")
            if path.stat().st_mtime >= started
        )
        if generated:
            return generated
        return sorted(
            path for path in set(self.output_dir.rglob("*")) - before if path.is_file()
        )

    def _ensure_pf_temperatures(self, temperatures_k: list[float]) -> None:
        """Add explicitly interpolated in-range Q values to the staged copy.

        PyExoCross 1.1.9 uses exact ``isin`` matching rather than interpolation.
        The project requirement specifies linear interpolation and forbids
        extrapolation, so only the disposable staged copy is augmented.
        """
        frame = pd.read_csv(
            self.staged.pf_path, sep=r"\s+", header=None, names=["T", "Q"]
        )
        existing = set(frame["T"].astype(float))
        additions: list[dict[str, float]] = []
        for temperature in temperatures_k:
            value = float(temperature)
            if value in existing:
                continue
            if value < frame["T"].min() or value > frame["T"].max():
                raise ValueError(
                    f"cannot extrapolate partition function to {value:g} K"
                )
            additions.append({
                "T": value,
                "Q": float(np.interp(value, frame["T"], frame["Q"])),
            })
        if additions:
            frame = pd.concat([frame, pd.DataFrame(additions)], ignore_index=True)
            frame = frame.sort_values("T").drop_duplicates("T")
            # Write a new inode and move it into place rather than truncating the
            # staged file. os.replace rebinds only this name, so the augmented
            # copy can never propagate through a link to a shared original.
            scratch = self.staged.pf_path.with_suffix(".pf.tmp")
            np.savetxt(scratch, frame[["T", "Q"]].to_numpy(), fmt="%12.6f %15.8e")
            os.replace(scratch, self.staged.pf_path)
            self.calls.append({
                "operation": "linear_pf_interpolation_for_staged_input",
                "status": "success",
                "temperatures_K": [row["T"] for row in additions],
                "extrapolation": False,
            })

    def calculate_cross_section(
        self,
        temperature_k: float,
        pressure_bar: float,
        min_wavenumber: float,
        max_wavenumber: float,
        profile: str = "Gaussian",
        bin_size: float = 0.2,
        cutoff: float | None = None,
    ) -> list[Path]:
        """Run an explicitly exploratory profile-based cross-section."""
        before = set(self.output_dir.rglob("*"))
        self._invoke(
            "cross_sections",
            temperatures=[temperature_k],
            pressures=[pressure_bar],
            min_range=min_wavenumber,
            max_range=max_wavenumber,
            wn_wl="WN",
            wn_wl_unit="cm-1",
            profile=profile,
            bin_size=bin_size,
            cutoff=cutoff,
            abs_emi="Ab",
            nlte_method="L",
            plot=False,
        )
        generated = sorted(self.output_dir.rglob("*.xsec"))
        if generated:
            return generated
        return sorted(
            path for path in set(self.output_dir.rglob("*")) - before if path.is_file()
        )

    def low_level_partition_function(
        self, temperatures_k: np.ndarray
    ) -> np.ndarray:
        """Call PyExoCross's published calculation function as a precise check."""
        try:
            module = importlib.import_module(
                "pyexocross.calculation.calculate_partition_func"
            )
        except ImportError as exc:
            raise PyExoCrossUnavailable(str(exc)) from exc
        states = read_states(self.staged.states_path)
        values = np.array([
            module.cal_partition_func(
                states["energy"].to_numpy(float),
                states["g"].to_numpy(int),
                float(temperature),
            )
            for temperature in temperatures_k
        ])
        reference = reference_partition_function(
            states["energy"].to_numpy(float),
            states["g"].to_numpy(int),
            np.asarray(temperatures_k, dtype=float),
        )
        # PyExoCross 1.1.9 carries a few more digits of c2 than the requested
        # published sanity-check constant, so agreement is expected at ~1e-9,
        # not bit-for-bit.
        if not np.allclose(values, reference, rtol=1e-9, atol=1e-12):
            raise PyExoCrossExecutionError(
                "PyExoCross low-level PF disagrees with independent reference"
            )
        return values
