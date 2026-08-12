from __future__ import annotations

import importlib
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from validation.calculations import (
    calculate_lifetimes,
    interpolate_in_bounds,
    reference_partition_function,
)
from validation.exomol_validator import (
    read_states,
    validate_pf,
    validate_states,
    validate_transitions,
)
from validation.metadata import int_to_roman, parse_ion, roman_to_int
from validation.discovery import discover_ion_files
from validation.cross_source import match_states, read_nist_states
from validation.pyexocross_adapter import PyExoCrossUnavailable, pyexocross_info
from validation.run_pyexocross import _to_wavenumber_range

FIXTURES = Path(__file__).parent / "fixtures"
STATES = FIXTURES / "C_III__Kurucz.states"
TRANS = FIXTURES / "C_III__Kurucz.trans"
PF = FIXTURES / "C_III__Kurucz.pf"


def test_roman_stage_charge_mapping() -> None:
    assert roman_to_int("III") == 3
    assert int_to_roman(4) == "IV"
    metadata = parse_ion("C-III", STATES, charge=2)
    assert metadata.spectroscopic_stage == 3
    assert metadata.charge == 2
    assert metadata.electron_count == 4
    with pytest.raises(ValueError, match="not"):
        parse_ion("C-III", STATES, charge=3)


def test_states_parser_and_schema() -> None:
    result, states = validate_states(STATES)
    assert result.passed
    assert list(states.columns[:4]) == ["id", "energy", "g", "J"]
    assert states["configuration"].tolist() == ["2s2", "2s2p", "2s3d"]
    assert result.metrics["energy_unit"] == "cm^-1 (project schema)"


def test_transition_id_mapping_and_energy_difference() -> None:
    states = read_states(STATES)
    result = validate_transitions(TRANS, states)
    assert result.passed
    assert result.metrics["referenced_state_missing_count"] == 0
    assert result.metrics["wavenumber_max_error_cm-1"] == 0


def test_missing_transition_reference(tmp_path: Path) -> None:
    path = tmp_path / "bad.trans"
    path.write_text("9 1 1.0e1 1000.0\n", encoding="utf-8")
    result = validate_transitions(path, read_states(STATES))
    assert not result.passed
    assert result.metrics["referenced_state_missing_count"] == 1


def test_invalid_A_coefficient(tmp_path: Path) -> None:
    path = tmp_path / "bad.trans"
    path.write_text("2 1 -1.0 1000.0\n", encoding="utf-8")
    result = validate_transitions(path, read_states(STATES))
    assert not result.passed
    assert result.metrics["invalid_A_count"] == 1


def test_pf_temperature_sorting(tmp_path: Path) -> None:
    good, _ = validate_pf(PF, ground_degeneracy=1)
    assert good.passed
    path = tmp_path / "bad.pf"
    path.write_text("1000 2\n100 1\n", encoding="utf-8")
    bad, _ = validate_pf(path)
    assert not bad.passed
    assert any("strictly increasing" in error for error in bad.errors)


def test_reference_partition_function() -> None:
    q = reference_partition_function(
        np.array([0.0, 1000.0]), np.array([1, 3]), np.array([1000.0])
    )
    expected = 1 + 3 * np.exp(-1.4387768775039338)  # c2 = h*c/k_B, exact under SI 2019
    assert q[0] == pytest.approx(expected)


def test_interpolation_boundaries() -> None:
    assert interpolate_in_bounds(
        np.array([100.0, 200.0]), np.array([1.0, 3.0]), np.array([150.0])
    )[0] == pytest.approx(2.0)
    with pytest.raises(ValueError, match="outside"):
        interpolate_in_bounds(
            np.array([100.0, 200.0]), np.array([1.0, 3.0]), np.array([99.0])
        )


def test_vacuum_wavelength_range_conversion() -> None:
    assert _to_wavenumber_range((100.0, 200.0), "nm") == (50_000.0, 100_000.0)
    with pytest.raises(ValueError):
        _to_wavenumber_range((200.0, 100.0), "nm")


def test_lifetime_calculation_and_missing_channels() -> None:
    output = calculate_lifetimes(read_states(STATES), TRANS)
    by_id = output.set_index("id")
    assert by_id.loc[2, "calculated_lifetime_s"] == pytest.approx(0.1)
    assert by_id.loc[3, "calculated_lifetime_s"] == pytest.approx(0.05)
    assert by_id.loc[1, "status"] == "no_valid_downward_transition"
    assert np.isnan(by_id.loc[1, "calculated_lifetime_s"])


def test_pyexocross_unavailable_error(monkeypatch: pytest.MonkeyPatch) -> None:
    real_import = importlib.import_module

    def fake_import(name: str, *args: object, **kwargs: object):
        if name == "pyexocross":
            raise ImportError("missing")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(importlib, "import_module", fake_import)
    with pytest.raises(PyExoCrossUnavailable, match="install"):
        pyexocross_info()


def test_installed_pyexocross_version() -> None:
    info = pyexocross_info()
    assert info["version"] == "1.1.9"
    assert info["module_path"].endswith("pyexocross/__init__.py")


def test_minimal_data_are_mutually_consistent() -> None:
    states_result, states = validate_states(STATES)
    trans_result = validate_transitions(TRANS, states)
    pf_result, pf = validate_pf(PF, ground_degeneracy=1)
    assert states_result.passed and trans_result.passed and pf_result.passed
    q = reference_partition_function(
        states["energy"].to_numpy(float),
        states["g"].to_numpy(int),
        pf["temperature"].to_numpy(float),
    )
    assert np.allclose(q, pf["Q"], rtol=2e-7)


def test_generic_file_discovery_and_identity(tmp_path: Path) -> None:
    kurucz = tmp_path / "kurucz" / "Ti-IV"
    nist = tmp_path / "nist" / "nested"
    kurucz.mkdir(parents=True)
    nist.mkdir(parents=True)
    for suffix in ("states", "trans", "pf"):
        (kurucz / f"Ti_IV__Kurucz.{suffix}").touch()
        (nist / f"Ti_IV__NIST.{suffix}").touch()
    ion = discover_ion_files("Ti-IV", tmp_path / "kurucz", tmp_path / "nist")
    assert ion.charge == 3
    assert ion.electron_count == 19
    assert ion.kurucz_files.availability.value == "AVAILABLE"
    assert ion.nist_files.availability.value == "AVAILABLE"


def test_quantum_matching_does_not_use_database_ids(tmp_path: Path) -> None:
    nist_path = tmp_path / "Fe_III__NIST.states"
    nist_path.write_text(
        "99 0.0 9 4 0.1 3d6 5D +\n"
        "77 436.0 7 3 0.1 3d6 5D +\n",
        encoding="utf-8",
    )
    nist = read_nist_states(nist_path)
    kurucz = pd.DataFrame({
        "id": [1, 2], "energy": [0.02, 435.9], "g": [9, 7],
        "J": [4.0, 3.0], "configuration": ["d6", "d6"],
        "term": ["5D", "5D"],
    })
    matches = match_states(kurucz, nist)
    assert set(matches["nist_id"]) == {99, 77}
    assert set(matches["kurucz_id"]) == {1, 2}
    assert set(matches["classification"]) <= {
        "EXACT_QUANTUM_MATCH", "HIGH_CONFIDENCE"
    }
