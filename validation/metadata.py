"""Charge-aware atomic metadata and spectroscopic-stage helpers."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path
import re
import subprocess

from process_kurucz_atom import ELEMENTS


ROMAN_VALUES = {
    "I": 1,
    "V": 5,
    "X": 10,
    "L": 50,
    "C": 100,
    "D": 500,
    "M": 1000,
}

# Standard atomic weights (or representative mass number for elements without a
# standard atomic weight). The mass is used only by PyExoCross line profiles.
ATOMIC_MASSES = {
    "H": 1.008, "He": 4.002602, "Li": 6.94, "Be": 9.0121831,
    "B": 10.81, "C": 12.011, "N": 14.007, "O": 15.999, "F": 18.998403163,
    "Ne": 20.1797, "Na": 22.98976928, "Mg": 24.305, "Al": 26.9815385,
    "Si": 28.085, "P": 30.973761998, "S": 32.06, "Cl": 35.45,
    "Ar": 39.948, "K": 39.0983, "Ca": 40.078, "Sc": 44.955908,
    "Ti": 47.867, "V": 50.9415, "Cr": 51.9961, "Mn": 54.938044,
    "Fe": 55.845, "Co": 58.933194, "Ni": 58.6934, "Cu": 63.546,
    "Zn": 65.38, "Ga": 69.723, "Ge": 72.630, "As": 74.921595,
    "Se": 78.971, "Br": 79.904, "Kr": 83.798, "Rb": 85.4678,
    "Sr": 87.62, "Y": 88.90584, "Zr": 91.224, "Nb": 92.90637,
    "Mo": 95.95, "Tc": 98.0, "Ru": 101.07, "Rh": 102.90550,
    "Pd": 106.42, "Ag": 107.8682, "Cd": 112.414, "In": 114.818,
    "Sn": 118.710, "Sb": 121.760, "Te": 127.60, "I": 126.90447,
    "Xe": 131.293, "Cs": 132.90545196, "Ba": 137.327, "La": 138.90547,
    "Ce": 140.116, "Pr": 140.90766, "Nd": 144.242, "Pm": 145.0,
    "Sm": 150.36, "Eu": 151.964, "Gd": 157.25, "Tb": 158.92535,
    "Dy": 162.500, "Ho": 164.93033, "Er": 167.259, "Tm": 168.93422,
    "Yb": 173.045, "Lu": 174.9668, "Hf": 178.49, "Ta": 180.94788,
    "W": 183.84, "Re": 186.207, "Os": 190.23, "Ir": 192.217,
    "Pt": 195.084, "Au": 196.966569, "Hg": 200.592, "Tl": 204.38,
    "Pb": 207.2, "Bi": 208.98040, "Po": 209.0, "At": 210.0,
    "Rn": 222.0, "Fr": 223.0, "Ra": 226.0, "Ac": 227.0,
    "Th": 232.0377, "Pa": 231.03588, "U": 238.02891,
}


def roman_to_int(value: str) -> int:
    """Convert a canonical positive Roman numeral to an integer."""
    value = value.strip().upper()
    if not value or any(char not in ROMAN_VALUES for char in value):
        raise ValueError(f"invalid Roman numeral: {value!r}")
    total = 0
    previous = 0
    for char in reversed(value):
        number = ROMAN_VALUES[char]
        total += -number if number < previous else number
        previous = max(previous, number)
    if int_to_roman(total) != value:
        raise ValueError(f"non-canonical Roman numeral: {value!r}")
    return total


def int_to_roman(number: int) -> str:
    """Convert a positive integer to a canonical Roman numeral."""
    if number <= 0:
        raise ValueError("Roman numerals require a positive integer")
    pairs = (
        (1000, "M"), (900, "CM"), (500, "D"), (400, "CD"),
        (100, "C"), (90, "XC"), (50, "L"), (40, "XL"),
        (10, "X"), (9, "IX"), (5, "V"), (4, "IV"), (1, "I"),
    )
    result: list[str] = []
    for value, token in pairs:
        count, number = divmod(number, value)
        result.extend([token] * count)
    return "".join(result)


def code_version() -> str:
    """Return the current Git revision, without requiring a clean tree."""
    try:
        return subprocess.run(
            ["git", "rev-parse", "--short", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


@dataclass(frozen=True)
class IonMetadata:
    """Metadata whose stage/charge/electron relationships are validated."""

    element: str
    atomic_number: int
    spectroscopic_stage: int
    spectroscopic_label: str
    charge: int
    electron_count: int
    atomic_mass_da: float
    source_database: str
    source_file: str
    generation_timestamp: str
    code_version: str

    def as_dict(self) -> dict[str, object]:
        return asdict(self)


def parse_ion(
    ion: str,
    source_file: Path,
    charge: int | None = None,
    source_database: str = "Kurucz",
) -> IonMetadata:
    """Parse ``Fe-III``/``Fe_III`` and enforce stage = charge + 1."""
    match = re.fullmatch(r"([A-Z][a-z]?)[-_]([IVXLCDM]+)", ion.strip())
    if not match:
        raise ValueError(f"ion must look like Fe-III or Fe_III, got {ion!r}")
    element, numeral = match.groups()
    if element not in ELEMENTS:
        raise ValueError(f"unknown element symbol: {element}")
    stage = roman_to_int(numeral)
    inferred_charge = stage - 1
    if charge is not None and charge != inferred_charge:
        raise ValueError(
            f"{element} {numeral} has charge +{inferred_charge}, not +{charge}"
        )
    atomic_number = ELEMENTS[element]
    electron_count = atomic_number - inferred_charge
    if inferred_charge < 0 or inferred_charge >= atomic_number or electron_count <= 0:
        raise ValueError(
            f"invalid charge +{inferred_charge} for Z={atomic_number}"
        )
    filename_match = re.match(r"([A-Z][a-z]?)[_-]([IVXLCDM]+)__", source_file.name)
    if filename_match and filename_match.groups() != (element, numeral):
        raise ValueError(
            f"filename {source_file.name!r} disagrees with ion {element}-{numeral}"
        )
    parent = source_file.parent.name.replace("_", "-")
    if re.fullmatch(r"[A-Z][a-z]?-[IVXLCDM]+", parent) and parent != f"{element}-{numeral}":
        raise ValueError(f"directory {parent!r} disagrees with ion {element}-{numeral}")
    return IonMetadata(
        element=element,
        atomic_number=atomic_number,
        spectroscopic_stage=stage,
        spectroscopic_label=numeral,
        charge=inferred_charge,
        electron_count=electron_count,
        atomic_mass_da=ATOMIC_MASSES[element],
        source_database=source_database,
        source_file=str(source_file.resolve()),
        generation_timestamp=datetime.now(timezone.utc).isoformat(),
        code_version=code_version(),
    )
