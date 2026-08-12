"""Generic atomic-ion identity and cross-database file discovery."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from enum import Enum
from pathlib import Path

from .metadata import ATOMIC_MASSES, int_to_roman, parse_ion


class Availability(str, Enum):
    AVAILABLE = "AVAILABLE"
    PARTIAL = "PARTIAL"
    NOT_AVAILABLE = "NOT_AVAILABLE"
    INVALID = "INVALID"


@dataclass(frozen=True)
class IonFiles:
    states: Path | None = None
    trans: Path | None = None
    pf: Path | None = None

    @property
    def availability(self) -> Availability:
        present = [path is not None and path.is_file() for path in (self.states, self.trans, self.pf)]
        if all(present):
            return Availability.AVAILABLE
        if any(present):
            return Availability.PARTIAL
        return Availability.NOT_AVAILABLE

    def as_dict(self) -> dict[str, str | None]:
        return {
            "states": str(self.states.resolve()) if self.states else None,
            "trans": str(self.trans.resolve()) if self.trans else None,
            "pf": str(self.pf.resolve()) if self.pf else None,
            "availability": self.availability.value,
        }


@dataclass(frozen=True)
class ResolvedIon:
    element_symbol: str
    atomic_number: int
    spectroscopic_stage: int
    roman_stage: str
    charge: int
    electron_count: int
    display_name: str
    atomic_mass_da: float
    kurucz_files: IonFiles
    nist_files: IonFiles

    def as_dict(self) -> dict[str, object]:
        result = asdict(self)
        result["kurucz_files"] = self.kurucz_files.as_dict()
        result["nist_files"] = self.nist_files.as_dict()
        return result


def _find_component(root: Path, stem: str, source: str, extension: str) -> Path | None:
    filename = f"{stem}__{source}.{extension}"
    direct_candidates = (
        root / stem.replace("_", "-") / filename,
        root / filename,
    )
    for candidate in direct_candidates:
        if candidate.is_file():
            return candidate
    matches = sorted(root.rglob(filename)) if root.exists() else []
    return matches[0] if matches else None


def discover_ion_files(
    ion: str,
    kurucz_root: Path,
    nist_root: Path,
) -> ResolvedIon:
    """Resolve an arbitrary ``Element-ROMAN`` ion under both data roots."""
    probe = Path(f"{ion.replace('-', '_')}__Kurucz.states")
    identity = parse_ion(ion, probe)
    roman = int_to_roman(identity.spectroscopic_stage)
    stem = f"{identity.element}_{roman}"

    def files(root: Path, source: str) -> IonFiles:
        return IonFiles(**{
            extension: _find_component(root, stem, source, extension)
            for extension in ("states", "trans", "pf")
        })

    return ResolvedIon(
        element_symbol=identity.element,
        atomic_number=identity.atomic_number,
        spectroscopic_stage=identity.spectroscopic_stage,
        roman_stage=roman,
        charge=identity.charge,
        electron_count=identity.electron_count,
        display_name=f"{identity.element} {roman}",
        atomic_mass_da=ATOMIC_MASSES[identity.element],
        kurucz_files=files(kurucz_root, "Kurucz"),
        nist_files=files(nist_root, "NIST"),
    )


def discover_complete_kurucz_ions(root: Path) -> list[str]:
    """Return every ion with a complete Kurucz states/trans/pf triplet."""
    ions: list[str] = []
    for states in root.rglob("*__Kurucz.states"):
        prefix = states.name.removesuffix("__Kurucz.states")
        if (
            states.with_name(prefix + "__Kurucz.trans").is_file()
            and states.with_name(prefix + "__Kurucz.pf").is_file()
        ):
            ions.append(prefix.replace("_", "-"))
    return sorted(set(ions))
