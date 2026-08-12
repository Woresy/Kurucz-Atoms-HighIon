import argparse
import bisect
import csv
import gzip
import itertools
import math
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from io import StringIO
from pathlib import Path

import pandas as pd


BASE_URL = "http://kurucz.harvard.edu/atoms"

ELEMENTS = {
    "H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5, "C": 6, "N": 7, "O": 8,
    "F": 9, "Ne": 10, "Na": 11, "Mg": 12, "Al": 13, "Si": 14, "P": 15,
    "S": 16, "Cl": 17, "Ar": 18, "K": 19, "Ca": 20, "Sc": 21, "Ti": 22,
    "V": 23, "Cr": 24, "Mn": 25, "Fe": 26, "Co": 27, "Ni": 28, "Cu": 29,
    "Zn": 30, "Ga": 31, "Ge": 32, "As": 33, "Se": 34, "Br": 35, "Kr": 36,
    "Rb": 37, "Sr": 38, "Y": 39, "Zr": 40, "Nb": 41, "Mo": 42, "Tc": 43,
    "Ru": 44, "Rh": 45, "Pd": 46, "Ag": 47, "Cd": 48, "In": 49, "Sn": 50,
    "Sb": 51, "Te": 52, "I": 53, "Xe": 54, "Cs": 55, "Ba": 56, "La": 57,
    "Ce": 58, "Pr": 59, "Nd": 60, "Pm": 61, "Sm": 62, "Eu": 63, "Gd": 64,
    "Tb": 65, "Dy": 66, "Ho": 67, "Er": 68, "Tm": 69, "Yb": 70, "Lu": 71,
    "Hf": 72, "Ta": 73, "W": 74, "Re": 75, "Os": 76, "Ir": 77, "Pt": 78,
    "Au": 79, "Hg": 80, "Tl": 81, "Pb": 82, "Bi": 83, "Po": 84, "At": 85,
    "Rn": 86, "Fr": 87, "Ra": 88, "Ac": 89, "Th": 90, "Pa": 91, "U": 92,
}
NUMBER_TO_ELEMENT = {number: symbol for symbol, number in ELEMENTS.items()}

LEVEL_RE = re.compile(r"^\d+\.\d+(?:EVE|ODD|ERz|ORz|EPo|OPo)$")

# Kurucz match tables are laid out as six fixed-width cells per 72-character line.
MAPPING_CELL_WIDTH = 12


@dataclass
class FileGroup:
    stem: str
    gam: str
    lines: str | None
    agafgf: str | None


@dataclass
class IonDiscovery:
    element: str
    charge: int
    code: str
    ion_name: str
    url: str
    exists: bool
    groups: list[FileGroup]
    life: str | None
    pf: str | None
    reason: str = ""
    # Deduplicated list of agafgf-type files used to build transitions. Transitions
    # are built from these files alone (each contains the full transition endpoints),
    # decoupled from the gam groups so a single combined agafgf serving several gam
    # suffixes is processed exactly once.
    trans_sources: list[str] = field(default_factory=list)
    # .lines files usable to top up an incomplete agafgf. Per-gam-stem files come
    # first; the code-level gf{code}.lines is a fallback for ions that publish the
    # merged line list under the bare code while the .gam files carry suffixes.
    lines_sources: list[str] = field(default_factory=list)


def roman(number: int) -> str:
    values = [
        (1000, "M"), (900, "CM"), (500, "D"), (400, "CD"), (100, "C"),
        (90, "XC"), (50, "L"), (40, "XL"), (10, "X"), (9, "IX"),
        (5, "V"), (4, "IV"), (1, "I"),
    ]
    result = []
    for value, symbol in values:
        while number >= value:
            result.append(symbol)
            number -= value
    return "".join(result)


def ion_name(element: str, charge: int) -> str:
    return f"{element}-{roman(charge + 1)}"


def ion_code(element: str, charge: int) -> str:
    return f"{ELEMENTS[element]:02d}{charge:02d}"


def code_to_element_charge(code: str) -> tuple[str, int] | None:
    if not re.fullmatch(r"\d{4}", code):
        return None
    atomic_number = int(code[:2])
    charge = int(code[2:])
    element = NUMBER_TO_ELEMENT.get(atomic_number)
    if element is None:
        return None
    return element, charge


# Records every file the Kurucz server returns HTTP 403 (Forbidden) for. A 403
# means the file exists but is access-restricted (worth asking an admin to unlock),
# unlike 404 which means it genuinely does not exist. Written out as a CSV at the
# end of a run via write_blocked().
BLOCKED: list[dict[str, object]] = []

# Records every ion whose resolved transition source carries far fewer records
# than its own .gam header declares. Kurucz publishes several products per ion
# (gf{code}, gf{code}y/z, gfemq{code}, .posagafgf subsets ...) and only some of
# them ship an .agafgf; when none matches, discovery falls back to whatever
# agafgf-shaped file exists, which can be a tiny subset of a different product.
# That yields a well-formed but near-empty .trans, so the shortfall is recorded
# here and written out by write_incomplete().
INCOMPLETE: list[dict[str, object]] = []

# A resolved transition source below this fraction of the declared count is
# reported. Legitimate losses (unmappable states, filtered rows) are small; the
# observed failures are three or more orders of magnitude.
TRANSITION_COMPLETENESS_THRESHOLD = 0.9

# These ions are published by Kurucz as complementary letter-suffixed blocks.
# The original ExoAtom workflow explicitly merged them; treating an unsuffixed
# file as authoritative would silently discard part of the line list.
SPLIT_ION_CODES = {"1100", "1900", "2001", "3000", "3901"}


def record_blocked(element: str, charge: int, code: str, ion: str, role: str, filename: str, url: str) -> None:
    BLOCKED.append({
        "element": element, "charge": charge, "code": code, "ion": ion,
        "role": role, "filename": filename, "url": url, "http_status": 403,
    })


def record_blocked_file(discovery: "IonDiscovery", filename: str, role: str) -> None:
    record_blocked(
        discovery.element, discovery.charge, discovery.code, discovery.ion_name,
        role, filename, urllib.parse.urljoin(discovery.url, filename),
    )


def check_transition_completeness(
    discovery: "IonDiscovery",
    declared_lines: dict[str, int],
    trans_rows: int,
    dropped: int,
) -> None:
    """Compare records read from the transition source against the .gam header.

    The level file states how many transitions its own computation produced, so
    a resolved source that supplies far fewer records is reading the wrong or a
    partial product rather than losing rows to normal filtering. Compare against
    records *read* (written plus dropped) so legitimately unmappable rows are
    not counted as missing data.
    """
    declared_total = sum(declared_lines.values())
    if not declared_total:
        return
    records_read = trans_rows + dropped
    ratio = records_read / declared_total
    if ratio >= TRANSITION_COMPLETENESS_THRESHOLD:
        return
    print(
        f"  WARNING: transitions incomplete -- read {records_read} of "
        f"{declared_total} declared ({ratio:.4%}); sources: "
        f"{', '.join(discovery.trans_sources) or 'none'}"
    )
    INCOMPLETE.append({
        "element": discovery.element,
        "charge": discovery.charge,
        "code": discovery.code,
        "ion": discovery.ion_name,
        "gam_files": ";".join(sorted(declared_lines)),
        "declared_transitions": declared_total,
        "records_read": records_read,
        "transitions_written": trans_rows,
        "dropped": dropped,
        "completeness": f"{ratio:.6f}",
        "trans_sources": ";".join(discovery.trans_sources),
    })


def is_gz(filename: str | None) -> bool:
    return filename is not None and (filename.endswith("-gz") or filename.endswith(".gz"))


def pick_first(files: set[str], candidates: list[str]) -> str | None:
    for candidate in candidates:
        if candidate in files:
            return candidate
    return None


def read_url(url: str, timeout: int) -> str:
    with urllib.request.urlopen(url, timeout=timeout) as response:
        return response.read().decode("utf-8", errors="replace")


def list_directory(url: str, timeout: int) -> list[str] | None:
    try:
        html = read_url(url, timeout)
    except urllib.error.HTTPError as exc:
        if exc.code == 404:
            return None
        raise
    return [urllib.parse.unquote(item) for item in re.findall(r'href="([^"]+)"', html)]


def discover_ion(element: str, charge: int, timeout: int, base_url: str) -> IonDiscovery:
    code = ion_code(element, charge)
    name = ion_name(element, charge)
    url = f"{base_url.rstrip('/')}/{code}/"
    try:
        names = list_directory(url, timeout)
    except urllib.error.HTTPError as exc:
        if exc.code == 403:
            record_blocked(element, charge, code, name, "directory", "", url)
            return IonDiscovery(element, charge, code, name, url, False, [], None, None, "directory forbidden (403)")
        raise
    if names is None:
        return IonDiscovery(element, charge, code, name, url, False, [], None, None, "directory not found")

    files = {name for name in names if not name.startswith("?") and not name.endswith("/")}
    gam_files = sorted(name for name in files if re.fullmatch(rf"gf{code}[A-Za-z0-9]*\.gam", name))

    simple = f"gf{code}.gam"
    if code in SPLIT_ION_CODES:
        split_gam_files = [
            filename for filename in gam_files
            if re.fullmatch(rf"gf{code}[A-Za-z]+\.gam", filename)
        ]
        if split_gam_files:
            gam_files = split_gam_files
    elif simple in gam_files:
        gam_files = [simple]

    # Drop digit-suffixed gam files (e.g. gf03006.gam, gf050010.gam) when a
    # non-digit variant exists. These are alternate *computations* of the same
    # level set (identical J/label, slightly different energies) rather than the
    # complementary parity/block splits carried by letter suffixes (w/y/z). Merging
    # them would double-count every level with a different energy, inflating the
    # states. Letter-suffixed files (complementary) are kept and still merged.
    def _gam_suffix(name: str) -> str:
        return name[len(f"gf{code}"):-len(".gam")]
    non_digit = [g for g in gam_files if not _gam_suffix(g).isdigit()]
    if non_digit and len(non_digit) < len(gam_files):
        gam_files = non_digit

    groups = []
    for gam in gam_files:
        stem = gam[:-4]
        # Per-group resolution of the transition source files, kept at the ORIGINAL
        # priority so any ion that already resolves an agafgf keeps the exact same
        # file (no regression to previously-produced outputs). The standard ".lines"
        # / ".agafgf" names cover most ions; gzip-compressed (-gz) and ".all" /
        # ".allagafgf" variants share the exact column layout, so only the file name
        # and decompression differ. Code-level singleton fallbacks (".all" etc.) only
        # apply to the single-gam case to avoid double-counting across suffix groups.
        lines_candidates = [f"{stem}.lines", f"{stem}.lines-gz"]
        agafgf_candidates = [f"{stem}.agafgf", f"{stem}.agafgf-gz"]
        if stem == simple[:-4]:
            lines_candidates += [f"{stem}.all", f"gfemq{code}.all"]
            agafgf_candidates += [f"{stem}.allagafgf", f"gfemq{code}.allagafgf"]
        agafgf = pick_first(files, agafgf_candidates)
        # A .lines file is only row-aligned with an agafgf from the SAME product
        # (same name before the first dot). Pairing e.g. gfemq1600.allagafgf with
        # gf1600.lines zips two unrelated line lists and corrupts every row; an
        # unpaired extended agafgf is parsed standalone from its own endpoints.
        if agafgf:
            agafgf_product = agafgf.split(".", 1)[0]
            lines_candidates = [c for c in lines_candidates if c.split(".", 1)[0] == agafgf_product]
        lines = pick_first(files, lines_candidates)
        groups.append(FileGroup(stem, gam, lines, agafgf))

    # Transitions are built from the deduplicated set of per-group agafgf files. For
    # ions that already resolved an agafgf above this is identical to the per-group
    # files (the dedup is a no-op when they are distinct), so existing outputs are
    # unchanged.
    trans_sources = list(dict.fromkeys(g.agafgf for g in groups if g.agafgf))
    if not trans_sources:
        # Only-when-nothing-found fallback: some (mostly higher) ions store the agafgf
        # under a name not tied to the gam stem -- the gfemq variant, a single combined
        # gf{code}.agafgf shared by several gam suffixes, or a .posagafgf subset. This
        # branch never runs when a per-group agafgf was found, so it cannot affect any
        # ion that currently produces transitions.
        fallback = pick_first(files, [
            f"gf{code}.agafgf", f"gf{code}.agafgf-gz",
            f"gfemq{code}.agafgf", f"gfemq{code}.agafgf-gz",
            f"gf{code}.allagafgf", f"gfemq{code}.allagafgf",
            f"gf{code}.posagafgf", f"gfemq{code}.posagafgf",
        ])
        if fallback:
            trans_sources = [fallback]

    lines_sources = list(dict.fromkeys(g.lines for g in groups if g.lines))
    if not lines_sources:
        code_level = pick_first(files, [f"gf{code}.lines", f"gf{code}.lines-gz"])
        if code_level:
            lines_sources = [code_level]

    life = f"life{code}.dat" if f"life{code}.dat" in files else None
    pf = f"partfn{code}.dat" if f"partfn{code}.dat" in files else None

    reason = ""
    if not groups:
        reason = "no gf*.gam files found"
    return IonDiscovery(
        element, charge, code, name, url, True, groups, life, pf, reason,
        trans_sources=trans_sources, lines_sources=lines_sources,
    )


def discover_all_ions(
    timeout: int,
    base_url: str,
    offset: int = 0,
    limit: int | None = None,
    verbose: bool = False,
) -> list[IonDiscovery]:
    names = list_directory(f"{base_url.rstrip('/')}/", timeout)
    if names is None:
        raise RuntimeError(f"Could not list {base_url}")
    codes = sorted({
        match.group(1)
        for name in names
        for match in [re.match(r"^(\d{4})/?$", name)]
        if match
    })
    codes = codes[offset:]
    if limit is not None:
        codes = codes[:limit]
    discoveries = []
    for index, code in enumerate(codes, start=1):
        parsed = code_to_element_charge(code)
        if parsed is None:
            continue
        element, charge = parsed
        if verbose:
            print(f"discovering [{index}/{len(codes)}] {element}-{roman(charge + 1)} ({code})")
        discoveries.append(discover_ion(element, charge, timeout, base_url))
    return discoveries


def fetch_file(discovery: IonDiscovery, filename: str, timeout: int) -> str:
    url = urllib.parse.urljoin(discovery.url, filename)
    with urllib.request.urlopen(url, timeout=timeout) as response:
        data = response.read()
    if is_gz(filename):
        data = gzip.decompress(data)
    return data.decode("utf-8", errors="replace")


def fetch_optional_file(discovery: IonDiscovery, filename: str, timeout: int, role: str = "file") -> str | None:
    try:
        return fetch_file(discovery, filename, timeout)
    except urllib.error.HTTPError as exc:
        if exc.code in {403, 404}:
            print(f"  warning: optional {filename} unavailable ({exc.code}); continuing without it")
            if exc.code == 403:
                record_blocked_file(discovery, filename, role)
            return None
        raise


def extract_mapping_parts(gam_text: str, first_level_type: str) -> tuple[dict[str, str], dict[str, str]]:
    lines = gam_text.splitlines()
    # Mapping tables, when present, are the leading block before the first level
    # record.  Some files (notably Ni II/gf2801.gam) omit the textual "level type"
    # separator entirely.  Stopping only at that marker lets the parser consume the
    # complete level table as if it were a mapping table and corrupts encoded labels.
    end = next(
        (
            i for i, line in enumerate(lines)
            if (
                "level type" in line
                or line.lstrip().startswith("ELEM ")
                or LEVEL_RE.match(line[0:8].strip())
            )
        ),
        len(lines),
    )
    table_lines = lines[:end]
    start = next((i for i, line in enumerate(table_lines) if line.strip().startswith("1 ")), None)
    if start is None:
        return {}, {}

    table_lines = table_lines[start:]
    # The match table has two sections (even parity then odd parity), each starting
    # with a "1 ..." line as the key alphabet restarts. Split at the SECOND "1 "
    # line. (An earlier approach split on a "~ 0" marker, but that symbol row only
    # appears for large ions whose alphabet overflows past 'z' into symbols; small
    # ions like B II have no "~", so the odd section was lost and its match-table
    # letters stayed undecoded. The second "1 " line immediately follows the "~ 0"
    # row when present, so this is equivalent for large ions and also fixes small ones.)
    split = next(
        (i for i, line in enumerate(table_lines) if i > 0 and line.strip().startswith("1 ")),
        None,
    )
    if split is None:
        return {}, {}

    first_part = "\n".join(table_lines[:split])
    second_part = "\n".join(table_lines[split:])
    first_map = parse_mapping_table(first_part)
    second_map = parse_mapping_table(second_part)

    if first_level_type in {"ODD", "ORz", "OPo"}:
        return second_map, first_map
    return first_map, second_map


def parse_mapping_table(part: str) -> dict[str, str]:
    """Decode one parity section of a Kurucz match table.

    The table is a fixed-width grid: six ``MAPPING_CELL_WIDTH``-character cells per
    line, each holding a one-character key followed by its configuration. Splitting
    on whitespace instead fails in two ways.

    A configuration that fills its cell exactly (``d1 4s2 10p``) leaves no separator
    before the next cell, so the following key is swallowed into the value and lost
    as a key of its own -- levels using it then fall back to the raw label character.
    Slicing by column keeps neighbouring cells apart.

    A key can also appear twice in one section: once carrying a configuration, and
    again among the trailing rows that enumerate the still-unused characters. The
    assignment is authoritative, so a blank cell never overwrites one. Keys that are
    only ever blank are still recorded, because label decoding distinguishes a key
    that maps to nothing from one that is absent entirely.
    """
    mapping: dict[str, str] = {}
    for line in part.splitlines():
        for start in range(0, len(line), MAPPING_CELL_WIDTH):
            cell = line[start:start + MAPPING_CELL_WIDTH]
            key = cell[:1]
            if not key.strip():
                continue
            value = cell[1:].strip()
            if value or key not in mapping:
                mapping[key] = value
    return mapping


def gam_is_alternate_layout(gam_text: str) -> bool:
    """Detect the alternate Kurucz .gam layout used by a few ions (e.g. Cr VIII).

    The standard layout is ELEM, Index, E, J, label, g_lande; the alternate one
    swaps J and Index (ELEM, J, Index, E, ...) and ships no match table. They are
    told apart by the first data row's [8:12] field: in the standard layout this is
    the integer Index (no dot), in the alternate one it is the half/integer J value
    (always written with a decimal point). Returns False for any standard file, so
    the standard parsing path is never altered.
    """
    for line in gam_text.splitlines():
        if LEVEL_RE.match(line[0:8].strip()):
            return "." in line[8:12]
    return False


def declared_line_count(gam_text: str) -> int | None:
    """Return the transition count the .gam header claims for its own product.

    Kurucz writes the authoritative line total into the first record of every
    level file, in one of two layouts:

        " 164228  LINES SAVED  1191 EVEN LEVELS  1090 ODD LEVELS"
        "24.07    692274 lines saved  35 positive lines saved  1317 even ..."

    The leading count is the full list; the later "positive lines saved" figure
    is a small special-purpose subset, so only the first match is authoritative.
    Returns None when the header carries no such declaration.
    """
    header = gam_text.splitlines()[0] if gam_text else ""
    match = re.search(r"(\d+)\s+lines\s+saved", header, re.IGNORECASE)
    return int(match.group(1)) if match else None


def parse_gam(gam_text: str) -> pd.DataFrame:
    if gam_is_alternate_layout(gam_text):
        specs = [(0, 8), (8, 12), (12, 16), (16, 29), (29, 40), (40, 47)]
        names = ["ELEM", "J", "Index", "E", "label", "g_lande"]
    else:
        specs = [(0, 8), (8, 12), (12, 24), (24, 29), (29, 40), (40, 47)]
        names = ["ELEM", "Index", "E", "J", "label", "g_lande"]
    df = pd.read_fwf(StringIO(gam_text), colspecs=specs, header=None)
    df.columns = names
    df = df[["ELEM", "Index", "E", "J", "label", "g_lande"]]
    df = df[df["ELEM"].astype(str).str.match(LEVEL_RE, na=False)].copy()
    df["E"] = pd.to_numeric(df["E"], errors="coerce")
    df["J"] = pd.to_numeric(df["J"], errors="coerce")
    df["g_lande"] = pd.to_numeric(df["g_lande"], errors="coerce")
    # Avoid serialising IEEE negative zero as "-0.000000".  It is numerically
    # identical to zero and the sign carries no physical information here.
    df.loc[df["g_lande"] == 0, "g_lande"] = 0.0
    return df.dropna(subset=["E", "J", "label"])


def parse_life(life_text: str) -> pd.DataFrame:
    specs = [(0, 8), (8, 12), (12, 24), (24, 29), (29, 40), (40, 49), (49, 59), (59, 74)]
    df = pd.read_fwf(StringIO(life_text), colspecs=specs, header=None)
    df.columns = ["ELEM", "Index", "E", "J", "label", "SUM_A", "Life1", "Life(ns)"]
    df = df[df["ELEM"].astype(str).str.match(LEVEL_RE, na=False)].copy()
    df["E"] = pd.to_numeric(df["E"], errors="coerce")
    df["J"] = pd.to_numeric(df["J"], errors="coerce")
    df["Life(ns)"] = pd.to_numeric(df["Life(ns)"], errors="coerce")
    df = df.dropna(subset=["E", "J", "label"])
    df = df.drop_duplicates(subset=["E", "J", "label"])
    df["Life(s)"] = df["Life(ns)"] / 1e9
    return df[["E", "J", "label", "Life(s)"]]


def label_to_config_term(label: str, elem_type: str, mapping_even: dict[str, str], mapping_odd: dict[str, str]) -> tuple[str, str]:
    # A few fixed-width Kurucz labels lose the opening parenthesis at the start of
    # the 11-character label field (e.g. Cl I "3P)5s" and Co II "5D)4sp").  The
    # closing parenthesis and LS parent term make the intended form unambiguous.
    label = re.sub(r"^(\d+[A-Z]\))", r"(\1", str(label))
    parts = str(label).split()
    configuration = "unknown"
    term = "unknown"

    if len(parts) == 2:
        config, possible_term = parts
        if config.endswith("nd"):
            pass
        elif possible_term.isdigit():
            configuration = config[:-2]
            term = config[-2:]
            if configuration:
                if elem_type in ["EVE", "ERz", "EPo"]:
                    configuration = mapping_even.get(configuration[0], configuration[0]) + configuration[1:]
                elif elem_type in ["ODD", "ORz", "OPo"]:
                    configuration = mapping_odd.get(configuration[0], configuration[0]) + configuration[1:]
        else:
            configuration = config
            term = possible_term
    elif len(parts) == 3:
        config_1, config_2, possible_term = parts
        configuration = config_1 + config_2
        term = possible_term
    else:
        if len(str(label)) >= 4 and str(label)[-4].isupper():
            configuration = str(label)[:-5]
            term = str(label)[-5:-3]
            if configuration:
                if elem_type in ["EVE", "ERz", "EPo"]:
                    configuration = mapping_even.get(configuration[0], configuration[0]) + configuration[1:]
                elif elem_type in ["ODD", "ORz", "OPo"]:
                    configuration = mapping_odd.get(configuration[0], configuration[0]) + configuration[1:]
        elif "?" in str(label)[-4:-1]:
            configuration, term = str(label).split("?")[0], str(label).split("?")[1]

    if len(term) == 3 and re.match(r"^[a-zA-Z]\d[a-zA-Z]$", term):
        term = f"{term[0]}({term[1:]})"
    return configuration.replace(" ", ""), term


def has_lowercase_orbital_term(term: str) -> bool:
    """Reject malformed terms such as ``v(3p)``; orbital symbols are uppercase."""
    return re.search(r"\(\d+[spdfghiklmnoq]\)$", str(term)) is not None


def format_energy(value: float) -> str:
    """Render an energy into exactly 12 characters, trading decimals for digits.

    The .states layout gives E a 12-character field, but ``{:12.Nf}`` treats 12 as
    a minimum rather than a maximum, so an over-wide value silently pushes every
    later column right. Levels only reach seven integer digits above 1e6 cm^-1,
    which neutral and singly ionized species never do -- the shortfall first
    appears in the higher ionization stages this pipeline added.
    """
    int_part = len(str(abs(int(value))))
    decimals = 6 if int_part <= 5 else max(0, 11 - int_part)
    text = f"{value:12.{decimals}f}"
    # Rounding can carry into a new integer digit (999999.999999 -> 1000000.0).
    while len(text) > 12 and decimals > 0:
        decimals -= 1
        text = f"{value:12.{decimals}f}"
    return text


def write_states(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if df.empty:
        raise ValueError("no states to write")
    if float(df["J"].iloc[0]).is_integer():
        format_str = "{:>12d} {:>12} {:>6d} {:>7d} {:>12.6f} {:>12.4e} {:>10.6f} {:<12} {:<7} {:>2}\n"
    else:
        format_str = "{:>12d} {:>12} {:>6d} {:>7.1f} {:>12.6f} {:>12.4e} {:>10.6f} {:<12} {:<7} {:>2}\n"

    with path.open("w") as handle:
        for _, row in df.iterrows():
            j_value = int(row["J"]) if float(row["J"]).is_integer() else row["J"]
            handle.write(format_str.format(
                int(row["Index"]), format_energy(row["E"]), int(row["g_j"]), j_value,
                row["Uncertainty"], row["Life(s)"], row["g_lande"],
                row["Configuration"], row["Term"], row["Abbr"],
            ))


def build_states(
    raw_frames: list[pd.DataFrame],
    life: pd.DataFrame | None,
    ion_dir: Path,
    output_path: Path,
    save_intermediate: bool,
) -> pd.DataFrame:
    states = pd.concat(raw_frames, ignore_index=True)
    states = states.drop(columns=["Index"])
    states = states.drop_duplicates(subset=["ELEM", "E", "J", "label"])
    states["g_j"] = (2 * states["J"] + 1).astype(int)
    states["Uncertainty"] = 0.1

    if life is not None and not life.empty:
        combined = states.merge(life, on=["E", "J", "label"], how="left")
    else:
        combined = states.copy()
        combined["Life(s)"] = math.nan

    combined["Abbr"] = combined["E"].apply(lambda value: "CA" if value < 0 else "NI")
    combined["E"] = combined["E"].abs()

    configs = []
    terms = []
    for _, row in combined.iterrows():
        config, term = label_to_config_term(row["label"], str(row["ELEM"])[-3:], row["mapping_even"], row["mapping_odd"])
        configs.append(config)
        terms.append(term)
    combined["Configuration"] = configs
    combined["Term"] = terms
    combined = combined[(combined["Configuration"] != "unknown") & (combined["Term"] != "unknown")]
    malformed_terms = combined["Term"].map(has_lowercase_orbital_term)
    if malformed_terms.any():
        print(f"    removed {int(malformed_terms.sum())} state(s) with lowercase orbital term symbols")
        combined = combined[~malformed_terms]
    combined = combined.sort_values(by="E").reset_index(drop=True)
    if not combined.empty:
        combined.loc[0, "Life(s)"] = float("inf")

    combined.insert(0, "Index", range(1, len(combined) + 1))

    states_to_trans = combined.drop(columns=["ELEM", "Configuration", "Term", "mapping_even", "mapping_odd"])
    if save_intermediate:
        ion_dir.mkdir(parents=True, exist_ok=True)
        states_to_trans.to_csv(ion_dir / "States_Final.csv", index=False)

    out = combined.drop(columns=["label", "ELEM", "mapping_even", "mapping_odd"])
    out = out[["Index", "E", "g_j", "J", "Uncertainty", "Life(s)", "g_lande", "Configuration", "Term", "Abbr"]]
    write_states(out, output_path)
    return states_to_trans


def parse_agafgf_transition(line: str) -> tuple[float, float, str, float, float, str, float, float] | None:
    """Parse one .agafgf / .allagafgf row.

    The Kurucz agafgf files carry, in addition to the wavenumber and log A, the
    full transition endpoints (E, J, label for both states) in fixed columns (the
    extended layout from the paper's Note 2). Because every value needed for a
    .trans row lives in this single file, transitions are built from it alone and
    no positional join with the .lines file is required -- which is what makes the
    output immune to the .lines/.agafgf row-misalignment problem.

    Returns (E1, J1, label1, E2, J2, label2, wn, A) or None for header/garbage.
    """
    try:
        wn = float(line[12:24])
        log_a = float(line[45:52])
        e1 = float(line[65:77])
        j1 = float(line[77:82])
        label1 = line[82:93].strip()
        e2 = float(line[93:106])
        j2 = float(line[106:111])
        label2 = line[111:121].strip()
    except (ValueError, IndexError):
        return None
    if not label1 or not label2:
        return None
    return e1, j1, label1, e2, j2, label2, abs(wn), 10 ** log_a


def parse_line_transition(line: str) -> tuple[float, float, str, float, float, str] | None:
    """Parse the two transition endpoints from one standard Kurucz .lines row."""
    try:
        e1 = float(line[24:36])
        j1 = float(line[36:41])
        label1 = line[41:52].strip()
        e2 = float(line[52:64])
        j2 = float(line[64:69])
        label2 = line[69:80].strip()
    except (ValueError, IndexError):
        return None
    if not label1 or not label2:
        return None
    return e1, j1, label1, e2, j2, label2


def parse_agafgf_value(line: str) -> tuple[float, float] | None:
    """Parse wavenumber and Einstein A from one standard Kurucz .agafgf row."""
    try:
        wn = float(line[12:24])
        log_a = float(line[45:52])
    except (ValueError, IndexError):
        return None
    return abs(wn), 10 ** log_a


def iter_parsed_lines(lines):
    for line in lines:
        record = parse_line_transition(line)
        if record is not None:
            yield record


def iter_parsed_agafgf_values(lines):
    for line in lines:
        record = parse_agafgf_value(line)
        if record is not None:
            yield record


def transition_source_pairs(discovery: IonDiscovery) -> list[tuple[str, str | None]]:
    """Return each AGAFGF source with its corresponding .lines source, if any."""
    line_by_agafgf = {
        group.agafgf: group.lines
        for group in discovery.groups
        if group.agafgf
    }
    return [(source, line_by_agafgf.get(source)) for source in discovery.trans_sources]


def iter_agafgf_records(text: str):
    for line in text.splitlines():
        record = parse_agafgf_transition(line)
        if record is not None:
            yield record


# A_ul = 6.6702e15 * gf / (g_u * lambda_A^2); with lambda_A = 1e8 / wn[cm-1] the
# constant folds to 6.6702e15 / 1e16.
EINSTEIN_A_FROM_GF = 0.66702


def parse_lines_transition(line: str) -> tuple[float, float, str, float, float, str, float, float] | None:
    """Parse one Kurucz .lines row into the same 8-tuple the agafgf parser yields.

    The .lines files carry both endpoints and log(gf) but no Einstein A, so A is
    derived from log(gf) with the standard relation above. The wavenumber comes
    from the endpoint energies themselves (|E| because negative marks a predicted
    level), which is what the .agafgf wavenumber column also holds.

    Returns (E1, J1, label1, E2, J2, label2, wn, A) or None for header/garbage.
    """
    endpoints = parse_line_transition(line)
    if endpoints is None:
        return None
    e1, j1, label1, e2, j2, label2 = endpoints
    try:
        log_gf = float(line[11:18])
    except (ValueError, IndexError):
        return None
    wn = abs(abs(e2) - abs(e1))
    if wn <= 0:
        return None
    # The upper state carries the degeneracy in the A relation; |E| decides which
    # endpoint that is, never the column order.
    j_upper = j2 if abs(e2) > abs(e1) else j1
    g_upper = 2.0 * j_upper + 1.0
    if g_upper <= 0:
        return None
    return e1, j1, label1, e2, j2, label2, wn, EINSTEIN_A_FROM_GF * (10.0 ** log_gf) * wn * wn / g_upper


def transition_pairs_in(path: Path) -> set[tuple[int, int]]:
    """Read back the (upper, lower) pairs already written to a .trans file."""
    pairs: set[tuple[int, int]] = set()
    with path.open() as handle:
        for row in handle:
            parts = row.split()
            if len(parts) >= 2:
                pairs.add((int(parts[0]), int(parts[1])))
    return pairs


def map_records_to_rows(
    records, index: dict, seen_pairs: set | None = None
) -> tuple[list[tuple], int, int]:
    """Resolve raw transition records to (upper, lower, A, wn) output rows.

    When ``seen_pairs`` is given, records resolving to a (upper, lower) pair already
    in it are skipped and every accepted pair is added, so a second source can top
    up a first one without duplicating the lines they share.
    """
    rows = []
    dropped = 0
    misaligned = 0
    for e1, j1, label1, e2, j2, label2, wn, a_value in records:
        if not endpoints_match_wavenumber(e1, e2, wn):
            misaligned += 1
            continue
        index1 = index.get(state_key(e1, j1, label1))
        index2 = index.get(state_key(e2, j2, label2))
        if index1 is None or index2 is None:
            dropped += 1
            continue
        upper, lower = order_transition(index1, index2, e1, e2)
        if seen_pairs is not None:
            if (upper, lower) in seen_pairs:
                continue
            seen_pairs.add((upper, lower))
        rows.append((upper, lower, a_value, wn))
    return rows, dropped, misaligned


def build_transitions(records: list[tuple], states_to_trans: pd.DataFrame, output_path: Path) -> tuple[int, int]:
    index = build_state_index(states_to_trans)
    rows, dropped, misaligned = map_records_to_rows(records, index)
    if misaligned:
        dropped += misaligned
        print(
            f"    warning: {misaligned} rows dropped: agafgf wavenumber disagrees with "
            f"the endpoint energies (misaligned pairing or inconsistent source rows)"
        )

    if not rows:
        return 0, dropped

    # No de-duplication: the agafgf file is the authoritative line list and the
    # streaming path keeps every row, so the pandas path must match it (it only
    # adds a wavenumber sort). Removing duplicates here would drop genuine lines
    # and disagree with both the streaming output and the published counts.
    result = pd.DataFrame(rows, columns=["Index1", "Index2", "A", "wn"])
    result = result.sort_values(by="wn")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as handle:
        for _, row in result.iterrows():
            handle.write("{:>12d}{:>1}{:>12d}{:>1}{:>10.4e}{:>1}{:>15.6e}\n".format(
                int(row["Index1"]), "", int(row["Index2"]), "", row["A"], "", row["wn"]
            ))
    return len(result), dropped


def build_transitions_supplemented(
    discovery: "IonDiscovery",
    records: list[tuple],
    states_to_trans: pd.DataFrame,
    output_path: Path,
    declared_total: int,
    timeout: int,
) -> tuple[int, int]:
    """Write transitions from the agafgf records, topping up from .lines if short.

    Kurucz ships several products per ion and only some carry an .agafgf. When the
    resolved one is a small subset (typically a gfemq POS/METAPOS export, which has
    the predicted lines deleted), the agafgf rows are kept -- they are the
    authoritative log A values and the only source of the forbidden M1/E2 lines --
    and the missing E1 lines are added from the .lines file, whose A is derived
    from log(gf). Sources are merged on the endpoint pair, agafgf winning.
    """
    index = build_state_index(states_to_trans)
    # Tracking pairs costs as much memory as the records themselves, so only do it
    # when the agafgf is already too small to be the whole line list.
    shortfall = bool(declared_total) and len(records) < TRANSITION_COMPLETENESS_THRESHOLD * declared_total
    seen_pairs: set | None = set() if shortfall else None
    rows, dropped, misaligned = map_records_to_rows(records, index, seen_pairs)
    if misaligned:
        dropped += misaligned
        print(
            f"    warning: {misaligned} rows dropped: agafgf wavenumber disagrees with "
            f"the endpoint energies (misaligned pairing or inconsistent source rows)"
        )

    if shortfall and discovery.lines_sources:
        print(
            f"    transitions are short of the {declared_total} declared by the .gam "
            f"header; supplementing from {', '.join(discovery.lines_sources)}"
        )
        for lines_source in discovery.lines_sources:
            try:
                stream = iter_response_lines(discovery, lines_source, timeout)
                supplement = (
                    record for line in stream
                    if (record := parse_lines_transition(line)) is not None
                )
                added, added_dropped, _ = map_records_to_rows(supplement, index, seen_pairs)
            except urllib.error.HTTPError as exc:
                if exc.code in {403, 404}:
                    print(f"    warning: {lines_source} unavailable ({exc.code}); skipping")
                    if exc.code == 403:
                        record_blocked_file(discovery, lines_source, "lines")
                    continue
                raise
            print(f"    {lines_source}: added={len(added)} dropped={added_dropped}")
            rows.extend(added)
            dropped += added_dropped
    elif shortfall:
        print("    warning: transitions are short of the declared count and no .lines source is available")

    if not rows:
        return 0, dropped
    result = pd.DataFrame(rows, columns=["Index1", "Index2", "A", "wn"]).sort_values(by="wn")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as handle:
        for upper, lower, a_value, wn in result.itertuples(index=False):
            handle.write("{:>12d}{:>1}{:>12d}{:>1}{:>10.4e}{:>1}{:>15.6e}\n".format(
                int(upper), "", int(lower), "", a_value, "", wn))
    return len(result), dropped


def endpoints_match_wavenumber(e1: float, e2: float, wn: float, tolerance: float = 0.5) -> bool:
    """True when the wavenumber agrees with |E2 - E1| (both files rounded to
    ~0.05 cm-1). A violation means the paired .lines/.agafgf rows describe two
    different transitions, i.e. the files are not row-aligned."""
    return abs(abs(abs(e2) - abs(e1)) - wn) <= tolerance


def order_transition(index1: int, index2: int, e1: float, e2: float) -> tuple[int, int]:
    """Return (upper, lower) state indices for a .trans row.

    The ExoMol format requires the upper state in the first column, but Kurucz
    raw rows carry their endpoints in arbitrary order. Compare on |E| because
    negative energies mark predicted levels.
    """
    if abs(e2) > abs(e1):
        return index2, index1
    return index1, index2


def state_key(energy: float, j_value: float, label: str) -> tuple[float, float, str]:
    return (round(abs(float(energy)), 6), round(float(j_value), 3), str(label).strip())


class StateIndex:
    """(E, J, label) -> state Index lookup with an energy-tolerance fallback.

    Kurucz regenerates the agafgf line lists more often than the gam level
    files, shifting level energies by ~0.001-0.1 cm-1 between generations
    (e.g. Mn-II). An exact-energy dictionary would drop every transition that
    references a recalibrated level, so a miss falls back to the nearest level
    with the same J and label within TOLERANCE.
    """

    TOLERANCE = 0.1

    def __init__(self, states: pd.DataFrame) -> None:
        self.exact: dict[tuple[float, float, str], int] = {}
        by_j: dict[float, list[tuple[float, str, int]]] = {}
        for _, row in states.iterrows():
            key = state_key(row["E"], row["J"], row["label"])
            index = int(row["Index"])
            self.exact.setdefault(key, index)
            by_j.setdefault(key[1], []).append((key[0], key[2], index))
        self.by_j = {j: sorted(items) for j, items in by_j.items()}
        self.energies_by_j = {j: [energy for energy, _, _ in items] for j, items in self.by_j.items()}

    def get(self, key: tuple[float, float, str]) -> int | None:
        hit = self.exact.get(key)
        if hit is not None:
            return hit
        items = self.by_j.get(key[1])
        if not items:
            return None
        energies = self.energies_by_j[key[1]]
        low = bisect.bisect_left(energies, key[0] - self.TOLERANCE)
        high = bisect.bisect_right(energies, key[0] + self.TOLERANCE)
        # Prefer a same-label candidate; the newer line lists relabel levels
        # whose eigenvector composition changed, so a differing label with an
        # essentially identical energy is still the same physical level.
        best = None
        for energy, label, index in items[low:high]:
            rank = (label != key[2], abs(energy - key[0]))
            if best is None or rank < best[0]:
                best = (rank, index)
        return best[1] if best else None


def build_state_index(states_to_trans: pd.DataFrame) -> StateIndex:
    states = states_to_trans[["Index", "E", "J", "label"]].drop_duplicates(subset=["E", "J", "label"])
    return StateIndex(states)


def iter_response_lines(discovery: IonDiscovery, filename: str, timeout: int):
    url = urllib.parse.urljoin(discovery.url, filename)
    with urllib.request.urlopen(url, timeout=timeout) as response:
        if is_gz(filename):
            with gzip.GzipFile(fileobj=response) as stream:
                for raw_line in stream:
                    yield raw_line.decode("utf-8", errors="replace")
        else:
            for raw_line in response:
                yield raw_line.decode("utf-8", errors="replace")


def build_transitions_stream(
    discovery: IonDiscovery,
    sources: list[tuple[str, str | None]],
    states_to_trans: pd.DataFrame,
    output_path: Path,
    timeout: int,
    declared_total: int = 0,
) -> tuple[int, int]:
    index = build_state_index(states_to_trans)
    written = 0
    dropped = 0
    output_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = output_path.with_name(f".{output_path.name}.tmp")

    try:
        with temporary_path.open("w") as handle:
            for source, lines_source in sources:
                try:
                    source_written = 0
                    source_parsed = 0
                    source_dropped = 0
                    source_misaligned = 0
                    # Probe the agafgf head: when it carries the extended endpoint
                    # columns it is parsed standalone. Positional pairing with a
                    # .lines file is kept only for standard-layout agafgf files --
                    # the two site files can fall out of sync (e.g. Mn-II, where
                    # gf2501.agafgf gained ~10k rows that gf2501.lines lacks), and
                    # a desynced zip corrupts every row after the first insertion.
                    source_lines = iter_response_lines(discovery, source, timeout)
                    head = list(itertools.islice(source_lines, 200))
                    source_lines = itertools.chain(head, source_lines)
                    extended = any(parse_agafgf_transition(line) is not None for line in head)
                    if extended:
                        if lines_source:
                            print(f"    note: {source} carries endpoint columns; parsing standalone (ignoring {lines_source})")
                        paired_records = (
                            (record[:6], record[6:])
                            for line in source_lines
                            if (record := parse_agafgf_transition(line)) is not None
                        )
                    elif lines_source:
                        endpoint_records = iter_parsed_lines(
                            iter_response_lines(discovery, lines_source, timeout)
                        )
                        value_records = iter_parsed_agafgf_values(source_lines)
                        paired_records = itertools.zip_longest(endpoint_records, value_records)
                    else:
                        paired_records = iter(())

                    for endpoints, values in paired_records:
                        if endpoints is None or values is None:
                            source_dropped += 1
                            continue
                        source_parsed += 1
                        e1, j1, label1, e2, j2, label2 = endpoints
                        wn, a_value = values
                        if not endpoints_match_wavenumber(e1, e2, wn):
                            source_misaligned += 1
                            continue
                        index1 = index.get(state_key(e1, j1, label1))
                        index2 = index.get(state_key(e2, j2, label2))
                        if index1 is None or index2 is None:
                            source_dropped += 1
                            continue
                        upper, lower = order_transition(index1, index2, e1, e2)
                        handle.write("{:>12d}{:>1}{:>12d}{:>1}{:>10.4e}{:>1}{:>15.6e}\n".format(
                            upper, "", lower, "", a_value, "", wn
                        ))
                        source_written += 1
                    written += source_written
                    dropped += source_dropped + source_misaligned
                    partner = "extended endpoints" if extended else (lines_source or "no endpoint source")
                    print(
                        f"    {source} + {partner}: parsed={source_parsed} "
                        f"written={source_written} dropped={source_dropped}"
                    )
                    if source_misaligned:
                        print(
                            f"    warning: {source_misaligned} rows dropped: agafgf wavenumber disagrees with "
                            f"the endpoint energies (misaligned pairing or inconsistent source rows)"
                        )
                except urllib.error.HTTPError as exc:
                    if exc.code in {403, 404}:
                        print(f"    warning: transition source {source} unavailable ({exc.code}); skipping")
                        if exc.code == 403:
                            record_blocked_file(discovery, source, "agafgf")
                        continue
                    raise

            # Kurucz ships several products per ion and only some carry an .agafgf.
            # When the resolved one is a small subset (typically a gfemq POS/METAPOS
            # export, which has the predicted lines deleted), keep its rows -- they
            # are the authoritative log A values and the only source of the forbidden
            # M1/E2 lines -- and add the missing E1 lines from .lines, whose A is
            # derived from log(gf). The pairs are read back from what was just
            # written, which is cheap precisely because the file is short.
            if declared_total and written + dropped < TRANSITION_COMPLETENESS_THRESHOLD * declared_total:
                if not discovery.lines_sources:
                    print("    warning: transitions are short of the declared count and no .lines source is available")
                else:
                    print(
                        f"    transitions are short of the {declared_total} declared by the "
                        f".gam header; supplementing from {', '.join(discovery.lines_sources)}"
                    )
                    handle.flush()
                    seen_pairs = transition_pairs_in(temporary_path)
                    for lines_source in discovery.lines_sources:
                        try:
                            stream = iter_response_lines(discovery, lines_source, timeout)
                            added = added_dropped = 0
                            for line in stream:
                                record = parse_lines_transition(line)
                                if record is None:
                                    continue
                                rows, row_dropped, _ = map_records_to_rows([record], index, seen_pairs)
                                added_dropped += row_dropped
                                for upper, lower, a_value, wn in rows:
                                    handle.write(
                                        "{:>12d}{:>1}{:>12d}{:>1}{:>10.4e}{:>1}{:>15.6e}\n".format(
                                            upper, "", lower, "", a_value, "", wn))
                                    added += 1
                        except urllib.error.HTTPError as exc:
                            if exc.code in {403, 404}:
                                print(f"    warning: {lines_source} unavailable ({exc.code}); skipping")
                                if exc.code == 403:
                                    record_blocked_file(discovery, lines_source, "lines")
                                continue
                            raise
                        print(f"    {lines_source}: added={added} dropped={added_dropped}")
                        written += added
                        dropped += added_dropped

        if written:
            temporary_path.replace(output_path)
        else:
            temporary_path.unlink(missing_ok=True)
            print("    warning: no valid transitions were produced; existing output was preserved")
    except BaseException:
        temporary_path.unlink(missing_ok=True)
        raise
    return written, dropped


def parse_pf(pf_text: str) -> pd.DataFrame:
    rows = []
    for line in pf_text.splitlines():
        parts = line.split()
        if len(parts) >= 4:
            try:
                float(parts[0])
                temp = float(parts[2])
                value = float(parts[3])
            except ValueError:
                continue
            rows.append((temp, value))
    return pd.DataFrame(rows, columns=["T", "Value"])


def write_pf(pf: pd.DataFrame, csv_path: Path | None, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if csv_path is not None:
        csv_path.parent.mkdir(parents=True, exist_ok=True)
        pf.to_csv(csv_path, index=False)
    with output_path.open("w") as handle:
        for _, row in pf.iterrows():
            handle.write("{:>8.1f}{:>1}{:>15.4f}\n".format(row["T"], "", row["Value"]))


def ion_paths(args: argparse.Namespace, discovery: IonDiscovery) -> tuple[Path, Path, Path, Path]:
    ion_root = Path(args.data_root) / discovery.ion_name
    raw_dir = ion_root / "raw"
    intermediate_dir = ion_root / "intermediate"
    exomol_dir = ion_root / "exomol"
    return ion_root, raw_dir, intermediate_dir, exomol_dir


def exomol_name(discovery: IonDiscovery, suffix: str) -> str:
    """Return an ExoMol filename such as ``Zr_III__Kurucz.states``."""
    species = discovery.ion_name.replace("-", "_", 1)
    return f"{species}__Kurucz.{suffix}"


def expected_outputs(args: argparse.Namespace, discovery: IonDiscovery) -> list[Path]:
    _, _, _, exomol_dir = ion_paths(args, discovery)
    outputs = [exomol_dir / exomol_name(discovery, "states")]
    if discovery.trans_sources:
        outputs.append(exomol_dir / exomol_name(discovery, "trans"))
    if discovery.pf:
        outputs.append(exomol_dir / exomol_name(discovery, "pf"))
    return outputs


def is_complete(args: argparse.Namespace, discovery: IonDiscovery) -> bool:
    outputs = expected_outputs(args, discovery)
    return bool(outputs) and all(path.exists() and path.stat().st_size > 0 for path in outputs)


def save_raw(raw_dir: Path, filename: str, text: str) -> None:
    raw_dir.mkdir(parents=True, exist_ok=True)
    (raw_dir / filename).write_text(text)


def process_ion(discovery: IonDiscovery, args: argparse.Namespace) -> None:
    print(f"\n{discovery.ion_name} ({discovery.code})")
    if not discovery.exists:
        print(f"  skipped: {discovery.reason}")
        return
    if not discovery.groups:
        print(f"  skipped: {discovery.reason}")
        return
    if not args.overwrite and is_complete(args, discovery):
        print("  skipped: outputs already exist; use --overwrite to regenerate")
        return

    _, raw_dir, intermediate_dir, exomol_dir = ion_paths(args, discovery)
    if args.save_raw:
        raw_dir.mkdir(parents=True, exist_ok=True)
    if args.save_intermediate:
        intermediate_dir.mkdir(parents=True, exist_ok=True)
    exomol_dir.mkdir(parents=True, exist_ok=True)
    raw_state_frames = []
    declared_lines: dict[str, int] = {}

    life = None
    if discovery.life:
        print(f"  downloading {discovery.life}")
        life_text = fetch_optional_file(discovery, discovery.life, args.timeout, "life")
        if life_text is not None:
            if args.save_raw:
                save_raw(raw_dir, discovery.life, life_text)
            life = parse_life(life_text)
            if args.save_intermediate:
                life.to_csv(intermediate_dir / "LIFE.csv", index=False, header=False)
    else:
        print("  no lifetime file; Life(s) will be NaN except ground state")

    for group in discovery.groups:
        print(f"  group {group.stem}")
        try:
            gam_text = fetch_file(discovery, group.gam, args.timeout)
        except urllib.error.HTTPError as exc:
            if exc.code in {403, 404}:
                if exc.code == 403:
                    record_blocked_file(discovery, group.gam, "gam")
                print(f"    {group.gam} unavailable ({exc.code}); skipping group")
                continue
            raise
        if args.save_raw:
            save_raw(raw_dir, group.gam, gam_text)
        states = parse_gam(gam_text)
        if states.empty:
            print(f"    no state rows in {group.gam}; skipped")
            continue
        first_type = str(states["ELEM"].iloc[0])[-3:]
        mapping_even, mapping_odd = extract_mapping_parts(gam_text, first_type)
        states["mapping_even"] = [mapping_even] * len(states)
        states["mapping_odd"] = [mapping_odd] * len(states)
        if args.save_intermediate:
            states.drop(columns=["mapping_even", "mapping_odd"]).to_csv(intermediate_dir / f"{group.stem}_GAM.csv", index=False, header=False)
        raw_state_frames.append(states)
        declared = declared_line_count(gam_text)
        if declared is not None:
            declared_lines[group.stem] = declared
        print(f"    states={len(states)}" + (f", {declared} lines declared" if declared else ""))

    if not raw_state_frames:
        print("  skipped: no usable states")
        return

    states_path = exomol_dir / exomol_name(discovery, "states")
    states_to_trans = build_states(raw_state_frames, life, intermediate_dir, states_path, args.save_intermediate)
    print(f"  wrote {states_path} ({len(states_to_trans)} states)")

    if args.states_only:
        print("  --states-only: leaving existing .trans/.pf untouched")
        return

    if discovery.trans_sources:
        trans_path = exomol_dir / exomol_name(discovery, "trans")
        print(f"  transitions from: {', '.join(discovery.trans_sources)}")
        source_pairs = transition_source_pairs(discovery)
        if args.stream_transitions:
            trans_rows, dropped = build_transitions_stream(
                discovery, source_pairs, states_to_trans, trans_path, args.timeout,
                sum(declared_lines.values()),
            )
            if trans_rows:
                print(f"  wrote {trans_path} ({trans_rows} transitions, dropped {dropped}, streamed)")
            else:
                print(f"  no transition output written ({dropped} records dropped)")
        else:
            records = []
            for source, lines_source in source_pairs:
                agafgf_text = fetch_optional_file(discovery, source, args.timeout, "agafgf")
                if agafgf_text is None:
                    continue
                if args.save_raw:
                    save_raw(raw_dir, source, agafgf_text)
                agafgf_lines = agafgf_text.splitlines()
                # Same extended-layout probe as the streaming path: standalone
                # parsing whenever the agafgf carries its own endpoints, because
                # positional pairing corrupts every row once the two site files
                # fall out of sync.
                extended = any(parse_agafgf_transition(line) is not None for line in agafgf_lines[:200])
                if extended and lines_source:
                    print(f"    note: {source} carries endpoint columns; parsing standalone (ignoring {lines_source})")
                if not extended and lines_source:
                    lines_text = fetch_optional_file(discovery, lines_source, args.timeout, "lines")
                    if lines_text is None:
                        continue
                    if args.save_raw:
                        save_raw(raw_dir, lines_source, lines_text)
                    endpoints = iter_parsed_lines(lines_text.splitlines())
                    values = iter_parsed_agafgf_values(agafgf_lines)
                    for endpoint, value in itertools.zip_longest(endpoints, values):
                        if endpoint is None or value is None:
                            continue
                        records.append((*endpoint, *value))
                else:
                    records.extend(iter_agafgf_records(agafgf_text))
            trans_rows, dropped = build_transitions_supplemented(
                discovery, records, states_to_trans, trans_path,
                sum(declared_lines.values()), args.timeout,
            )
            if trans_rows:
                print(f"  wrote {trans_path} ({trans_rows} transitions, dropped {dropped})")
            else:
                print(f"  no transition output written ({dropped} records dropped)")
        check_transition_completeness(discovery, declared_lines, trans_rows, dropped)
    else:
        print("  no transition source; transitions skipped")

    if discovery.pf:
        print(f"  downloading {discovery.pf}")
        pf_text = fetch_optional_file(discovery, discovery.pf, args.timeout, "pf")
        if pf_text is not None:
            if args.save_raw:
                save_raw(raw_dir, discovery.pf, pf_text)
            pf = parse_pf(pf_text)
            if not pf.empty:
                pf_path = exomol_dir / exomol_name(discovery, "pf")
                pf_csv_path = intermediate_dir / "PF.csv" if args.save_intermediate else None
                write_pf(pf, pf_csv_path, pf_path)
                print(f"  wrote {pf_path} ({len(pf)} rows)")
            else:
                print("  partition function file had no parseable rows")
    else:
        print("  no partition function file")


def print_discovery(discovery: IonDiscovery) -> None:
    print(f"{discovery.ion_name:8} {discovery.code} ", end="")
    if not discovery.exists:
        print(f"missing ({discovery.reason})")
        return
    if not discovery.groups:
        print(f"exists, skipped ({discovery.reason})")
        return
    groups = ", ".join(
        f"{group.stem}[gam={bool(group.gam)},lines={bool(group.lines)},agafgf={bool(group.agafgf)}]"
        for group in discovery.groups
    )
    print(f"life={bool(discovery.life)} pf={bool(discovery.pf)} groups={groups}")


def write_manifest(discoveries: list[IonDiscovery], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "element", "charge", "ion", "code", "url", "exists", "processable",
            "life", "pf", "groups", "reason",
        ])
        for discovery in discoveries:
            groups = ";".join(
                f"{group.stem}:gam={bool(group.gam)},lines={bool(group.lines)},agafgf={bool(group.agafgf)}"
                for group in discovery.groups
            )
            writer.writerow([
                discovery.element,
                discovery.charge,
                discovery.ion_name,
                discovery.code,
                discovery.url,
                discovery.exists,
                bool(discovery.groups),
                bool(discovery.life),
                bool(discovery.pf),
                groups,
                discovery.reason,
            ])


def write_blocked(path: Path) -> int:
    """Write the accumulated 403 (Forbidden) files to CSV, merged with any rows
    from a previous run so batched runs (--offset/--limit) keep accumulating into
    one table. De-duplicated by (code, filename, role). Returns the total count."""
    headers = ["element", "charge", "code", "ion", "role", "filename", "url", "http_status"]
    merged: dict[tuple[str, str, str], dict[str, str]] = {}
    if path.exists():
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                merged[(row["code"], row["filename"], row["role"])] = row
    for entry in BLOCKED:
        row = {key: str(value) for key, value in entry.items()}
        merged[(row["code"], row["filename"], row["role"])] = row

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers)
        writer.writeheader()
        writer.writerows(sorted(merged.values(), key=lambda r: (r["code"], r["filename"])))
    return len(merged)


def read_manifest(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    with path.open(newline="") as handle:
        return {row["ion"]: row for row in csv.DictReader(handle)}


def truthy(value: str | bool | None) -> bool:
    return str(value).lower() == "true"


def validate_outputs(data_root: Path, manifest_path: Path, report_dir: Path) -> None:
    manifest = read_manifest(manifest_path)
    report_dir.mkdir(parents=True, exist_ok=True)

    rows_by_category = {
        "complete": [],
        "missing_trans": [],
        "missing_pf": [],
        "missing_states": [],
        "failed_empty": [],
        "unusual": [],
    }

    report_dir_resolved = report_dir.resolve()
    for ion_dir in sorted(path for path in data_root.iterdir() if path.is_dir()):
        if ion_dir.resolve() == report_dir_resolved:
            continue
        ion = ion_dir.name
        if manifest and ion not in manifest:
            continue
        exomol = ion_dir / "exomol"
        files = sorted(path.name for path in exomol.iterdir() if path.is_file()) if exomol.exists() else []
        have = {Path(name).suffix.lstrip(".") for name in files}
        manifest_row = manifest.get(ion, {})

        source_has_trans = "lines=True,agafgf=True" in manifest_row.get("groups", "")
        source_has_pf = truthy(manifest_row.get("pf"))
        source_has_states = truthy(manifest_row.get("processable"))

        missing = sorted({"states", "trans", "pf"} - have)
        reasons = []
        if "states" in missing:
            reasons.append("states_missing_after_processing" if source_has_states else "source_not_processable")
        if "trans" in missing:
            reasons.append("source_missing_lines_or_agafgf" if not source_has_trans else "trans_missing_after_processing")
        if "pf" in missing:
            reasons.append("source_missing_partfn" if not source_has_pf else "pf_missing_after_processing")

        row = {
            "ion": ion,
            "code": manifest_row.get("code", ""),
            "have": ",".join(sorted(have)),
            "missing": ",".join(missing),
            "reason": ";".join(reasons),
            "files": ",".join(files),
            "source_groups": manifest_row.get("groups", ""),
            "source_pf": manifest_row.get("pf", ""),
        }

        if not files:
            rows_by_category["failed_empty"].append(row)
        elif not missing:
            rows_by_category["complete"].append(row)
        else:
            if "states" in missing:
                rows_by_category["missing_states"].append(row)
            if "trans" in missing:
                rows_by_category["missing_trans"].append(row)
            if "pf" in missing:
                rows_by_category["missing_pf"].append(row)
            if len(have) not in {1, 2, 3}:
                rows_by_category["unusual"].append(row)

    headers = ["ion", "code", "have", "missing", "reason", "files", "source_groups", "source_pf"]
    for category, rows in rows_by_category.items():
        with (report_dir / f"{category}.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=headers)
            writer.writeheader()
            writer.writerows(rows)

    print(f"Validation reports written to {report_dir}")
    for category, rows in rows_by_category.items():
        print(f"  {category}: {len(rows)}")


def process_discoveries(discoveries: list[IonDiscovery], args: argparse.Namespace) -> None:
    for index, discovery in enumerate(discoveries, start=1):
        print(f"\n[{index}/{len(discoveries)}]", end="")
        try:
            process_ion(discovery, args)
        except Exception as exc:
            print(f"\n  failed: {type(exc).__name__}: {exc}")
            if not args.keep_going:
                raise
        if args.delay and index < len(discoveries):
            time.sleep(args.delay)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Discover and process Kurucz atom/ion data into ExoMol-style files.")
    parser.add_argument("--element", help="Element symbol, for example Li, N, Fe. Omit when using --all.")
    parser.add_argument("--all", action="store_true", help="Discover every atom/ion directory under the Kurucz atoms index.")
    parser.add_argument("--validate", action="store_true", help="Validate local exomol outputs and write CSV reports. Does not access the network.")
    parser.add_argument("--charge", type=int, help="Process one charge state: 0=I, 1=II, 2=III.")
    parser.add_argument("--start-charge", type=int, default=0, help="First charge state when using --max-charge.")
    parser.add_argument("--max-charge", type=int, default=1, help="Last charge state when --charge is not set.")
    parser.add_argument("--list-only", action="store_true", help="Only discover available files; do not download/process data.")
    parser.add_argument("--base-url", default=BASE_URL)
    parser.add_argument("--data-root", default="Kurucz-data", help="Root folder for raw, intermediate, and final ion data.")
    parser.add_argument("--manifest", help="Write discovery results to this CSV file.")
    parser.add_argument("--manifest-path", default="reports/kurucz-discovery/discovery.csv", help="Manifest CSV used by --validate.")
    parser.add_argument("--report-dir", default="reports/kurucz-validation", help="Directory for --validate CSV reports.")
    parser.add_argument("--blocked-report", default="reports/kurucz-blocked/blocked-403.csv", help="CSV listing files the server returns 403 (Forbidden) for, to request access. Merged across runs.")
    parser.add_argument("--incomplete-report", default="reports/kurucz-incomplete/incomplete-transitions.csv", help="CSV listing ions whose transition source supplied fewer records than the .gam header declares. Merged across runs.")
    parser.add_argument("--offset", type=int, default=0, help="Skip this many atom/ion directories in --all mode.")
    parser.add_argument("--limit", type=int, help="Limit the number of atom/ion directories discovered in --all mode.")
    parser.add_argument("--no-save-raw", dest="save_raw", action="store_false", help="Do not keep downloaded raw Kurucz files on disk.")
    parser.add_argument("--no-save-intermediate", dest="save_intermediate", action="store_false", help="Do not keep parsed CSV intermediate files on disk.")
    parser.add_argument("--pandas-transitions", dest="stream_transitions", action="store_false", help="Use the old pandas transition builder instead of streaming. This needs much more RAM but sorts transitions by wavenumber.")
    parser.add_argument("--states-only", action="store_true", help="Regenerate only the .states file; leave existing .trans/.pf untouched. State indices are unaffected by states-only fixes (e.g. label decoding), so existing transitions stay valid.")
    parser.add_argument("--delay", type=float, default=0.0, help="Seconds to sleep between ions during batch processing.")
    parser.add_argument("--overwrite", action="store_true", help="Regenerate outputs even if they already exist.")
    parser.add_argument("--keep-going", action="store_true", help="Continue batch processing if one ion fails.")
    parser.add_argument("--yes", action="store_true", help="Required for --all processing without --list-only.")
    parser.add_argument("--timeout", type=int, default=60)
    parser.set_defaults(save_raw=True, save_intermediate=True, stream_transitions=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.validate:
        validate_outputs(Path(args.data_root), Path(args.manifest_path), Path(args.report_dir))
        return

    if args.all:
        if args.charge is not None:
            raise SystemExit("--charge cannot be combined with --all")
        if not args.list_only and not args.yes:
            raise SystemExit("Refusing full-site processing without --yes. Run --all --list-only first, then add --yes when ready.")
        discoveries = discover_all_ions(args.timeout, args.base_url, args.offset, args.limit, verbose=True)
    else:
        if not args.element:
            raise SystemExit("--element is required unless --all is used")
        args.element = args.element[:1].upper() + args.element[1:].lower()
        if args.element not in ELEMENTS:
            raise SystemExit(f"Unknown element symbol: {args.element}")
        charges = [args.charge] if args.charge is not None else list(range(args.start_charge, args.max_charge + 1))
        discoveries = [discover_ion(args.element, charge, args.timeout, args.base_url) for charge in charges]

    if args.manifest:
        write_manifest(discoveries, Path(args.manifest))

    if args.list_only:
        for discovery in discoveries:
            print_discovery(discovery)
        report_blocked(args)
        return

    process_discoveries(discoveries, args)
    report_blocked(args)
    report_incomplete(args)


def report_blocked(args: argparse.Namespace) -> None:
    if not BLOCKED:
        return
    total = write_blocked(Path(args.blocked_report))
    print(f"\n{len(BLOCKED)} forbidden (403) file(s) hit this run; {total} total recorded in {args.blocked_report}")


def write_incomplete(path: Path) -> int:
    """Write ions whose transitions fell short of the declared count, merged with
    any rows from a previous run so batched runs keep accumulating into one
    table. De-duplicated by ion code. Returns the total count."""
    headers = [
        "element", "charge", "code", "ion", "gam_files", "declared_transitions",
        "records_read", "transitions_written", "dropped", "completeness",
        "trans_sources",
    ]
    merged: dict[str, dict[str, str]] = {}
    if path.exists():
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                merged[row["code"]] = row
    for entry in INCOMPLETE:
        row = {key: str(value) for key, value in entry.items()}
        merged[row["code"]] = row

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers)
        writer.writeheader()
        for row in sorted(merged.values(), key=lambda item: item["code"]):
            writer.writerow(row)
    return len(merged)


def report_incomplete(args: argparse.Namespace) -> None:
    if not INCOMPLETE:
        return
    total = write_incomplete(Path(args.incomplete_report))
    print(
        f"\n{len(INCOMPLETE)} ion(s) produced incomplete transitions this run; "
        f"{total} total recorded in {args.incomplete_report}"
    )


if __name__ == "__main__":
    main()
