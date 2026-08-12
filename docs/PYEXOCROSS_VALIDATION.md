# PyExoCross validation pipeline

## Scope

PyExoCross is an independent consumer of the Kurucz-derived `.states`,
`.trans`, and `.pf` outputs. It does not replace `process_kurucz_atom.py` and
does not modify source data. The pipeline checks format compatibility, internal
consistency, selected physical relationships, partition functions, radiative
lifetimes, and LTE stick spectra. A successful PyExoCross run does **not** prove
that the line list is complete or universally correct.

The implementation was inspected and tested against Python 3.10.12 and
PyExoCross 1.1.9. The installed public signatures are keyword-driven:

- `pyexocross.partition_functions(inp_filepath=None, **kwargs)`
- `pyexocross.lifetimes(inp_filepath=None, **kwargs)`
- `pyexocross.stick_spectra(inp_filepath=None, **kwargs)`
- `pyexocross.cross_sections(inp_filepath=None, **kwargs)`

The adapter creates a disposable standard ExoAtom `atom/dataset` tree and
`.adef.json`, then calls those public APIs. It never patches PyExoCross.

## Installation

```bash
python3 -m pip install -r requirements-validation.txt
```

Confirm the actual environment:

```bash
python3 -c "import pyexocross; print(pyexocross.__version__, pyexocross.__file__)"
```

## Input schemas

Files are headerless and whitespace/fixed-width compatible.

`.states` columns:

1. state ID (integer)
2. energy (cm-1)
3. total degeneracy `g`
4. `J`
5. uncertainty (cm-1)
6. lifetime (s; `nan`/`inf` may occur)
7. Landé g factor
8. configuration
9. term
10. `NI`/`CA` source flag

The current converter explicitly defines `g = 2J + 1`, so this relationship is
enforced by default. Disable it in configuration only for a different,
documented schema.

`.trans` columns are upper state ID, lower state ID, Einstein A (s-1), and
wavenumber (cm-1). Upper/lower semantics are checked from state energies, never
from ID ordering.

`.pf` columns are temperature (K) and dimensionless Q. Non-uniform grids are
valid and reported. Linear interpolation is used only inside the original
temperature bounds; extrapolation is rejected.

## Single-ion run

This verified Sc III command selects three temperatures within the input PF
range:

```bash
python3 -m validation.run_pyexocross \
  --ion Sc-III \
  --charge 2 \
  --states Kurucz-Nist-Overlap-data/Sc-III/Sc_III__Kurucz.states \
  --trans Kurucz-Nist-Overlap-data/Sc-III/Sc_III__Kurucz.trans \
  --pf Kurucz-Nist-Overlap-data/Sc-III/Sc_III__Kurucz.pf \
  --temperatures 1000 10000 50000 \
  --range 1000 100000 \
  --range-unit cm-1 \
  --output reports/Sc-III \
  --config config/pyexocross_validation.yaml
```

If temperatures are omitted, the CLI chooses low/mid/high points from the
actual `.pf` coverage. Wavelength output uses
`lambda_nm = 1e7 / wavenumber_cm-1` and is explicitly vacuum wavelength.
Input wavelength intervals are also supported with `--range MIN MAX
--range-unit nm` (or `um`); they are converted to an ordered wavenumber
interval without extrapolation or an unstated medium change.

## Batch run

```bash
python3 -m validation.run_pyexocross \
  --input-root Kurucz-Nist-Overlap-data \
  --all-ions \
  --output reports/all \
  --config config/pyexocross_validation.yaml
```

Batch mode discovers only directories containing all three required files.
Large transition files are parsed in configurable chunks.

Useful switches include `--skip-spectrum`, `--skip-lifetime`, `--skip-pf`,
`--strict`, `--cross-section`, repeated `--temperature`, and an explicit
`--charge` consistency check. Invalid arguments or failed validation produce a
non-zero exit code.

## Output

Each ion report contains:

- `validation_summary.json`
- `validation_report.md`
- `states_issues.csv`
- `transitions_issues.csv`
- `pf_comparison.csv`
- `lifetime_comparison.csv`
- `strongest_lines.csv`
- `partition_function.png`
- `partition_function_log.png`
- `partition_function_relative_error.png`
- `stick_spectrum_wavenumber.png`
- `stick_spectrum_wavelength.png`
- `stick_spectrum.csv`
- `run.log`
- `pyexocross_raw/` with official PyExoCross outputs

Overall status is `PASS`, `PASS_WITH_WARNINGS`, or `FAIL`. It is based on input
validation and calculation outcomes, not merely process completion.

## Tolerances

Defaults in `config/pyexocross_validation.yaml` are:

- transition absolute wavenumber tolerance: 0.5 cm-1
- transition relative tolerance: 1e-5
- PF monotonic relative tolerance: 1e-8
- `g = 2J + 1`: enforced because the current converter defines it

A transition mismatch is reported when both its absolute and relative
tolerances are exceeded. PF differences include absolute difference, relative
difference, mean/RMS/maximum relative error, and log ratio.

## Lifetime interpretation

The calculated lifetime is `tau_i = 1 / sum_f(A_if)` for spontaneous channels
whose first transition column is the upper state. Output is seconds. States
without valid decay channels receive a status and NaN rather than division by
zero. Missing transitions reduce `sum(A)` and therefore make calculated
lifetimes too long. Metastable states are not automatically invalid.

## Cross-sections

Cross-section support is opt-in. With no reliable pressure-broadening data,
the default Gaussian calculation is labelled **exploratory /
format-validation result** and must not be cited as a precision physical
prediction. Atomic high-ionization data should not silently inherit molecular
pressure-broadening assumptions. Profile, pressure, bin size, range, and cutoff
are configurable.

## High-ionization cautions

`X III` means charge +2, not +3. The validator enforces:

```text
spectroscopic stage = charge + 1
electron count = Z - charge
```

The Kurucz four-digit directory code is confirmed by the current discovery
implementation as `xxyy`, with `xx = Z` and `yy = charge`, including higher
ions. This is stronger than extrapolating the older 00/01-only description.

Finite state lists can underestimate high-temperature Q. The pipeline does not
invent ionization-limit truncation, continuum lowering, or occupation
probability models. Agreement between a `.pf` and the same `.states` generation
path is internal/circular consistency, not independent physical validation.
External NIST comparisons must be identified separately.

## Known PyExoCross 1.1.9 behavior

The installed reader requires exact PF temperature matches and does not
interpolate. The adapter adds explicitly linearly interpolated, in-range points
to the disposable staged PF only. Multi-temperature stick calls can produce a
short Q array on sparse grids; the adapter therefore makes one documented
public API call per temperature. No third-party source is modified.

## Tests

```bash
python3 -m pytest -q
```

The fixture is small and reference-complete. Tests cover Roman numeral/charge
mapping, all parsers, state references, energy differences, PF ordering,
reference PF, interpolation boundaries, lifetimes, missing channels, invalid A,
missing PyExoCross, and a real public-API end-to-end smoke run.

## Common failures

- “ion must look like Fe-III”: use a canonical Roman numeral.
- charge mismatch: remember that stage III is charge +2.
- missing state reference: regenerate/filter transitions against the matching
  states file.
- PF temperature outside range: choose an in-range temperature; extrapolation
  is deliberately forbidden.
- no stick lines: check wavenumber range and transition/state consistency.
- `FAIL` despite PyExoCross output: inspect the preflight errors and runtime
  errors in `validation_summary.json`.
