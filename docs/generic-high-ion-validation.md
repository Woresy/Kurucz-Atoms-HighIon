# Generic higher-ionization validation

The reusable entry point is:

```bash
python3 -m validation.run --ion Fe-III
```

It discovers Kurucz and NIST files recursively from the configured data roots,
validates each available component, performs quantum-aware cross-source
matching, runs both datasets through PyExoCross with identical parameters, and
writes an ion-specific report.

## Common commands

```bash
# One ion
python3 -m validation.run \
  --ion Fe-III \
  --kurucz-root Kurucz-Nist-Overlap-data \
  --nist-root Nist-temp-data \
  --output-root validation/reports \
  --temperatures 1000 3000 6000 \
  --wavelength-range-nm 275 330

# Several ions, without editing validation source code
python3 -m validation.run \
  --ions Fe-III Fe-IV Ti-IV \
  --output-root validation/reports

# Every complete Kurucz ion
python3 -m validation.run --all-ions --output-root validation/reports

# Derive the complete wavelength interval from every positive transition
python3 -m validation.run \
  --ion Ti-IV \
  --all-wavelengths \
  --output-root validation/reports-full-spectrum

# Fast format/matching smoke test
python3 -m validation.run --ion Ti-IV --skip-pyexocross
```

Requested temperatures must be covered by both partition-function files when
both sources are to be run with identical parameters.

`--all-wavelengths` and `--wavelength-range-nm MIN MAX` are mutually exclusive.
The all-wavelength mode uses the union of the positive Kurucz and NIST
transition-wavenumber ranges. Very small fine-structure wavenumbers can produce
extremely large far-infrared wavelengths, so full-spectrum CSVs and plots may
be large or visually compressed.

## Matching rules

State IDs are never compared across databases. Matching uses normalized
configuration, normalized term, J, available parity, and energy proximity.
Results are classified as:

- `EXACT_QUANTUM_MATCH`
- `HIGH_CONFIDENCE`
- `AMBIGUOUS`
- `UNMATCHED`

Transition comparison is then restricted to the resulting matched-state map.
It reports wavelength, wavenumber, Einstein A and
`log10(A_Kurucz) - log10(A_NIST)`.

## Missing data

Each source is classified as `AVAILABLE`, `PARTIAL`, `NOT_AVAILABLE`, or
`INVALID`. An unavailable NIST component skips only the dependent comparison;
Kurucz validation and PyExoCross processing continue.

An empty metric is represented by `NaN` or a documented availability status.
The framework does not invent reference values.

## Interpretation

- Kurucz-only data may be predicted data or a coverage difference.
- NIST-only data may be observed reference coverage absent from Kurucz.
- Neither case alone proves a conversion error.
- PyExoCross success proves executable compatibility, not scientific
  correctness.
- The NIST state export has no lifetime column, so lifetime comparison is
  explicitly reported as `NOT_AVAILABLE` rather than fabricated.

## Current generic assumptions

- Configuration punctuation and principal-shell digits are normalized for
  matching while the original labels remain in output.
- Missing parity weakens a match instead of preventing it.
- PyExoCross 1.1.9 requires ten whitespace-delimited state columns. A disposable
  NIST compatibility view supplies unavailable lifetime/Landé fields and
  replaces internal configuration whitespace with dots; original data are not
  modified.
