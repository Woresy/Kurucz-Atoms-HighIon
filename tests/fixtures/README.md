# Minimal ExoAtom fixture

This fixture follows the exact ten-column `.states`, four-column `.trans`, and
two-column `.pf` schemas emitted by `process_kurucz_atom.py`. It is a deliberately
small, physically self-consistent C III subset-style dataset: every transition
references an included state, each wavenumber equals the energy difference, and
the lifetimes are `1 / sum(A)`. It is not presented as an observed line list and
is used only for deterministic software tests.
