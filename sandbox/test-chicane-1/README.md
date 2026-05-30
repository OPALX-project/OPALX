# test-chicane-1: analytic four-SBEND C-chicane

This directory contains the executable OPALX and OPAL 2022.1 decks for the
symmetric four-`SBEND` C-chicane.

The chicane survey formulas, small-angle `R56` estimate, OPALX/OPAL comparison
numbers, and plots are documented in `../BenTest/bend-tests.tex` and
`../BenTest/bend-tests.pdf`.  Keep this README limited to local run mechanics.

## Files

- `test-chicane-1.in`: OPALX deck with symbolic survey variables and a
  five-particle momentum scan for `R56`.
- `test-chicane-1_opal.in`: OPAL 2022.1 comparison deck.
- `test-chicane-1_opal_distribution.txt`: old OPAL `FROMFILE` distribution for
  the same momentum scan.
- `check-r56.py`: fits the exit-monitor H5 data and reports the transport
  `R56` convention used in the bend-test document.

## Run OPALX and export geometry

```sh
/Users/adelmann/git/opalx/build/src/opalx test-chicane-1.in
(cd data && /Users/adelmann/.venv-h6/bin/python test-chicane-1_ElementPositions.py --export-vtk)
(cd data && /Users/adelmann/.venv-h6/bin/python test-chicane-1_ElementPositions.py --export-web)
/Users/adelmann/.venv-h6/bin/python check-r56.py test-chicane-1_exit.h5
```

The generated `data/test-chicane-1_ElementPositions.txt` table is the primary
placement check.  Use the `BEGIN`/`END` rows and `ENTRY EDGE`/`EXIT EDGE` rows.
Generated `FIELD BEGIN`/`FIELD END` rows for drifts are not currently a
reliable visual reference for this input.

## Run old OPAL comparison

```sh
source /Users/adelmann/OPAL-2022.1/etc/profile.d/opal.sh
opal test-chicane-1_opal.in
/Users/adelmann/.venv-h6/bin/python check-r56.py test-chicane-1_opal_exit.h5
```

The old OPAL deck writes `test-chicane-1_opal_exit.h5`.  The OPALX and old OPAL
`R56` values should be compared through `check-r56.py`, not by absolute monitor
`s` positions, because the old OPAL deck carries a reference-line offset.
