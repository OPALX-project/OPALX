# test-rbend-1: RBEND placement check

This directory contains the executable input decks and helper script for the
single-particle analytic `RBEND` placement case.

The rectangular-bend geometry, OPAL-compatible reference-path convention,
OPALX/OPAL comparison plots, and current numerical results are documented in
`../BenTest/bend-tests.tex` and `../BenTest/bend-tests.pdf`.  Keep this README
limited to local run mechanics.

## Files

- `test-rbend-1.in`: OPALX deck using explicit `X`, `Y`, `Z`, `THETA`, `PHI`,
  and `PSI` placement.
- `test-rbend-1_opal.in`: OPAL 2022.1 comparison deck using `ELEMEDGE`.
- `run_opal_opalx_compare.sh`: runs OPALX and old OPAL in separate temporary
  directories and compares `*_DesignPath.dat` plus `.stat` diagnostics.

## Run OPALX only

```sh
/Users/adelmann/git/opalx/build/src/opalx test-rbend-1.in
```

## Compare with old OPAL

```sh
./run_opal_opalx_compare.sh
```

The helper defaults to `/tmp/opalx-rbend-compare`, writes CSV/Markdown
summaries under `comparison/`, and generates PNG diagnostics for the design
path, momentum, and sampled `B_y`.  The comparison intentionally ignores `.h5`
and `.lbal`.
