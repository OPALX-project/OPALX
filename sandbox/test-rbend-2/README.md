# test-rbend-2: RBEND vertical focusing comparison

This directory contains the vertical-offset `RBEND` comparison against old OPAL
2022.1.

The rectangular-bend geometry, edge-focusing setup, OPALX/OPAL plots, vertical
motion figure, and current numerical comparison are documented in
`../BenTest/bend-tests.tex` and `../BenTest/bend-tests.pdf`.  Keep this README
limited to local run mechanics.

## Files

- `test-rbend-2.in`: OPALX deck with explicit rectangular-bend placement and a
  vertically offset particle.
- `test-rbend-2_opal.in`: old OPAL comparison deck using `ELEMEDGE`.
- `run_opal_opalx_compare.sh`: runs both codes in isolated work directories
  and compares `*_DesignPath.dat` plus `.stat` diagnostics.

## Run OPALX only

```sh
/Users/adelmann/git/opalx/build/src/opalx test-rbend-2.in
```

## Compare with old OPAL

```sh
./run_opal_opalx_compare.sh
```

The helper defaults to `/tmp/opalx-rbend-2-compare`, runs OPALX and old OPAL in
separate directories, and writes CSV/Markdown summaries plus PNG diagnostics in
`comparison/`.  The plots include the vertical particle motion from the `.stat`
`mean_y` column.  The comparison intentionally ignores `.h5` and `.lbal`.
