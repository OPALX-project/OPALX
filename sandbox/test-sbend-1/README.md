# test-sbend-1: SBEND placement check

This directory contains the executable input decks and helper scripts for the
single-particle analytic `SBEND` placement case.

The authoritative geometry, formulas, OPALX/OPAL comparison plots, and current
numerical results are documented in `../BenTest/bend-tests.tex` and
`../BenTest/bend-tests.pdf`.  Keep this README limited to local run mechanics.

## Files

- `test-sbend-1.in`: OPALX deck using explicit `X`, `Y`, `Z`, `THETA`, `PHI`,
  and `PSI` placement.
- `test-sbend-1_opal.in`: OPAL 2022.1 comparison deck using `ELEMEDGE`.
- `run_opal_opalx_compare.sh`: runs OPALX and old OPAL in separate temporary
  directories and compares `*_DesignPath.dat` plus `.stat` diagnostics.

## Run OPALX only

```sh
/Users/adelmann/git/opalx/build/src/opalx test-sbend-1.in
```

The main placement diagnostic is `test-sbend-1_Monitors.stat`.  The `.loss`
file can be header-only for this one-particle monitor case.

## Compare with old OPAL

```sh
./run_opal_opalx_compare.sh
```

The helper defaults to `/tmp/opalx-sbend-compare`, runs OPALX in `opalx/`, runs
OPAL 2022.1 in `opal-2022/`, and writes CSV/Markdown summaries plus PNG plots
under `comparison/`.  The comparison intentionally ignores `.h5` and `.lbal`.

## Local acceptance

Use the monitor row with the smallest local in-plane residual.  In the current
monitor statistics format, `mean_x`, `mean_y`, and `mean_s` are columns 15, 16,
and 17:

```sh
awk '$1=="test-sbend-1_exit_monitor" {
  r=sqrt($15*$15 + $17*$17)
  if (n++==0 || r<best) { best=r; mx=$15; my=$16; ms=$17; t=$3 }
}
END {
  print "best monitor residual =", best, \
        "mean_x =", mx, "mean_y =", my, "mean_s =", ms, "t =", t
}' test-sbend-1_Monitors.stat
```

The practical checks are: OPALX exits with status 0, the monitor stat file is
written, no `nan` or OPALX error appears, and the best
`sqrt(mean_x^2 + mean_s^2)` residual stays below `5.0e-5 m`.
