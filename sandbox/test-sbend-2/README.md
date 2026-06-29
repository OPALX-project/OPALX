# test-sbend-2: SBEND vertical focusing comparison

This directory contains the vertical-offset `SBEND` comparison against old
OPAL 2022.1.

The bend model, fringe-field formulas, marker locations, OPALX/OPAL plots, and
current numerical comparison are documented in `../BenTest/bend-tests.tex` and
`../BenTest/bend-tests.pdf`.  Keep this README limited to local run mechanics.

## Files

- `test-sbend-2.in`: OPALX deck with explicit placement and a vertically
  offset particle.
- `test-sbend-2_opal.in`: old OPAL comparison deck using `ELEMEDGE`.
- `run_opal_opalx_compare.sh`: runs both codes in isolated work directories
  and compares `*_DesignPath.dat` plus `.stat` diagnostics.

## Run and compare

```sh
./run_opal_opalx_compare.sh /tmp/opalx-sbend-2-compare
```

The script writes `summary.md`, `geometry_markers.csv`, CSV comparison tables,
and PNG plots under `/tmp/opalx-sbend-2-compare/comparison`.  The plot
`vertical_particle_motion.png` compares `.stat` `mean_y` versus `s` from the
line start; because the input is a replicated monoenergetic particle, this is
the vertical particle trajectory.

The comparison intentionally uses only `*_DesignPath.dat` and `.stat`; HDF5
and load-balance output are ignored.
