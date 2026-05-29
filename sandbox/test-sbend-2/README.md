# test-sbend-2: vertical pole-face focusing comparison

This case extends `test-sbend-1` from pure placement to a first optics check.
It deliberately launches one particle with a vertical offset and uses the entry
pole-face focusing of the analytic `SBEND` to compare OPALX against old OPAL
2022.1.

## Analytic design

The OPALX input uses explicit placement. The old OPAL input uses `ELEMEDGE`.
Both decks use the same sector-bend chart as `test-sbend-1`:

```text
R       = L / (2 sin(ANGLE/2))
L_ref   = R ANGLE
r_body  = r_entry + (R(cos(ANGLE/2)-1), 0, R sin(ANGLE/2))
r_exit  = r_entry + (R(cos(ANGLE)-1),   0, R sin(ANGLE))
t_exit  = (-sin(ANGLE), 0, cos(ANGLE))
```

The `0.4 m` upstream drift is long enough that the negative-`E1` entrance
fringe support is present inside the tracked region. OPALX places the SBEND
body from derived `REAL` variables instead of hard-coded coordinates.

| quantity | value |
| --- | ---: |
| upstream drift | `0.4 m` |
| `L` | `1.0 m` |
| `ANGLE` | `45 deg` |
| `E1` | `-49.99 deg` |
| `E2` | `0 deg` |
| `HGAP` | `0.02 m` |
| `FINT` | `0.5` |
| OPAL full `GAP` | `0.04 m` |
| Enge half width `a = 5 GAP` | `0.2 m` |
| entry support `a / |cos(E1)|` | `0.311079876857156 m` |
| exit support `a / cos(E2)` | `0.2 m` |

This separates three longitudinal locations that can otherwise look like a
placement offset in plots:

| marker | `s` from line start |
| --- | ---: |
| upstream Enge support start | `0.088919934618024 m` |
| nominal entrance edge | `0.4 m` |
| nominal exit edge | `1.426172152977031 m` |

The field is mathematically nonzero at the upstream support start, but it is
tiny there. With the comparison script's `1e-6` peak-field support threshold,
the visible design-path `B_y` fringe starts at about `0.240 m` in both OPALX
and old OPAL. The nominal bend edge remains at `0.4 m` in both decks.

The vertical edge force is applied with an Enge-derivative weight over the full
fringe ramp, matching the longitudinal centroid of old OPAL's map. In the
current comparison OPALX ends at `mean_y = 2.168092e-3 m` and old OPAL at
`2.161588e-3 m`.

The `1 mm` vertical offset is large enough to read cleanly in diagnostics while
remaining paraxial. The present `E1 = -49.99 deg` is kept from the historical
vertical-focusing case so that changes in the current OPALX fringe/placement
implementation are visible against old OPAL.

## Run and compare

Run:

```sh
./run_opal_opalx_compare.sh /tmp/opalx-sbend-2-compare
```

The script runs OPALX and old OPAL in separate directories and writes
`summary.md`, `geometry_markers.csv`, CSV comparison tables, and PNG plots under
`/tmp/opalx-sbend-2-compare/comparison`. The plot
`vertical_particle_motion.png` compares `.stat` `mean_y` versus `s` from the
line start; because both distributions are monoenergetic replicated particles,
this is the vertical particle trajectory. The comparison intentionally uses
only `*_DesignPath.dat` and `.stat`; HDF5 and load-balance output are ignored.
