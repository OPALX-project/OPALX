# test-rbend-2: single analytic RBEND vertical edge-focusing comparison

This directory mirrors `sandbox/test-sbend-2` for `RBEND`.  A particle is
launched with a `1 mm` vertical offset and zero vertical angle, then OPALX is
compared against old OPAL 2022.1 using `DesignPath.dat` and `.stat` diagnostics.

## Analytic configuration

| quantity | value |
| --- | ---: |
| upstream drift | `0.4 m` |
| `L` | `1.0 m` rectangular body length |
| `ANGLE` | `45 deg` |
| `E1` | `0 deg` |
| `E2 = ANGLE - E1` | `45 deg` |
| `HGAP` | `0.02 m` |
| `FINT` | `0.5` |
| initial vertical offset | `1 mm` |

For `RBEND`, the exit edge angle is fixed by `ANGLE - E1`.  This test therefore
keeps `E1 = 0` and uses the `45 deg` exit edge as the main vertical-focusing
diagnostic.

The OPALX deck derives the explicit rectangular body pose instead of inserting
raw numbers:

```text
ENTRY_LOCAL_X = -0.5 * L * tan(ANGLE / 4)
ENTRY_LOCAL_Z =  0.5 * L
BODY_THETA = -ANGLE / 2
```

The hard-edge exit point is

```text
EXIT_X = -L * tan(ANGLE / 2)
EXIT_Z = DRIFT_LENGTH + L
EXIT_THETA = -ANGLE
```

## Run

```sh
/Users/adelmann/git/opalx/build/src/opalx test-rbend-2.in
```

The OPALX/old OPAL comparison helper is:

```sh
./run_opal_opalx_compare.sh
```

It defaults to `/tmp/opalx-rbend-2-compare`, runs OPALX and old OPAL in
separate directories, and writes CSV/Markdown summaries plus PNG diagnostics in
`comparison/`.  The plots include the vertical particle motion from the `.stat`
`mean_y` column.

