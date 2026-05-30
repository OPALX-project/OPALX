# test-rbend-1: single analytic RBEND placement consistency test

This directory mirrors `sandbox/test-sbend-1` for the OPALX-native analytic
`RBEND`.  The OPALX input uses explicit `X`, `Y`, `Z`, `THETA`, `PHI`, and
`PSI` placement; the old OPAL input uses legacy `ELEMEDGE` placement.

## Analytic configuration

| quantity | value |
| --- | ---: |
| upstream drift | `0.3 m` |
| `L` | `1.0 m` rectangular body length |
| `ANGLE` | `45 deg` |
| `E1` | `0 deg` |
| `E2 = ANGLE - E1` | `45 deg` |
| `HGAP` | `0.02 m` |
| `FINT` | `0.5` |

For an `RBEND`, `L` is the straight rectangular body length.  OPAL uses the
pole-face chord radius and design reference path length

```text
R_CH = L / (sin(E1) + sin(E2))
L_ref = abs(ANGLE) * abs(R_CH)
```

and the hard-edge exit point for an entrance edge at `(0, 0, DRIFT_LENGTH)` is

```text
EXIT_X = -L * tan(ANGLE / 2)
EXIT_Z = DRIFT_LENGTH + L
EXIT_THETA = -ANGLE
```

The OPALX body pose is the rectangular body-frame origin.  The input derives it
from the local entry-edge vector

```text
ENTRY_LOCAL_X = -0.5 * L * tan(ANGLE / 4)
ENTRY_LOCAL_Z =  0.5 * L
BODY_THETA = -ANGLE / 2
```

by rotating that vector into the lab and subtracting it from the requested
entry-edge location.

## Run

```sh
/Users/adelmann/git/opalx/build/src/opalx test-rbend-1.in
```

The OPALX/old OPAL comparison runs each code in a separate temporary directory:

```sh
./run_opal_opalx_compare.sh
```

The helper defaults to `/tmp/opalx-rbend-compare`, writes CSV/Markdown summaries
under `comparison/`, and generates PNG diagnostics for the design path,
momentum, and sampled `B_y`.  The comparison intentionally ignores `.h5` and
`.lbal` and uses only `*_DesignPath.dat` and `.stat`.
