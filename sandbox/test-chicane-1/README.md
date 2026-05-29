# test-chicane-1: analytic four-SBEND C-chicane placement seed

This case is the first OPALX input for the symmetric four-bend C-chicane
outlined in `sandbox/BenTest/bend-tests.tex`.

The design is intentionally hard-edge and analytic:

- `L_b = 1.0 m`
- `theta = 0.1 rad`
- `D = 1.5 m`
- `C = 2.0 m`
- physical bend signs: `(+theta, -theta, -theta, +theta)`

The OPALX input keeps the survey symbolic.  The top-level `REAL` block defines
`LB`, `TH`, `D_CH`, `C_CH`, the OPAL-compatible chord radius
`R_CH = LB / (2 * sin(TH / 2))`, and derived survey points such as `B1_X`,
`B1_Z`, `M1_X`, `M1_Z`, ..., `MEXIT_X`, `MEXIT_Z`.  The element definitions
then use these symbols for `L`, `ANGLE`, `X`, `Y`, `Z`, and `THETA` rather
than frozen decimal coordinates.  This is the executable counterpart of the
chicane formulas in
`sandbox/BenTest/bend-tests.tex`.

The entrance is centered at `(X, Y, Z) = (0, 0, 0)` and the exit is centered at
`X = MEXIT_X = 0`, `Z = MEXIT_Z = 8.980013537482710 m`, parallel to the
entrance.

The later optics target is the small-angle school formula

```text
R56 ~= -theta^2 * (4 L_b / 3 + 2 D) = -4.333333333333334e-2 m
```

Run from this directory:

```sh
/Users/adelmann/git/opalx/build/src/opalx test-chicane-1.in
(cd data && /Users/adelmann/.venv-h6/bin/python test-chicane-1_ElementPositions.py --export-vtk)
(cd data && /Users/adelmann/.venv-h6/bin/python test-chicane-1_ElementPositions.py --export-web)
```

The generated `data/test-chicane-1_ElementPositions.txt` table is the primary
placement check.  The generated VTK and HTML files are useful visual inspection
artifacts for the bend placement.

The input also carries a five-particle momentum scan for an R56 check:

| Particle | delta | pz = beta gamma |
| --- | ---: | ---: |
| 1 | -2.0e-3 | 1954.035026333805 |
| 2 | -1.0e-3 | 1955.992977261995 |
| 3 | 0.0 | 1957.950928190185 |
| 4 | 1.0e-3 | 1959.908879118375 |
| 5 | 2.0e-3 | 1961.866830046566 |

With `OPTION, ASCIIDUMP = FALSE`, the temporal monitors write H5 files.  Use
`test-chicane-1_exit.h5` and fit the final local monitor `z` coordinate against
the initial distribution deltas above.

```text
R56 = -d z_OPALX_exit / d delta
R56_expected ~= -theta^2 * (4 L_b / 3 + 2 D) = -4.333333333333334e-2 m
```

The minus sign is the convention conversion from the OPALX local monitor `z`
coordinate to the usual transport `R56` sign for this chicane.  To extract the
number after a run:

```sh
/Users/adelmann/.venv-h6/bin/python check-r56.py
```

With the current input and `DT = 1.0e-11 s`, the extractor gives
`R56 = -4.334628765343e-2 m`, about `-1.30e-5 m` away from the small-angle
formula.

## Old OPAL R56 deck

`test-chicane-1_opal.in` is the OPAL 2022 comparison deck.  It uses old
OPAL-T `PARALLEL-T` tracking and reference-line `ELEMEDGE` placement:

- `B1, B2, B3, B4` are placed at `s = 1.0, 3.5, 6.5, 9.0 m`.
- The bend signs are `(+theta, -theta, -theta, +theta)`.
- `test-chicane-1_opal_distribution.txt` carries the same five-particle
  momentum scan as the OPALX input, but in old OPAL `FROMFILE` format.
- The `MEXIT` temporal monitor writes `test-chicane-1_opal_exit.h5`.

Old OPAL's analytic SBEND requires the default `1DPROFILE1-DEFAULT` field map.
The deck sets `GAP = 1.0e-6 m` so the finite fringe is negligible for this
hard-edge R56 probe.

Run from this directory:

```sh
source /Users/adelmann/OPAL-2022.1/etc/profile.d/opal.sh
opal test-chicane-1_opal.in
/Users/adelmann/.venv-h6/bin/python check-r56.py test-chicane-1_opal_exit.h5
```

The old OPAL run gives
`R56 = -4.337546378253e-2 m`, about `-4.21e-5 m` away from the small-angle
formula.  The monitor is placed downstream of a `S0 = 1.0 m` offset, but the
fit uses the local monitor `z`, so the absolute `s` offset does not enter the
R56 slope.

The current placement table gives the expected C-chicane survey:

| Point | Z [m] | X [m] |
| --- | ---: | ---: |
| B1 entry | 0.000000000 | 0.000000000 |
| B1 exit / D12 begin | 0.998750260 | 0.049979169 |
| B2 entry | 2.491256508 | 0.199729294 |
| B2 exit / D23 begin | 3.490006769 | 0.249708464 |
| B3 entry | 5.490006769 | 0.249708464 |
| B3 exit / D34 begin | 6.488757029 | 0.199729294 |
| B4 entry | 7.981263277 | 0.049979169 |
| B4 exit / MEXIT | 8.980013537 | 0.000000000 |

Use the `BEGIN`/`END` rows and `ENTRY EDGE`/`EXIT EDGE` rows for this
placement check.  The generated `FIELD BEGIN`/`FIELD END` rows for drifts are
currently not a reliable visual reference for this input.
