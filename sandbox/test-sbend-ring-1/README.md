# test-sbend-ring-1: ten-SBEND ring seed

This directory is a first OPALX-only sandbox input for a closed, zero-order
ring made from ten analytic `SBEND` elements and ten `0.25 m` drifts.  It is
based on the single-particle setup in `../test-sbend-2`: explicit placement,
one vertically offset electron, no space charge, and the same OPAL-style fringe
parameters.

## Zero-order geometry

The first deck represents one turn as an ordered line:

```text
MSTART, B1, D1, B2, D2, ..., B10, D10, MCLOSE
```

The last drift closes the polygon back to the starting plane.  This is not yet
a periodic tracking setup; it is a one-turn placement and tracking seed.

Adopted design values:

| quantity | value |
| --- | ---: |
| number of bends `N` | `10` |
| chord length per `SBEND`, `L_c` | `1.0 m` |
| drift length between bends, `L_d` | `0.25 m` |
| bend angle per cell, `theta = 2 pi / N` | `0.628318530717959 rad = 36 deg` |
| OPAL-compatible sector radius, `rho = L_c / (2 sin(theta/2))` | `1.618033988749895 m` |
| reference arc length per bend, `L_ref = rho theta` | `1.016640738463052 m` |
| zero-order circumference, `C0 = N (L_ref + L_d)` | `12.66640738463052 m` |

OPALX `SBEND` placement follows the chord-to-sector convention used in the
single-bend tests.  With a positive OPALX bend angle, the reference tangent
rotates by `-theta` in the global `X-Z` drawing plane.  For a bend entry point
`p_i = (x_i, z_i)` and entry tangent angle `a_i`, define

```text
R(a) (x, z) = (cos(a) x + sin(a) z,
              -sin(a) x + cos(a) z)

u(a) = (sin(a), cos(a))

b_half = (rho (cos(theta/2) - 1), rho sin(theta/2))
b_full = (rho (cos(theta)   - 1), rho sin(theta))
```

The zero-order recurrence used in the input file is

```text
body_i  = p_i + R(a_i) b_half
q_i     = p_i + R(a_i) b_full
a_{i+1} = a_i - theta
drift_i = q_i
p_{i+1} = q_i + L_d u(a_{i+1})
```

with `p_0 = (0, 0)` and `a_0 = 0`.  After ten cells,
`a_10 = -2 pi` and `p_10 = p_0` to roundoff in the zero-order geometry.
For explicit OPALX drifts in this sandbox, `X`, `Y`, `Z` are set to the drift
start, not the midpoint.

## Beam

The input launches one electron with kinetic energy `250 MeV`.  Using
`m_e c^2 = 0.510998950 MeV`,

```text
gamma      = 1 + 250 / 0.510998950 = 490.2377958897959
beta       = 0.9999979195519578
beta gamma = sqrt(gamma^2 - 1)     = 490.23677597553325
```

The FROMFILE distribution starts at the ring entrance with
`x = 0`, `y = 1 mm`, `z = 0`, `px = py = 0`, and
`pz = beta gamma`.

## Run

From this directory:

```sh
/Users/adelmann/git/opalx/build/src/opalx test-sbend-ring-1.in
```

Optional geometry exports after a successful run:

```sh
(cd data && /Users/adelmann/.venv-h6/bin/python test-sbend-ring-1_ElementPositions.py --export-vtk)
(cd data && /Users/adelmann/.venv-h6/bin/python test-sbend-ring-1_ElementPositions.py --export-web)
```

Because the closed ring is not monotonic in global `Z`, `ZSTOP` is deliberately
set far away and `MAXSTEPS` is the practical stop condition for this first
one-turn seed.
