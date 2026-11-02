# Scaling FFA placement and full-ring tracking

This milestone validates the FFA25 school cell's field/placement conventions and
adds aligned scaling FFA magnets to OPALX's ordinary spatial ring tracking and
explicit `RUN, TURNS=...` path. The subsequent COF extension and the 16-cell school
closed orbit are documented in [SCHOOL-COF.md](SCHOOL-COF.md). The complete notebook
port and built-in FFA spectral tunes remain separate work. Production tracking
tunes and independent cell optics are documented in [TUNES.md](TUNES.md).

## Placement and fields

`SCALINGFFAMAGNET` now uses the same geometrical entrance convention for an
explicit `X/Y/Z/THETA/PHI/PSI` pose and for `ELEMEDGE`. The field centre is half
the nominal magnet length after that entrance. Explicit poses that previously
compensated for the incorrect entrance-centred field must remove that compensation.

For radius `R0`, local nominal coordinates at path length `s` are
`x=R0*(cos(s/R0)-1), y=0, z=R0*sin(s/R0)`. Positive `R0` bends towards negative
x; negative `R0` bends towards positive x. Nominal arc length is positive for
both signs. `ELEMEDGE` FFA placement is absolute on this circle, relative to the
beamline frame, and advances the following element from the resolved exit.

The shared host/device support test uses radial limits, vertical height and
`phi - tan(delta)*log(r/abs(R0)) - phi_centre`. It retains overlapping spiral
fringe tails rather than substituting a straight z slab. Production fields use
the existing batched FFA/Enge/Tanh kernels and accumulate across all selected
occurrences. The simplified test-only device transport oracle explicitly rejects
FFA lattices; it is not the production tracker.

`SchoolScalingFFACellEntrancePosesAndSignedRadii` in
`unit_tests/Elements/TestOpalBeamlinePlacement.cpp` checks the actual school F/D
parameters (R0=4.42 m, field index 7.717022064259446, fourth-order expansion,
Enge coefficients -1.3, 3.42, -1.20, 0.23). It covers:

- Both placement modes, both radius signs, 45° and 55° spiral angles, and a
  translated/rotated beamline frame.
- Median-plane values against an independent polynomial/logistic Enge expression
  and the radial scaling law, including the fringes and off-reference radii.
- First-order off-plane fields against finite-difference median-plane gradients
  required by curl B=0; switching off the D magnet.
- A separate exit-frame test checks a downstream drift after a nonzero circular
  placement and a gap.

The placement test failed before the fixes and passes afterward. The descriptor
suite additionally checks support near radial, vertical and spiral boundaries
and accumulation of overlapping F/D fields off the median plane. The existing
18 field tests retain their previous Maxwell and expansion-convergence checks.

## Full-ring production regression

Run from the OPALX repository with the OPALX Python environment:

```sh
cmake --build build-omp --target opalx_exe TestOpalBeamlinePlacement TestDeviceExternalField TestScalingFFAMagnet TestOpalScalingFFAMagnet -j 6
~/.venv-h6/bin/python docs/validation/scaling-ffa/validate.py --output /tmp/ffa-ring-validation
```

Use a new output directory. `--check-only --output <existing-directory>` rechecks
saved runs without launching OPALX. Each run saves its input, distribution,
probe configuration, log and plane samples; the directory root receives
`summary.json`. `results.json` is the recorded result of this milestone.

The generated ring contains 16 adjacent FFA sectors with overlapping tanh tails,
field index zero and zero spiral angle. Their median-plane fields telescope to
a uniform 0.1 T field. The momentum is chosen independently from
`p[GeV/c] = c*1e-9*B[T]*R[m]` at R=4.42 m. Eight particles have distinct initial
horizontal displacements of -35 to +35 µm. They are tracked for three directed
turns through the real production reference and bunch paths, with no space charge.
The circular orbit is the independent oracle; this is not a test of the
simplified device transport helper.

The matrix contains 12 runs: both radius signs, explicit-pose and ELEMEDGE
placement, 200 ps and 100 ps timesteps on one rank, and 100 ps on two ranks.
Passive samples use the existing internal `OPALX_TEST_PASSIVE_RING_PROBES`
interface, which becomes eligible for aligned FFA rings through spatial support.
This internal interface remains experimental, not a new public PROBE syntax.

Every run must report completion of the third directed turn. Comparisons use the
first two complete returns of all eight particles: at the final reference stop,
a particle's third crossing can fall on either side of the last endpoint bracket.
IPPL assigns interleaved IDs by creation rank whereas FROMFILE distributes
contiguous input blocks, so serial/MPI matching uses the physical input row and
crossing count rather than assuming the same numeric ID across rank counts.

Recorded results:

| Quantity | Result |
|---|---:|
| Completed production runs | 12 / 12 |
| Maximum return-position error, 200 ps | 1.976 µm |
| Maximum return-position error, 100 ps | 0.495 µm |
| Coarse/fine error ratio | 3.998 |
| Largest serial/MPI position difference | 3.90e-12 m |
| Largest serial/MPI momentum difference | 1.18e-13 in p/(mc) |
| Largest serial/MPI crossing-time difference | 8e-20 s |
| Largest pose/ELEMEDGE position difference | 3.04e-12 m |

Assertions are new acceptance criteria, not changes to existing tolerances.
The 3 µm return bound covers the observed second-order Boris/interpolation error
at 200 ps. Refining the timestep must reduce it by at least 2.5× (the expected
asymptotic factor is 4). The 50 µm radial bound includes the imposed 35 µm offsets.
The 1e-7 normalized-momentum bound allows the O(dt²) interpolation error between
accepted endpoints; it is not a relaxation of Boris energy conservation. MPI
bounds are 1e-10 m in position, 1e-10 in normalized momentum and 1e-16 s in time.
Median-plane field checks use 2e-12 T for coordinate/libm rounding; off-plane
finite differences use 2e-7 T/m with 1 µm steps. Existing test tolerances are unchanged.

For this machine's Open MPI configuration, launch unit executables with:

```sh
OMP_NUM_THREADS=1 OMP_PROC_BIND=false mpiexec --host localhost:4 --map-by slot --bind-to none -n 2 build-omp/unit_tests/Elements/TestOpalBeamlinePlacement
```

Seven unit executables (placement, device fields, FFA fields/interface, ring,
passive probes and orbit threader) passed all 94 tests on one rank and again
on two ranks. Counts and exit codes are recorded in `unit-results.json`.

Default CTest launch mapping fails here before assertions; direct launches with
these flags work. The validation script supplies these flags and `--info 1` so
OPALX diagnostics are visible.

## Limits and next work

These checks establish single-cell physics consistency and aligned full-ring
production transport. They do not establish numerical equivalence to old OPAL
for the 12 MeV or fourfold school closed orbits and tunes. A controlled old-OPAL
baseline is still needed. The independent full-ring oracle has zero field index
and zero spiral angle; the nonzero school field index/spiral/off-plane behavior
is covered by the single-cell tests. Subsequent 50-turn school tracking and
independent linear-optics tune checks are documented in [TUNES.md](TUNES.md).

The existing ring restrictions on restart, emission and misalignment eligibility
remain. No claim of new collective-field validation is made. FFA COF support has since been added and validated in [SCHOOL-COF.md](SCHOOL-COF.md).
Built-in FFA spectral tunes remain separate work.

When converting the notebook, also rename the legacy asymmetric Enge attributes:
`X0_START/END` → `START/END_HALF_LENGTH`, `LAMBDA_START/END` →
`START/END_FRINGE_LENGTH`, and `COEFFICIENTS_START/END` →
`START/END_COEFFICIENTS`. Use the current FROMFILE header and source conventions;
FROMFILE momentum must not also be specified through BEAM PC/ENERGY/GAMMA.
