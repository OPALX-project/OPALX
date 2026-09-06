# Cyclotron sector and trim-coil one-particle plan

Date: 2026-09-06

## Source documentation review (2026-09-06)

Reviewed the new cyclotron sector, map reader, parser/trim-coil definition, and
directed turn counter against their implementations. Expanded Doxygen with SI
units, local/global coordinate conventions, field formulas and signs, half-open
sector boundaries, interpolation preconditions, legacy mirrored-coil limitations,
device-view lifetime expectations, and endpoint-based turn-counting semantics.
Map views are explicitly documented as read-only by convention, not by their types.
No physics or numerical behavior changed. Targeted Doxygen generation passed with
documentation warnings treated as errors; `git diff --check` passed. The review
does not constitute a documentation audit of the entire source tree.
Next step: retain these documentation improvements with the cyclotron implementation.

## Goal and scope

Read the PSI Ring median-plane magnetic field map and its trim-coil definition in
OPALX, place eight explicit `CYCLOTRONSECTOR` occurrences in a `RING`, and track
one coasting proton for one physical revolution.  Compare the result step by step
with OPAL 2022.1.  Tune extraction, space charge, and RF cavities are out of scope
for this first milestone.

## Repository and reference state

- Local OPALX `master` was fast-forwarded to `c8bfffb08`, the head of PR #550.
  PR #550 contains PR #548 in its ancestry, so both changes are present locally.
- Development branch: `545-add-cyclotron-coasting-beam`, fast-forwarded to the
  merged master. The user authorized pushing this branch, not master.
- Reference inputs are in `sandbox/cyclotron/opal/`.
- Old OPAL is `/Users/adelmann/OPAL-2022.1/bin/opal` after sourcing
  `/Users/adelmann/OPAL-2022.1/etc/profile.d/opal.sh`.
- OPAL 2022.1 source used for inspection is tag `OPAL-2022.1.0`, commit
  `0d1bf29f85c1799203116a6781db7fbb178b918c`.

## Baseline already established

A temporary, non-destructive variant of `cyclotron1.in` was run in
`/tmp/opalx-cyclotron-one-particle` with:

- kinetic energy: 72 MeV;
- initial radius: 2.1314 m;
- initial radial mechanical momentum: -0.000240 beta-gamma;
- initial tangential mechanical momentum computed by old OPAL: 0.39920181 beta-gamma;
- one particle from `dist2.dat`;
- 720 steps and `STEPSPERTURN=720`;
- magnetic PSI field map and `PSI-BFIELD-MIRRORED` trim coil;
- no RF and no space charge.

Old OPAL completed one turn in 118.46 ns and travelled 13.167 m.  Its final
reported state was approximately `r=2.13112 m`, `theta=-0.6182 deg` after step
720.  `cyclotron1-trackOrbit.dat` contains the initial state and steps 1--719;
the step-720 state is printed only in the log.  The comparison harness must account
for this output convention.

The supplied `ic.dat` ends near `r=4.3002 m`. Source inspection during implementation
corrected the original assumption about trim-coil support: `4.350--4.470 m` are
shape radii, not hard cutoffs. The mirrored model has radial tails (negligible at
72 MeV). Separate field tests near the shape radii exercise the substantial effect.

## Implementation state (2026-09-06)

- Added the PSI map reader, sector scalar/device evaluation, named mirrored trim
  definition, parser/visitor wiring, and directed reference turn counter.
- Added `TestCyclotronSector` and `opalx/coasting72.in` with the legacy launch state.
- Full build passed in the existing `omp-build` with `-j 10`.
- The initial compile exposed the Kokkos 5 mirror-type rename and IPPL expression
  conversion; both were corrected. Numeric (not whitespace-token) skipping is
  required for the PSI format's packed numeric header block.
- One 72 MeV proton now completes its first directed turn in 722 steps at
  DT=0.164527805199081 ns, t=118.789075354 ns, path=13.202970643 m.
- The tracked particle stays within 3.5e-13 m of the reference at this timestep;
  energy drift is below 3.2e-11 MeV.
- Old OPAL's default is RK4 (`TrackCmd.cpp`); explicit LF2 selects its Boris
  drift-kick-drift (`Steppers/LF2.hpp`). Both integrators were verified in the
  installed executable's runtime logs on one rank. The explicit RK4 reproducer
  produces a byte-identical orbit dump to the original baseline.
- Matching LF2 against OPALX at the same timestep gives max position difference
  8.53e-9 m and momentum difference 1.41e-9 beta-gamma, including old ASCII precision
  and rounded launch momentum. Against RK4 the position differences are 335, 76,
  and 36 micrometres at DT, DT/2, and DT/4.
- `TestCyclotronSector` (7 cases), `TestIndexMap`, `TestOrbitThreader`,
  `TestOpalBeamlinePlacement`, and `TestRing` passed.
- A separate two-rank comparison completed before the user's correction and was
  identical to one rank. The user subsequently specified that these runs use one
  rank; the reproducible comparison script now runs exclusively on one rank.
- Manuals were found in `/Users/adelmann/git/opalx-manual` (the paths previously
  listed in AGENTS.md do not exist). User and physics pages were added there;
  manual validation passed. Those changes are separate from the OPALX checkout.
- Directed-turn completion uses the first timestep endpoint past the return plane;
  its longitudinal overshoot is bounded by one timestep and must converge with DT.
- Explicit TURNS restart is rejected until counters/plane state are persisted.
- Reproducers: `reference72/coasting72-rk4.in`, `reference72/coasting72-lf2.in`,
  `opalx/coasting72.in`, and `compare_coasting72.py`. Run the old inputs after
  sourcing `/Users/adelmann/OPAL-2022.1/etc/profile.d/opal.sh` from `reference72`.
- Field map SHA256: `0bd65560cfe7c92df55e64a17058d6d856cc3b75d4737a2a100bcbbe97b0d305`.
- Final single-rank script passed in `results72-single-rank`, including the LF2
  agreement assertion. The relevant unit tests passed after the final code changes.
  Both new manual pages render successfully. Validation used the OpenMP CPU build;
  GPU compilation/execution has not been checked.
- Reviewed the source diff. Whitespace checks pass for code and authored inputs;
  the raw PSI map retains two original header lines with trailing spaces so its
  bytes and SHA256 remain unchanged.
- First milestone complete; next work is user review and selection of the next
  energy/initial condition. Tunes and RF remain future work.

## Integration checks already run

- `cmake --build omp-build -j 4`: passed after incorporating PRs #548 and #550.
  Use `-j 10` for subsequent builds.
- `TestIndexMap`: passed.
- `TestOrbitThreader`: passed.
- `TestOpalBeamlinePlacement`: passed.
- `TestRing`: passed.

## Agreed user model for the first milestone

The user accepted this design on 2026-09-06, including directed return-plane
crossings for turn counting. The first coasting milestone is implemented.

Keep `TRIMCOIL` as a named field-model definition, not a third beamline element.
Each explicitly placed sector references the definition:

```opal
TC1: TRIMCOIL, TYPE=PSI-BFIELD-MIRRORED,
     RMIN=4.350, RMAX=4.470, BMAX=1.4E-3, SLPTC=600.0;

SM0: CYCLOTRONSECTOR, FMAPFN="bfield.dat", SYM=8,
     RMIN=1.9, RMAX=4.7, VMIN=-0.05, VMAX=0.05,
     BSCALE=1.0, TRIMCOIL={TC1},
     X=..., Y=..., Z=..., THETA=..., PHI=0, PSI=0;
```

OPALX public units should be SI: metres, tesla, radians, seconds, and mechanical
momentum in beta-gamma internally.  The example `SLPTC` value above is expressed
in 1/m; old OPAL accepted `0.6 1/mm`.  The final attribute spelling for a list
must follow the existing OPALX string-array grammar.

For the initial full-azimuth mirrored coil, the same immutable trim-coil model is
referenced by all eight sectors.  Partial-azimuth trim coils need an explicit
decision later: their `PHIMIN/PHIMAX` must be defined in the enclosing RING frame,
not repeated independently in every sector-local frame.

## Implementation sequence

### 1. Freeze scalar field reference data

Create a small old-OPAL reference table at selected radii, sector angles, and
vertical offsets.  Include points on both sides of interpolation-cell boundaries,
the 0/45-degree seam, and map radial boundaries.  Generate paired trim-coil-on and
trim-coil-off values near 4.35, 4.41, and 4.47 m.

Acceptance: the reference data records coordinates, units, base-map field,
trim-coil contribution, and total field without relying only on trajectories.

### 2. Implement the PSI sector-map reader

Add a dedicated `CyclotronSectorFieldMap` rather than forcing the PSI NAR-like
format through the existing Cartesian `Fieldmap` factory.  On the host it should:

- parse `rmin`, `dr`, `theta_min`, `dtheta`, `nrad`, and `ntheta`;
- validate `RMIN/RMAX` against the file instead of cropping;
- convert millimetres to metres and kilogauss to tesla at the input boundary;
- reproduce the periodic seam explicitly;
- compute the old five-point radial and angular derivatives;
- retain immutable host data and Kokkos device views;
- cache raw map data by canonical filename so eight sectors do not allocate eight grids.

Keep `BSCALE` on the element, outside the shared immutable map.

Acceptance: parser/header, seam, derivative, bilinear interpolation, malformed-file,
and shared-allocation unit tests pass.  Selected median-plane values match old OPAL.

### 3. Port the trim-coil definition and first model

Add `OpalTrimCoil` plus a small immutable trim-coil model interface.  Port only
`PSI-BFIELD-MIRRORED` initially.  Preserve the old analytical profile while using
SI units internally.  Attach named trim-coil models to `OpalCyclotronSector` during
input update; device evaluation must use trivially copyable parameters, not host
polymorphic pointers.

The old model subtracts its analytic profile from Bz and its unmirrored derivative
times z from Br. Preserve that legacy derivative convention on both halves of the
profile; do not silently replace it with a differentiated mirrored function. The
proper rotation into OPALX is `(x,y,z)_old -> (x,-z,y)_new`, including field vectors.

Acceptance: analytic values, radial symmetry about the coil midpoint, azimuth
gate, zero-strength case, and on/off differences agree with old OPAL.

### 4. Implement `CYCLOTRONSECTOR`

Add:

- `OpalCyclotronSector` for parsing and validation;
- `CyclotronSector` for geometry, containment, transforms, and field application;
- a new `ElementType::CYCLOTRONSECTOR` and visitor plumbing.

Use the issue #545 convention: OPALX ring plane X-Z, vertical Y, with each Mode-A
pose describing the minimum-azimuth entrance boundary.  In sector-local
coordinates, transform the stored ring centre once and compute `(r, theta, v)`.
Containment is an annular wedge with half-open angular intervals; the last radial
interpolation cell must never read beyond the grid.

Convert the reconstructed cylindrical field to local Cartesian components and let
`OpalBeamline` perform the existing local/lab rotations.  Implement both the
reference-particle scalar path and a Kokkos particle-container path.

Acceptance: exact boundary tests, eight-sector seam tests, rotation covariance,
host/device agreement, and field superposition with the trim coil.

### 5. Define a cyclotron-compatible RING stopping rule

Do not use the field map's midpoint-radius circumference as the physical orbit
length.  The merged `TRACK,TURNS` implementation currently sets its stopping path
from `Ring::getLength()`, but cyclotron equilibrium-orbit circumference varies with
energy and need not equal `2*pi*(RMIN+RMAX)/2`.

Count physical directed return-plane crossings, with hysteresis to avoid counting
the initial position or numerical chatter as another turn. This is the agreed
stopping rule; measured path length is a diagnostic, not a substitute for crossing
counts. It must not require an arbitrary reference radius on every sector.

Acceptance: `TURNS=1` ends at the first directed return crossing for both inner and
outer cyclotron orbits, independent of the sector-map midpoint radius.

### 6. Map the 72 MeV initial condition into OPALX

At old-OPAL `PHIINIT=0`, map the state into the OPALX X-Z ring plane as:

```text
R = (r0, 0, 0)
P = (pr0, 0, +sqrt((beta*gamma)^2-pr0^2))
```

The RING source frame must be placed at that global position with its longitudinal
axis along the initial tangential momentum.  The one-particle distribution then has
zero offsets; `RINIT`, `PRINIT`, and `PHIINIT` do not become sector attributes.

Acceptance: OPALX prints the same initial global Cartesian state and momentum as
the transformed old-OPAL row before any integration step.

### 7. One-particle regression

Run old OPAL and OPALX with the same initial state and matched physical time step.
Compare in global Cartesian coordinates at every common time sample:

- position and mechanical momentum components;
- radius and azimuth;
- sampled magnetic-field components;
- accumulated path length;
- energy/momentum-norm conservation;
- final directed return-plane state.

First compare the base map with trim-coil strength zero, then enable the coil.  Run
a one-rank time-step convergence study. Following the user's 2026-09-06
instruction, this run family uses one rank exclusively.

No tolerance should be chosen until the scalar field comparison separates parser,
interpolation, coordinate/sign, and integrator differences.

## Expected source areas

```text
src/Elements/OpalCyclotronSector.{h,cpp}
src/AbsBeamline/CyclotronSector.{h,cpp}
src/Fields/CyclotronSectorFieldMap.{h,cpp}
src/TrimCoils/TrimCoilModel.*
src/TrimCoils/TrimCoilMirrored.*
src/TrimCoils/OpalTrimCoil.*
src/AbsBeamline/ElementBase.*
src/AbsBeamline/BeamlineVisitor.*
src/OpalConfigure/Configure.cpp
src/ValueDefinitions/StringConstant.cpp
unit_tests/Fields/
unit_tests/TrimCoils/
unit_tests/AbsBeamline/
```

## Deferred work

- tune extraction and the two-particle SEO convention;
- cyclotron `RFCAVITY` thin-gap kicks;
- other old trim-coil fit types;
- partial-azimuth trim-coil placement semantics;
- space charge and neighbouring bunches;
- continuous 3D RF field maps.
