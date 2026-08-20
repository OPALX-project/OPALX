# OPALX 12-Particle BeamBeam Experiment

This directory contains the OPALX version of the CAIN 12-particle comparison.

## Exact-timed comparison

The current comparison lives in `timed/`.  It replaces the older workaround
of generating one deck per pair with one physical clock and one OPALX deck.
The CAIN data describe six artificial electron--positron pairs moving in the
positive longitudinal direction.  They are born at

```text
ct = -0.9000, -0.5994, -0.3006, 0.0000, 0.3006, 0.5994 mm
```

and are sampled every `1.8 um/c = 6.004153713566736 fs` through
`ct = 1.8 mm`.  The tracker begins one interval before the first birth so that
the mirrored primary source is active before a witness is present.  Every
birth and OPALX output lies on the CAIN grid and the final global step is 1500,
so no time interpolation is used.  At an exact floating-point step boundary a
newborn may first appear in H5 at reference step zero or one; the comparison
uses the exact emitted-file row whenever H5 starts at step one.

The emitted files use the explicit named format

```text
x y z px py pz birth_time
```

where positions are offsets from the emission-source `R0`, momenta are
normalized by `m_e c`, and `birth_time` is an offset in seconds from the
source `T0`.  Electron and positron files intentionally contain the same
numeric initial phase-space states; their separate beam containers determine
the particle species and charge sign.

The physical primary and its reflection about the BeamBeam midpoint supply the
co-propagating and oncoming source fields.  The 12 witnesses gather those
fields but are excluded from deposition by `WITNESS_CONTAINERS = "1,2"`.
`BBRIGID = TRUE` holds the sampled source distribution fixed for this
reduced-order comparison.  `APERTURE = "RECTANGLE(2.4e-3, 2.4e-4)"` fixes a
field domain large enough for the CAIN electron excursions.
`DELETEONTRANSVERSEEXIT = FALSE` prevents particle deletion, but it does not
make the particle layout open: an over-kicked particle that reaches the edge
wraps periodically, which the analyzer reports explicitly.

The production parameters are 400000 deterministic primary macroparticles
(NumPy PCG64 seed 20260629) and a `1024 x 128 x 128` mesh.  The transverse cell
sizes are approximately `2.34 um` in x and `1.88 um` in y; substantially
coarser meshes are parser/workflow smoke tests, not trajectory validation.

Prepare and run locally with the OpenMP build:

```sh
cd sandbox/track12particles/opalx/timed
./run_timed_track12.sh
```

On Merlin6, first stage a prepared `production/` case and then submit the
four-A100 smoke and production dependency chain with:

```sh
./merlin/submit_track12_a100.sh submit \
  --run-root /path/to/run \
  --opalx /path/to/cuda-build/src/opalx
```

The smoke uses the full production mesh for three steps.  The production job
starts only after that allocation/parser/MPI check succeeds.  Multi-rank OPALX
is launched with `mpiexec`, with one rank per A100 and
`--kokkos-map-device-id-by=mpi_rank`; plain multi-task `srun` would create MPI
singletons on Merlin.

Fetch a completed production run and make all plots locally:

```sh
~/.venv-h6/bin/python timed/fetch_and_compare_a100.py \
  --remote-dir /path/to/run/production
```

The comparison reconstructs absolute particle coordinates by adding all three
components of the H5 `RefPartR` attribute to the stored `x`, `y`, and `z`
offsets.  Raw H5 IDs are not stable under MPI redistribution, so pair identity
is reconstructed from birth order and one-step phase-space continuity.  The
fixed particle layout is periodic transversely; wrapped H5 positions are
unwrapped for continuous plots while the raw positions are retained in the
pointwise CSV.  It matches all 13012 CAIN samples exactly by species, pair, and
grid step.  For births whose initial state precedes the first H5 record, the
exact generated emitted-file row is included as the step-zero sample.

## Four-A100 production result

The exact case completed on four A100-SXM4-40GB GPUs as Merlin job `353971` in
`1:43:15`.  All twelve witnesses were emitted in the intended order, no
witness contributed to the field solve, and all 13012 CAIN rows were matched.
The input/output mechanics are therefore validated, but this discretization
does **not** reproduce the CAIN trajectories.

| species | samples | RMSE x [um] | relative L2 x | RMSE y [um] | RMSE s [um] |
|---|---:|---:|---:|---:|---:|
| electron | 6506 | 511.68 | 1.12695 | 251.14 | 30.52 |
| positron | 6506 | 5.04 | 0.94299 | 1.08 | 0.78 |

The first transverse kick identifies the failure before long-time integration.
For the central pair, CAIN gives `px = 2.2261e-2`, whereas OPALX gives
`2.6816e-4`, only `1.20%` of the CAIN kick.  Pairs 3 and 6 begin with the wrong
transverse sign.  The equal-and-opposite electron/positron OPALX kicks confirm
that species/container assignment is correct; the source field sampled at the
witness is not.

This is now retained as the **pre-reference-offset-fix baseline**.  The weak
and sign-changing kicks were subsequently traced to witness-coordinate
conversion, not to a non-convergent field solve.  Four electron trajectories
also cross the periodic y boundary (nine wrap events in total); from the first
wrap onward their field sampling is not a physical open-boundary trajectory.

Generated production reports and plots are in the ignored local directory
`timed/a100_400k_1024x128x128/results/`.  In particular:

- `track12_cain_vs_opalx.png`: slide-scale `x(s)` comparison;
- `track12_cain_vs_opalx_wide.png`: full-range `x(s)` comparison;
- `track12_cain_vs_opalx_y.png`: unwrapped `y(s)` comparison;
- `track12_first_kick_summary.csv`: all twelve first-kick values and ratios;
- `track12_comparison.json`: aggregate, per-species, H5-identity, and wrap diagnostics.

## Witness reference offset and corrected first kick

Each particle container stores particle `R` relative to its own reference
particle.  In the timed track12 case the witness reference has
`RefPartR.x = sigma_x`, while its particle-local `R.x` is nearly zero.  The
BeamBeam witness gather previously translated only the longitudinal path
offset.  It therefore sampled an offset witness container near the primary
axis even though H5 correctly reported its absolute position at `x=sigma_x`.

The gather now translates the witness coordinates into the source container's
reference frame using the transverse difference

```text
RefPartR_witness - RefPartR_source
```

and retains `s_witness - s_source` as the authoritative longitudinal offset.
This changes only where a passive witness samples the already-computed source
field.  It does not change source deposition, the Poisson solve, field units,
charge normalization, or MPI reduction order.

`timed/scan_first_kick_fields.py` provides the persistent regression and mesh
diagnostic.  It compares the pair-4 field at `x=sigma_x` with the untruncated
rigid-Gaussian manufactured solution.  Both the normal seeded Gaussian sample
and a deterministic equal-probability tensor source converge monotonically:

| full domain and mesh | source | x cell [um] | OPALX/analytic `|E|` |
|---|---|---:|---:|
| 40 um, 32 x 32 | seeded random | 1.290 | 0.90184 |
| 40 um, 64 x 64 | seeded random | 0.635 | 0.95871 |
| 40 um, 128 x 128 | seeded random | 0.315 | 0.98426 |
| 40 um, 32 x 32 | tensor | 1.290 | 0.89321 |
| 40 um, 64 x 64 | tensor | 0.635 | 0.95323 |
| 40 um, 128 x 128 | tensor | 0.315 | 0.97868 |

The scan names these cases `square20` for their 20 um half-width.  Doubling the
full domain from 40 to 80 um at fixed cell size changes the field by less than 0.1%,
and the random/tensor difference is approximately 1% or less.  The field
direction cosine exceeds 0.99994 at the two finest random meshes.  The CTest
`TestBeamBeamWitnessReferenceOffset` runs the 64 x 64 seeded case and would
fail against the pre-fix on-axis gather.

The exact fixed-source production probe ran on four A100s as Merlin job
`353974`.  It used the same 400000-particle source, 1024 x 128 x 128 mesh, and
2.4 mm x 0.24 mm transverse domain as the trajectory baseline.  It completed
in 23.14 s and gave:

| quantity | value |
|---|---:|
| OPALX `|E|` | 6.569905e9 V/m |
| OPALX/analytic `|E|` | 0.678097 |
| field direction cosine | 0.999995 |
| field-predicted OPALX `Delta px` | 0.023143 |
| CAIN first `Delta px` | 0.022261 |
| OPALX/CAIN first kick | 1.039622 |

The analytic denominator is the untruncated continuum Gaussian, whereas the
OPALX result uses the fixed finite source and the comparatively coarse CAIN
production mesh.  The meaningful trajectory decision is that the corrected
OPALX field predicts the CAIN central first kick within 3.96%, instead of the
pre-fix trajectory's 1.20% residual.  A complete corrected 12-particle run is
therefore justified; it has not yet been performed.

The four-A100 probe can be submitted from a prepared portable run directory
with:

```sh
./timed/merlin/submit_first_kick_field_probe_a100.sh submit \
  --run-root /path/to/prepared/run \
  --opalx /path/to/cuda-build/src/opalx
```

## Files

- `generate_track12_opalx_inputs.py` converts `sandbox/TestParticleOrbit.dat`
  into OPALX `FROMFILE` witness distributions.
- `track12_electrons.fromfile` and `track12_positrons.fromfile` are generated
  inputs with columns `x y z px py pz`.
- `track12_witness_metadata.csv` records the particle names and source-file
  timestamp, since OPALX `FROMFILE` does not carry particle labels.
- `track12_beambeam.in` is the intended comparison input with 400000 primary
  source macroparticles.
- `track12_beambeam_smoke.in` is a cheaper parser/tracking smoke run with 4000
  primary source macroparticles and a coarser time step.
- `track12_one_bunch_electrons.in` is the first validation run: one read-in
  electron bunch only, no BeamBeam interaction.
- `track12_primary_one_beam_smoke.in` tracks one 4000-macroparticle primary
  Gaussian electron beam in a drift, with no BeamBeam interaction.
- `track12_primary_one_beam.in` is the same source-only run at 400000
  macroparticles.
- `track12_primary_with_pairs_smoke.in` tracks the primary beam plus the read-in
  six-electron and six-positron CAIN test particles in a drift, still without a
  BeamBeam element.
- `track12_beambeam_one_source_quick.in` is a fast one-source BeamBeam
  diagnostic. It enables `WITNESS_CONTAINERS` but leaves `COPY_TIME = 0`, so no
  copied counter-propagating source is included.
- `track12_beambeam_copied_source_1fs_probe.in` is a staged-timestep diagnostic
  with a 1 fs fine segment through the copied-source BeamBeam window.
- `generate_timing_pair_inputs.py` creates the pair-wise timing-correct OPALX
  decks under `timing_pairs/`. These decks keep the numeric TestParticleOrbit
  witness coordinates but shift the primary source centroid per pair so the
  first field evaluation uses the insertion timing implied by
  `TestParticleOrbitSimulation.pptx`. They are retained as historical
  diagnostics; use `timed/` for the exact one-clock comparison.

## Coordinate Convention

The CAIN columns are mapped as

```text
x_OPALX = x_CAIN
y_OPALX = y_CAIN
z_OPALX = s_CAIN
px_OPALX = Px_CAIN / (m_e c)
py_OPALX = Py_CAIN / (m_e c)
pz_OPALX = Ps_CAIN / (m_e c)
```

The witness emission sources apply `R0Z = bb_ip_s`, so the CAIN `s` coordinate
is interpreted as a longitudinal offset around the BeamBeam interaction point.

The PowerPoint source slides define the witness `ct` values as insertion times.
For an oncoming primary bunch centered at the IP at `ct = 0`, the source center
at insertion is `s_c = -ct_i`. Because the TestParticleOrbit initial rows use
`s_i = ct_i`, each pair needs a different primary source offset. The generated
pair decks use

```text
primary_source_r0z = bb_ip_s - ct_i
witness_r0z = bb_ip_s
```

so `s_i - s_c = 2 ct_i`; pair 1 starts at the `-3 sigma_z` source edge and
pair 4 starts at the bunch center. This fixes timing only; it intentionally does
not tune source direction, species mapping, field normalization, or pusher
details.

## Commands

Regenerate witness files:

```sh
~/.venv-h6/bin/python generate_track12_opalx_inputs.py
```

Generate pair-wise slide-timed BeamBeam inputs:

```sh
~/.venv-h6/bin/python generate_timing_pair_inputs.py
```

Run the smoke input from this directory:

```sh
/Users/adelmann/git/opalx-beambeam/build_openmp/src/opalx track12_one_bunch_electrons.in
```

Run the source-only primary-beam smoke input:

```sh
/Users/adelmann/git/opalx-beambeam/build_openmp/src/opalx track12_primary_one_beam_smoke.in
```

Run the primary-plus-e-/e+ no-BeamBeam smoke input:

```sh
/Users/adelmann/git/opalx-beambeam/build_openmp/src/opalx track12_primary_with_pairs_smoke.in
```

Run the BeamBeam smoke input from this directory:

```sh
/Users/adelmann/git/opalx-beambeam/build_openmp/src/opalx track12_beambeam_smoke.in
```

Run the 400000-particle comparison:

```sh
/Users/adelmann/git/opalx-beambeam/build_openmp/src/opalx track12_beambeam.in
```
