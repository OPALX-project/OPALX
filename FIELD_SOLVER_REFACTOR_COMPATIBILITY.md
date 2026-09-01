# Field Solver Refactor Compatibility Report

## Integration baseline

This report records the semantic merge of `master` into `refactor-take-two` and the resulting
self-field parity work.

- Refactor parent: `631876866d1e46433a06430ba4d118e9917ec31d`
- Master parent: `a5c8592401c9ff3b2d7f383cfc621ffa734f95d3`
- Common ancestor: `a2adf3b070c16394f5f68f3a7a9b8c9f11bf7310`
- Pinned IPPL source: `65172536f5d2693d50c0c1fb7bc57ac7ac8a346c`
- Host verification: Apple Silicon, Clang 21.1.8, OpenMPI 5.0.9, Release/SERIAL

Both merge parents and the integrated tree were built against that same IPPL checkout. This avoids
attributing dependency changes to the solver integration.

## Result

The merged implementation keeps the refactor's single tracker-facing boundary:

```text
ParallelTracker -> SelfFieldSystem::solve(SolveContext&)
                -> selected SelfFieldAlgorithm
                   |- Pic3DSolver
                   |  |- FFT periodic
                   |  |- OPEN
                   |  `- P3M long range + short range
                   `- Pic2d5Solver
```

The deleted `FieldSolver`, `BinnedFieldSolver`, `FieldContainer`, `LoadBalancer`,
`ImageChargeScatterController`, and `ippl::PicManager` orchestration were not restored. Parser
objects are converted once to immutable configuration. `PartBunch` owns particle storage and a
neutral particle-layout descriptor, while each concrete algorithm owns its field workspace,
backend, frame policy, and algorithm-specific state.

`TrackRun` preserves master's construction and destruction constraints. Configuration and particle
storage are created before sampling or restore, the solver system is created after particle state is
available, and the tracker is created last. Destruction therefore releases the tracker before the
solver and the solver before particle storage.

## Common interface changes

- `SelfFieldAlgorithmKind` now includes `Pic2d5`.
- `SelfFieldAlgorithmConfig` now contains `Pic3DConfig` or `Pic2d5Config`.
- `PoissonBackendKind` now includes `P3M`, with an immutable cutoff in `Pic3DConfig`.
- `SolverCapabilities` supplies `PrimaryOnly` or `AllTrackingActive` particle selection.
- `SelfFieldSystem` uses configuration visitors for requested physics, diagnostics, correction
  reporting, and reported bin count.
- `ParticleSetView` exposes controlled solve-selection mutation without exposing concrete storage.
- `ParticleLayoutConfig` describes spatial or overlap layout ownership, cutoff, and particle
  boundary behavior at construction time.

The public execution signatures remain unchanged:

```cpp
void SelfFieldAlgorithm::execute(SolveContext&, SelfFieldDiagnostics&);
void SelfFieldSystem::solve(SolveContext&);
```

## P3M parity

P3M is a typed Pic3D backend, not a separate legacy solver island. Its long-range path uses
`FFTTruncatedGreenSolver_t` with master's `alpha = 2 / RCUT`, force constant, GPU-aware FFT
parameters, and OPEN or PERIODIC boundary selection. A focused `P3MShortRangeInteraction` runs
`TruncatedGreenParticleInteraction` after grid gather with master's charge and normalization
conventions.

P3M particle storage uses an overlap layout. OPEN layouts use non-wrapping particle boundaries and
PERIODIC layouts use periodic boundaries. Layout update and migration calls are independent of the
concrete layout type.

Configuration rejects:

- non-positive `RCUT`;
- binning;
- mixed or unsupported field boundary conditions;
- source-plane corrections;
- neutralizing-background subtraction.

## FFT2D5 parity

Production FFT2D5 logic now lives under `src/SpaceCharge/Pic2d5/` and is independent of the deleted
solver hierarchy. `Pic2d5Solver` retains no `PartBunch*`, parser object, or borrowed Pic3D field
storage. The module contains separate reference-path, frame-policy, persistent-workspace, and
per-slice Poisson responsibilities.

The first solve lazily loads the configured or default design path after orbit threading. Resource
construction is idempotent and publishes initialized state only after the reference path, device
copy, persistent fields, and all slice solvers are valid.

The implementation preserves master's selected compatibility ordering and restores temporary
position and field transforms on exceptions. It also retains master's atomic charge deposition,
volume-density normalization, bilateral transverse gather, guarded longitudinal gradient,
single-slice handling, pipe modes, closed-ring behavior, optional longitudinal scattering, and
configurable reference-path filename.

FFT2D5 deliberately remains single-rank and rejects binning, corrections, and redistributed field
layouts.

## Other merged behavior

- `DataSink` combines master's checkpoint output rewind with explicit reported-bin-count input.
- `ParallelTracker` retains master's restart scheduling, preparation state, aperture handling,
  monitor suppression, and step reporting while making one common self-field call.
- New master distribution types and warning fixes are retained without distribution ownership of a
  field container.
- The old solver source and tests remain deleted. Their coverage is replaced by configuration,
  selection, Pic3D, P3M layout, short-range, Pic2d5 component, and FFT2D5 production tests.

## Verification evidence

### Builds and unit tests

- Clean Release/SERIAL build with `OPALX_ENABLE_UNIT_TESTS=OFF`: pass.
- Clean Release/SERIAL build with `OPALX_ENABLE_UNIT_TESTS=ON`: pass.
- Complete CTest suite after integration: 87/87 pass.
- Focused checkpoint, multi-container, particle-container, configuration, and Pic2d5 tests: pass.
- Checked-in `Pic2d5EndToEnd` CTest: pass, two tracking steps and 16 per-slice backend solves.
- End-to-end `TYPE=NONE` checkpoint restart: pass. A fresh 100-step run checkpointed at global
  step 90, advanced to step 100, then restored step 90, rewound `.stat` and HDF5 output, and
  reproduced the final ten steps.

Build trees:

- `/Users/aliemen/Desktop/opalx/build_master_merge_baseline`
- `/Users/aliemen/Desktop/opalx/build_refactor_merge_baseline`
- `/Users/aliemen/Desktop/opalx/build_refactor_merge_tests`

### P3M regressions

`Drift-3-p3m-periodic-fromfile` passes all eight checked-in assertions exactly or within floating
point roundoff on one and two ranks. Its artifacts are:

- `build_refactor_merge_baseline/smoke/Drift-3-p3m-periodic-fromfile/merge-parity-periodic-20260901T204736.155987Z`
- `build_refactor_merge_baseline/smoke/Drift-3-p3m-periodic-fromfile/merge-p3m-np2-20260901T210234.577041Z`

`dih-p3m-open` runs all 250 steps successfully. Against the repository reference, both the pinned
master control and integrated build miss the same two tight emittance tolerances. This is baseline
or platform drift, not an integration regression. Comparing the integrated result directly to the
pinned master control passes all eight assertions. The maximum relative integrated-to-master
difference among those checks is `1.15e-6`.

Control and parity artifacts:

- `build_master_merge_baseline/smoke/dih-p3m-open/master-control-open-20260901T205021.025442Z`
- `build_refactor_merge_baseline/smoke/dih-p3m-open/merge-vs-master-open-20260901T205038.536527Z`

### Post-merge MPI parity

The exact staged implementation also passes:

- `Drift-3-open-integrated-fromfile` on two ranks with shifted Green: 8/8 assertions;
- `SwissFEL-booster-SC-emittedfromfile` on two ranks with shifted Green and ORB migration: 9/9
  assertions.

Artifacts:

- `build_refactor_merge_baseline/smoke/Drift-3-open-integrated-fromfile/merge-shifted-green-np2-20260901T210206.112630Z`
- `build_refactor_merge_baseline/smoke/SwissFEL-booster-SC-emittedfromfile/merge-orb-np2-20260901T210217.765476Z`

### Retained pre-merge parity matrix

The refactor parent already passed the established one-rank FFT, OPEN, integrated Green, binning,
emission, and correction matrix, plus selected two-rank OPEN, shifted-Green, and ORB cases. This
merge keeps the same Pic3D implementation as the ownership base and recompiles it with master's
tracker and parser changes. The complete unit suite and new production gates pass after that
integration.

## Dependency audit

The common headers expose no IPPL field, Cartesian mesh, `FieldLayout`, ORB, Pic3D workspace,
binning-plan, or correction-plan type. `ParallelTracker` depends only on `SolveContext` and
`SelfFieldSystem` for self-field execution. Concrete IPPL types remain below `SpaceCharge/Ippl`,
`SpaceCharge/Pic3D`, and `SpaceCharge/Pic2d5`.

No production orchestration reference remains to the deleted legacy classes. A final
`git diff --check` and dependency scan are part of merge finalization.

## External acceptance still required

This Mac has no CUDA or HIP compiler or compatible GPU. The following production evidence must be
collected in the existing Linux accelerator CI environments:

1. clean `PLATFORMS=CUDA` and `PLATFORMS=HIP` production compiles;
2. P3M short-range runtime coverage on both backends;
3. shifted-Green communication and layout migration with GPU-aware MPI;
4. a multi-rank GPU-aware P3M run.

FFT2D5 distributed-memory support remains intentionally out of scope.
