# FFT2D5 Integration Audit

## Status

FFT2D5 is implemented as an independent `SelfFieldAlgorithm` under
`src/SpaceCharge/Pic2d5/`. This document records the resulting ownership and behavior. It replaces
the earlier future-extension audit.

The tracker contains no FFT2D5 branch. Its only self-field action is
`SelfFieldSystem::solve(SolveContext&)`, the same boundary used by Pic3D.

## Ownership

| Responsibility | Owner |
|---|---|
| immutable parser snapshot | `Pic2d5Config` |
| all-active particle selection | `SelfFieldSystem` using `SolverCapabilities` |
| design-path loading and device copy | `ReferencePath` |
| master-compatible outer frame ordering | `Pic2d5FramePolicy` |
| Frenet mapping and boost/deboost operations | `Pic2d5Solver` device helpers |
| persistent 3D staging fields and line density | `SliceWorkspace` |
| persistent per-slice 2D fields and OPEN solvers | `SliceWorkspace::Slice` |
| scatter, slice solves, gradient, and gather | `Pic2d5Solver` |

No component stores a `PartBunch*`, `FieldSolverCmd*`, or runtime parser reference. FFT2D5 does not
replace `PartBunch` field storage and does not inherit from `FieldSolver` or `BinnedFieldSolver`.

## Configuration and capabilities

`Pic2d5Config` snapshots:

- three-dimensional staging resolution;
- longitudinal field mode: OPEN, CIRCULAR, PLATES, or NONE;
- transverse pipe sizes and beam radius;
- closed-ring policy;
- longitudinal scattering selection;
- configured or default design-path filename.

Its capabilities select every tracking-active particle container. The common system applies that
selection before validating required readable and writable particle attributes. FFT2D5 reports no
binning or correction request and a compatibility bin count of one.

Parser and immutable configuration validation reject multiple MPI ranks, binning, source-plane
corrections, and parallel field-layout dimensions.

## Initialization contract

Construction does not require the design-path file because orbit threading has not produced it yet.
On the first execution, `ensureInitialized()`:

1. loads and validates at least two distinct reference-path points;
2. computes the path length;
3. transfers the path to persistent device storage;
4. creates the staging mesh, layouts, fields, density views, and per-slice backends;
5. publishes initialized state only after every resource is valid.

A failed attempt leaves the solver uninitialized and retryable. Later solves reuse every allocation
and backend.

## Execution sequence

For one common solve, `Pic2d5Solver`:

1. enters the solver-owned compatibility frame policy;
2. transforms selected particles to Frenet coordinates and boosts momenta;
3. deposits charge atomically into the persistent staging density;
4. applies cell-volume and time normalization;
5. handles open or closed longitudinal boundaries;
6. copies each slice into its persistent 2D backend and solves transverse Poisson;
7. constructs line density and its guarded longitudinal gradient;
8. performs the bilateral transverse gather and optional longitudinal field calculation;
9. deboosts fields and maps them back to lab coordinates;
10. restores temporary particle and field transforms on normal and exceptional exits.

The bilateral gather and compatibility frame ordering intentionally match current master, including
its selected historical conventions. A cleaner frame contract remains a possible later change, but
is not part of this parity merge.

## Physics behavior retained from master

- atomic charge deposition;
- charge-to-volume-density normalization;
- optional two-dimensional or three-dimensional longitudinal scatter;
- bilateral transverse gather from adjacent longitudinal slices;
- guarded gradient indexing and single-slice handling;
- OPEN, CIRCULAR, PLATES, and NONE longitudinal modes;
- closed-ring density, field, and gradient wrapping;
- configurable reference-path filename;
- field and position restoration when execution throws.

## GPU and lifetime audit

Reference-path and workspace data remain device resident after initialization. Device kernels
capture views and scalar policy values rather than parser or owning runtime objects. CUDA-visible
device helpers are public or standalone as required by the repository's CUDA compilation rule.

Workspace fields are declared before the solver array so reverse destruction releases backends
before the fields and layouts they borrow. Particle views are obtained for each execution instead
of being retained across storage changes.

The module uses no distributed communication. Its one-rank validation is therefore a deliberate
correctness boundary, not an accidental decomposition choice.

## Verification

`TestPic2d5Solver` covers lazy initialization, all supported pipe modes, persistent slice creation,
and per-slice solve diagnostics through the common API. Configuration tests cover algorithm kind,
neutral layout selection, and the all-tracking-active capability.

`Pic2d5EndToEnd` is a checked-in CTest fixture that runs the production executable for two tracking
steps. Orbit threading generates `data/Pic2d5EndToEnd_DesignPath.dat`; the first common solve loads
it lazily; eight persistent slice backends execute on each step. The test passes only when the run
completes and reports all 16 backend solves.

The Release/SERIAL unit-test build and complete 87-test suite pass after integration.

## Deferred work

- distributed-memory FFT2D5;
- batched or pooled slice FFT execution;
- a device-resident line-density reduction that avoids sequential slice orchestration;
- a redesigned reference-frame contract that intentionally changes master's ordering;
- CUDA and HIP production compile and hardware runtime evidence.
