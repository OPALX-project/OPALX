# Future 2.5D Field Solver Extension Audit

## Scope and result

This audit compares the refactored runtime boundary with the prototype in
`25dsolver/src/PartBunch/Solve2d5.h` and `Solve2d5.hpp`. It does not add a 2.5D
implementation or change its tests.

The runtime boundary passes the extension audit. A future `Pic2d5Solver` can be
implemented directly as a `SelfFieldAlgorithm`, consume the existing
collection-shaped `SolveContext`, and keep the tracker free of a concrete 2.5D
branch. The extension requires localized configuration, factory, and setup
work, but it does not require changing `SelfFieldAlgorithm::execute`,
`SolveContext`, or the `ParallelTracker` call site.

## Prototype responsibilities and future ownership

| Current 2.5D prototype responsibility | Future owner |
|---|---|
| reference-path loading and interpolation | `Pic2d5Solver` reference-path component |
| lab/Frenet coordinate conversion | 2.5D frame policy |
| slice mesh, charge and field arrays, line density | 2.5D workspace |
| array of 2D open Poisson solvers | 2.5D backend component |
| scatter, slice solve, gather, lab-field restore | `Pic2d5Solver::execute` orchestration |

The prototype currently combines these responsibilities in `Solve2d5`:

- `orbitThreadersReady()` loads the generated design path, creates the 3D
  staging fields, and allocates the slice mesh, fields, and solver array.
- `doRunSolver()` sequences scatter, slice solves, line-density calculation,
  and gather.
- the device helpers convert particles between lab and Frenet coordinates,
  boost and deboost fields, and calculate the longitudinal field.
- the class inherits `BinnedFieldSolver` and stores a `PartBunch*` plus live
  field pointers.

Only the algorithm and mathematics are reference material. The inheritance,
`PartBunch*` storage, low-level IPPL scatter/gather calls, global parser/output
lookup, and current no-partition slice layout must not be copied into the new
design.

## Existing boundary used by `Pic2d5Solver`

`ParallelTracker` presents particle positions in lab/reference coordinates.
The current 3D implementation performs its beam-frame transform inside
`Pic3DSolver::execute` through `FieldFramePolicy`, then restores positions and
fields before returning. A 2.5D implementation can ignore the supplied 3D
beam transform handles and apply its own reference-path/Frenet mapping.

`ParticleSetView` already supplies:

- all particle containers, an explicit primary index, solve-selection and
  tracking-active flags;
- stable handles for position, momentum, charge, mass, time step, invalid
  mask, and writable electric and magnetic fields;
- a storage generation used to prevent retaining stale Kokkos views.

The tracker currently marks only the primary container as selected, although
it exposes attribute handles and tracking-active state for every container.
The future extension should add a generic selection policy to
`SolverCapabilities` (`PrimaryOnly` for 3D PIC and `AllTrackingActive` for
2.5D). `SelfFieldSystem::solve` can apply that policy to the mutable particle
view before validation. This keeps the selection algorithm-neutral and avoids
a solver-name branch in the tracker.

The future solver must reacquire native device views from the selected handles
after migration, compaction, or reallocation. Quantities such as per-container
mean longitudinal momentum can be reduced from the exposed momentum attributes
and communicator instead of retaining a `PartBunch` or calling its statistics
methods.

`SolveContext::StepState` supplies the communicator, time step, reference
state, and current lab/reference frame handles. The context remains borrowed
for one call and must never be retained by the solver.

## Future integration steps

1. Add a validated immutable `Pic2d5Config` containing the transverse pipe
   dimensions, 3D staging resolution, longitudinal-field mode, beam radius,
   closed-ring policy, reference-path location, and 2D backend options.
2. Add `Pic2d5` to `SelfFieldAlgorithmKind` and `Pic2d5Config` to
   `SelfFieldAlgorithmConfig`. Extend `SelfFieldConfigBuilder` to snapshot the
   parser values; the runtime solver must retain no parser pointer.
3. Add one `SelfFieldFactory` branch that constructs `Pic2d5Solver` directly
   from immutable configuration and stable particle-container setup resources.
   It must not call `takePicWorkspace()` and the resulting runtime object must
   not store `PartBunch*`.
4. Replace the current `Pic3DConfig` reads inside `SelfFieldSystem` with a
   configuration visitor or algorithm-specific request policy. For 2.5D the
   existing tracker-facing methods return no correction, no binning request,
   and a reported bin count of one. `ParallelTracker` remains unchanged.
5. Resolve setup-time particle-layout lifetime explicitly. `PartBunch`
   currently constructs the IPPL particle layouts from a `PicWorkspace`. For
   2.5D it may retain that object only as the layout lifetime anchor, or the
   future project may extract a neutral particle-layout owner. It must not
   treat the 3D PIC fields as the 2.5D physics workspace.
6. Add a generic particle-selection policy and a selected-state setter/helper
   to `ParticleSetView`; apply the policy in `SelfFieldSystem` before common
   attribute validation.
7. Place the implementation under `src/SpaceCharge/Pic2d5/`. Keep its
   reference path, frame policy, workspace, and backend types private to that
   module.

These changes are confined to setup and algorithm selection. They do not add
Cartesian fields, ORB, binning, or 3D correction records to the common solver
interface.

## Delayed initialization

The prototype cannot allocate its domain until the orbit threader has produced
the reference-path file. The future solver can handle this without restoring
`orbitThreadersReady()` or a concrete solver check in the tracker:

1. Resolve the configured or default output path while building the immutable
   configuration, but do not require the file to exist during construction.
2. At the start of the first `Pic2d5Solver::execute`, call an idempotent
   `ensureInitialized()` host method.
3. Load and validate the reference path, transfer it once through an explicit
   host mirror/deep copy, then construct the persistent 2.5D workspace and
   Poisson backend array.
4. Mark initialization complete only after every resource is valid so a failed
   attempt does not leave a partially initialized solver.

The first self-field call already occurs after orbit threading. Lazy,
solver-owned initialization therefore needs no tracker API change. If a later
solver needs an earlier preparation point, it should use a generic
`SelfFieldSystem` lifecycle hook rather than a solver-kind branch.

## Execute and diagnostics model

`Pic2d5Solver::execute` owns the complete sequence:

1. validate required particle attributes and initialize persistent state;
2. transform selected lab particles to Frenet coordinates without exposing
   Frenet types through `SolveContext`;
3. deposit charge into persistent slice staging fields;
4. run the solver-owned 2D Poisson array and calculate line-density gradients;
5. gather, add the longitudinal field, deboost, and write lab-frame particle
   electric and magnetic fields;
6. restore any temporary particle mutations on every exit path.

The frame contract needs one explicit check during that implementation. The
particle positions delivered to the solver and the points in
`_DesignPath.dat` must be proven to use the same reference coordinate system,
or the required mapping must be supplied as a checked host service. The
existing `trackerToSolve` and `solveToTracker` handles describe the 3D
reference/beam transform and must not be silently reinterpreted as a Frenet or
global-path transform.

The call remains inside the common solve, solve-unit, backend-solve, and field-
composition lifecycle where applicable. Reference-path loading, Frenet
mapping, slice copies, line-density work, and per-slice solves may publish
solver-specific `IpplTimings` labels without expanding the permanent common
diagnostics API.

## GPU and distributed-memory requirements

- Keep the reference path and reusable slice state device-resident after the
  one-time host load. Use explicit mirrors and `deep_copy` for initialization.
- Capture only views and trivial scalar/policy values in Kokkos lambdas. Do not
  capture `this`, shared pointers, diagnostics objects with host state, or
  parser/configuration objects.
- Put CUDA-visible device helpers in public or standalone policy functions as
  required by this repository's CUDA build constraint.
- Reacquire particle views after any storage generation change and preserve
  the GPU-aware MPI path; do not introduce unconditional host staging.
- Design the slice layout and 2D backend ownership for multiple ranks. The
  prototype's no-partition layout and one solver object per slice are not a
  distributed-memory contract.
- Allocate fields, line-density arrays, and solver work buffers once. Do not
  allocate them per step or per slice solve.
- Precompute segment lengths, cumulative path length, and stable transported
  frames. The prototype scans the full path and recomputes arc length for every
  particle, which scales as particle count times path-segment count and assumes
  a fixed planar normal.
- Replace the prototype's one host-synchronizing line-density reduction per
  slice with a device-resident reduction strategy, and isolate the sequential
  per-slice FFT loop so a pooled or batched backend can replace it.
- Check or wrap the longitudinal-gradient index before reading it on device;
  the prototype can index outside the line-density-gradient view.
- Declare mesh/layout and field storage before backends that borrow them so
  reverse destruction releases the Poisson solver array first.

## Dependency audit

The common headers `SelfFieldAlgorithm.h`, `SolveContext.h`,
`ParticleSetView.h`, `SolverCapabilities.h`, and `SelfFieldSystem.h` include no
IPPL field, Cartesian mesh, `FieldLayout`, ORB, Pic3D workspace, binning-plan,
or correction-plan type. `ParallelTracker` includes only the common
`SolveContext` and `SelfFieldSystem` boundary for self-field work. Runtime code
under `SpaceCharge/Pic3D` and `SpaceCharge/Ippl` stores no `PartBunch*`,
`FieldSolverCmd*`, `BinningCmd*`, or `EmissionSource*`.

The audit therefore leaves one independent implementation path:
`Pic2d5Solver` owns its geometry, workspace, and backend, while the common call
and tracker remain unchanged.
