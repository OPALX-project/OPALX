# Field Solver Refactor - Implementation Plan (Take Two)

## Starting point

- Worktree: `/Users/aliemen/Desktop/opalx/refactor-take-two`
- Branch: `refactor-take-two`
- Base commit: `a2adf3b070c16394f5f68f3a7a9b8c9f11bf7310`
- The worktree was created from the committed `HEAD` of `opalx`. Uncommitted changes in the original `opalx` worktree are intentionally not present.
- The implementation guide is `solver-discussion-presentation/main.pdf`, backed by `content.tex` and the current/proposed class diagrams.
- The primary local style reference is `src/PartBunch/Binning/`.
- The future compatibility reference is the separate `25dsolver/` worktree. No 2.5D implementation is part of this refactor.

This plan changes no solver code yet. Unit tests are explicitly out of scope for this pass. The production executable and existing regression tests are the verification path.

## Refactor goal

Replace the current control path

```text
ParallelTracker
  -> PartBunch mesh/frame/load-balancing work
  -> BinnedFieldSolver
       -> FieldSolver/IPPL backend
       -> bin iteration
       -> corrections
       -> field composition and diagnostics
```

with this ownership model:

```text
TrackRun owns PartBunch and SelfFieldSystem
  |
  +-> ParallelTracker builds one borrowed SolveContext per solve
        |
        +-> SelfFieldSystem validates capabilities and records diagnostics
              |
              +-> SelfFieldAlgorithm
                    |
                    +-> Pic3DSolver
                    |     +-> PicDomainManager
                    |     +-> PicWorkspace
                    |     +-> IterationPlan -> SolveUnit(s)
                    |     +-> CorrectionPlan -> SolvePass(es)
                    |     +-> PicScatterGather
                    |     +-> IpplPoissonAdapter
                    |     +-> FieldFramePolicy / FieldComposer
                    |
                    +-> future Pic2d5Solver, using the same common boundary
```

The common boundary must not expose Cartesian mesh, `FieldLayout`, ORB, IPPL field, or 3D correction types. Those belong to `Pic3DSolver`.

## Non-negotiable behavior

The refactor must preserve these current semantics until a separate change explicitly alters them:

1. The tracker clears fields and respects the current `MINBINEMITTED` guard.
2. 3D PIC uses the same beam-frame transform, mesh update order, scatter normalization, Poisson solve, Lorentz transform, accumulation, and gather order.
3. Binned solves still rebin/merge/sort, calculate the per-bin mean momentum and $\gamma$, stretch the longitudinal mesh spacing, solve each bin, accumulate fields, and gather once.
4. Image-charge and shifted-Green behavior remain explicit additional solve passes. Their signs, mirror order, solve-call accounting, and field-dump behavior must not change accidentally.
5. A field-layout change refreshes all dependent state in this order: fields/workspace, particle layout and migration, then backend caches. Particle views are reacquired after migration.
6. `NONE`, periodic FFT, open standard Green, and open integrated Green modes keep their current input syntax and numerical behavior.
7. Existing output names, `.stat` columns such as `nBins`, and diagnostic scheduling remain stable unless a stage explicitly documents a compatibility change.
8. No hot-path allocation is introduced merely to achieve the new class structure.

Preserve the current correction distinctions during this refactor: monolithic image mode deposits primary and image charge together, binned image mode performs a separate image pass, shifted Green is binned-only, and the two correction modes remain mutually exclusive. Unifying those physics semantics would be a separate change.

## Design decisions to settle before extraction

### Common interface and context

`SelfFieldAlgorithm` has a small host-only interface:

```cpp
virtual SolverCapabilities capabilities() const = 0;
virtual void execute(SolveContext& context, SelfFieldDiagnostics& diagnostics) = 0;
```

`SolveContext` is rebuilt for every call and contains only borrowed state:

- a collection-shaped particle view with explicit primary/active selection;
- references or accessors for `R`, `P`, `Q`, `dt`, writable `E`, and writable `B`;
- invalid-particle selection plus charge, mass, conserved/dependent, and container metadata needed by the current algorithms;
- time, step, communicator, reference state, and available frame transforms;
- emission progress and the requested dynamic physics state for the call;
- writable result/diagnostic sinks.

It must not own a solver, mesh, FFT/tree state, persistent scratch, parser command object, or `PartBunch`. Attribute handles may outlive a migration, but any Kokkos views obtained from them may not. The code must explicitly reacquire device views after redistribution or particle-layout updates.

### Ownership and lifetime

- `TrackRun` owns `PartBunch` and `SelfFieldSystem`.
- Member declaration/destruction order must destroy `SelfFieldSystem` before the particle containers it borrows.
- `ParallelTracker` controls the step and constructs `SolveContext`; it does not choose a concrete solver or run 3D PIC substeps.
- `SelfFieldSystem` owns exactly one immutable algorithm selected by a factory.
- `Pic3DSolver` owns the Cartesian domain, field workspace, bin state, backend adapter, correction state, and any 3D redistribution strategy.
- `PartBunch` remains the owner of particle containers. A private 3D particle-domain adapter gives `Pic3DSolver` only the container/layout operations it needs; it does not give the solver a `PartBunch&`.

### Configuration

`SelfFieldConfigBuilder` reads `FieldSolverCmd`, `BinningCmd`, emission-source settings, and related parser objects once during setup. It produces an immutable `SelfFieldConfig` with an algorithm-specific variant. Initially the variant contains only `Pic3DConfig`; a later project can add `Pic2d5Config` without changing tracker call sites.

Static configuration belongs in the config object. Step-dependent facts, such as whether a correction is active for the current source and the current reference/mirror state, belong in `SolveContext`. No runtime algorithm object stores parser pointers.

### Redistribution

A concrete algorithm may redistribute particles because the required partition is solver-specific: 3D PIC uses a Cartesian field decomposition/ORB path, while a future tree solver or 2.5D solver may require something different. The common interface therefore exposes no ORB or `FieldLayout` type. `PicDomainManager` decides when 3D redistribution is required and uses the private particle-domain adapter to perform it.

### Units and passes

- `IterationPlan` converts the particle selection into `SolveUnit` records. The initial plans are whole-bunch and adaptive/fixed binning.
- A `SolveUnit` describes selection plus per-unit values such as mean momentum and $\gamma$; it does not own particles or fields.
- `CorrectionPlan` expands one unit into ordered `SolvePass` records: primary, image, or shifted-Green as applicable.
- A `SolvePass` is a host-side value record describing deposit selection, charge/sign policy, backend options, mirror/composition policy, and diagnostic labels. It is not a virtual object used inside a device kernel.

This makes binned and corrected paths explicit without putting binning or cathode concepts into the common algorithm interface.

### Diagnostics

The stable lifecycle events are:

1. self-field solve;
2. domain update/redistribution;
3. solve unit;
4. solve pass;
5. backend solve;
6. field composition.

Kernel-level timers can remain named `IpplTimings` entries. Do not make every kernel callback part of a permanent public diagnostics interface.

## Coding and GPU rules from `PartBunch/Binning`

The new module should follow the strongest recurring Binning conventions:

- Put the subsystem in one named namespace, tentatively `opalx::spacecharge`, instead of adding more global names.
- Put declarations in `.h`, template/device implementation in `.tpp`, and include the `.tpp` at the bottom of the header. Use `.cpp` for non-template host orchestration and explicit instantiations.
- Keep the public API before private state, use local `using` aliases near the class top, and suffix data members with `_m`.
- Add concise Doxygen comments that state ownership, invalidation, coordinate frame, and host/device behavior.
- Give every Kokkos kernel a stable descriptive label.
- Never capture `this`, a virtual object, `PartBunch`, `std::shared_ptr`, or a parser object in a device lambda. Copy the needed views and scalar/policy values to local variables first.
- Keep device policies and records trivially copyable. Put reusable device helpers in `static KOKKOS_INLINE_FUNCTION` functions. CUDA-visible enclosing functions must be public where required by the repository's CUDA constraint.
- Make host/device transfers explicit with `DualView` state, `modify_host`/`modify_device`, `sync_host`/`sync_device`, host mirrors, and `deep_copy`. Never dereference a device view on the host.
- Use `OPALX_DEVICE_COMPILATION`, Kokkos backend macros, and `if constexpr` for genuinely different device/host implementations. Do not hide unsupported code behind a runtime branch that a GPU compiler still instantiates.
- Preserve persistent IPPL/Kokkos fields and reusable scratch in `PicWorkspace`; avoid per-step field allocation.
- Preserve the current device-resident/GPU-aware MPI mirror path. Do not replace it with unconditional host staging.
- Keep temporary particle mutations contained and ordered. The current deposit/correction paths temporarily scale `dt` or modify `R`/`Q`; the extracted operation must restore them before any other consumer can observe the attributes.
- Balance timers on every early return and throw `OpalException` at host API boundaries rather than adding new `abort()` paths.
- Use explicit template instantiation for the production `double, 3` PIC path where that reduces compile cost.

The Binning directory is a style guide, not a requirement to preserve every old comment or abort path found there.

## Handling the earlier refactor attempt

Do not cherry-pick the `field-solver-refactor-staged` stack wholesale. It is far behind the current base and its top-level design still routes a driver through both `BinnedFieldSolver` and `PartBunch`, keeps OPALX `FieldSolver` derived from the IPPL solver base, and leaves frame/domain/ORB orchestration in the tracker. Moving those classes into a `SpaceCharge` directory did not fix ownership.

Ideas worth porting manually, only when their corresponding stage is reached:

- diagnostics extraction;
- typed solve requests and backend capability records;
- immutable configuration snapshots;
- correction-pass value descriptions;
- field-operation/Lorentz accumulation helpers;
- the local `SpaceCharge` CMake layout.

Treat that branch as design evidence and a source of small leaf-level implementations, not as a commit sequence or target architecture.

## Implementation stages

Each stage is intended to be a reviewable change with a production build and the named regression gate. Do not combine stages merely to reduce commit count.

### Stage 0 - Freeze the baseline and add a safe smoke runner

Work:

- Configure a worktree-specific Release/SERIAL build with unit tests disabled.
- Run the baseline smoke matrix before changing code and retain generated `.stat` files under the build tree.
- Capture a two-rank SwissFEL baseline separately. Later two-rank results may be compared to this baseline, not silently to the checked-in one-rank reference.
- Add a small runner under `tools/` that copies an existing regression-test directory into `build_refactor_take_two/smoke/<test>/`, runs there, and applies the same `.rt` comparison logic as `.github/workflows/cpu-serial.yml`.
- The runner must never write into `RegressionTests/` or `reference/`, and must support an exit-only `--np 2` mode for redistribution checks.
- Record the input commit, OPALX executable path, MPI rank count, and command line beside each baseline result.

Gate: build plus the fast gate (`Drift-1-fromfile`, `Drift-3-open-fromfile`, `Drift-3-periodic-fromfile`, `Drift-4-open-bins-fromfile`, and `Drift-4-multi-emit-open`) and the two-rank SwissFEL baseline.

### Stage 1 - Introduce common value contracts and immutable configuration

Add the initial `src/SpaceCharge/` module and CMake entry with:

- `SolverCapabilities.h`
- `SelfFieldConfig.h`
- `SelfFieldConfigBuilder.h/.cpp`
- `ParticleSetView.h`
- `SolveContext.h`
- `SelfFieldDiagnostics.h`

Work:

- Separate static input configuration from per-call state.
- Make particle/container selection explicit and collection-shaped while preserving today's primary-container behavior.
- Document view invalidation and writable attributes in the types themselves.
- Add an algorithm-specific config variant with only `Pic3DConfig` for now.
- Keep all types host-side unless a small record is deliberately copied into a Kokkos kernel.

Gate:

- Production build only; no runtime behavior changes.
- A dependency scan confirms common headers do not include `PartBunch`, IPPL `Field`, `FieldLayout`, ORB, `BinnedFieldSolver`, or 3D correction headers.

### Stage 2 - Add `SelfFieldSystem` and cut over the call site through a legacy bridge

Add:

- `SelfFieldAlgorithm.h`
- `SelfFieldSystem.h/.cpp`
- `SelfFieldFactory.h/.cpp`
- a clearly temporary `LegacyPic3DAlgorithm` implementation.

Work:

- Make `TrackRun` own the system with correct destruction order.
- Make `ParallelTracker` construct one `SolveContext` and call `SelfFieldSystem::solve`.
- Initially let the private legacy algorithm invoke the current `PartBunch`/`BinnedFieldSolver` path. This temporary implementation may privately borrow legacy state, but neither the common interface nor context may expose it.
- Validate requested features against capabilities before dispatch.
- Preserve the current tracker frame/mesh order during this stage.

Gate: full baseline smoke matrix with byte-level input parity and `.rt` result parity.

### Stage 3 - Extract diagnostics without changing solver work

Work:

- Move solve counters, field-dump scheduling, bin-table reporting coordination, and lifecycle timing into `SelfFieldDiagnostics` or focused sinks owned by the system.
- Preserve current file names, dump cadence, counter increments per backend solve, and `nBins` output.
- Add scoped host-side events for solve, unit, pass, domain update, backend, and composition.
- Keep low-level Kokkos/IPPL timers where they are until their owning operation moves.

Gate: fast gate plus `Drift-4-open-bins-fromfile`; compare diagnostic file set and `nBins` history with baseline.

### Stage 4 - Isolate the IPPL Poisson backend

Add `SpaceCharge/Ippl/IpplPoissonAdapter` and host-side request/capability records.

Work:

- Move concrete IPPL solver construction, `setRhs`, `setLhs`, `solve`, solver parameter setup, and backend refresh into the adapter.
- Replace repeated string-based dispatch with construction-time selection and typed capabilities.
- Keep `NONE`, FFT, open standard, open integrated, and current error behavior for unsupported CG configurations.
- Preserve the current construction warm-up and post-warm-up solver-counter reset while ownership is still legacy; make their later destination explicit rather than losing this hidden solve during extraction.
- Give `refresh` an explicit place in the API so a domain/layout change cannot leave stale backend caches.
- Leave a registration seam for future backends, but do not import unrelated or uncommitted backend work into this refactor.

Gate: fast gate plus `Drift-3-open-integrated-fromfile`. Review standard versus integrated Green selection and FFT background handling explicitly.

### Stage 5 - Move field ownership into `PicWorkspace`

Work:

- Turn `FieldContainer` responsibilities into a 3D solver-owned `PicWorkspace`, initially with a compatibility facade for the old callers.
- Move charge, potential, electric field, accumulated electric/magnetic fields, shifted-mirror scratch, mesh, and field layout under one lifetime.
- Allocate persistent scratch on initialization/domain change, not in the per-bin/per-pass loop.
- Make coordinate frame and valid-domain assumptions explicit for every field.
- Update `BinnedFieldSolver` to borrow the workspace rather than own or discover fields indirectly.

Gate: fast gate, integrated Green, and binned drift. Inspect allocation profiles/logs to ensure no new per-step field allocation.

### Stage 6 - Extract `PicDomainManager` and solver-owned redistribution

Work:

- Move bounds calculation, emission-domain stretching, image-domain extension, longitudinal resize rules, mesh/layout rebuilding, ORB scheduling, and backend refresh coordination out of `PartBunch`/`ParallelTracker`/`LoadBalancer`.
- Introduce a private `PicParticleDomainAdapter` over borrowed particle containers. It may understand the 3D particle-layout bridge; the common interface may not.
- Preserve the exact refresh order: workspace/layout, particle layout and migration, reacquire views, backend refresh.
- Make redistribution a `Pic3DSolver` decision. Remove any `FieldSolver` reference from the load-balancing operation.
- Reacquire all Kokkos attribute views after migration before launching another kernel.

Gate:

- Fast gate and binned drift at one rank with `.rt` comparison.
- `SwissFEL-booster-SC-emittedfromfile` at two ranks as a run-to-completion redistribution smoke. Do not compare its two-rank output to a one-rank reference unless a rank-independent tolerance is first demonstrated.
- Compare its two-rank `.rt` metrics to the Stage 0 two-rank baseline where deterministic parity is available.

### Stage 7 - Extract scatter/gather and field composition

Add `PicScatterGather` and `FieldComposer`.

Work:

- Encapsulate CIC deposition/gather using the existing IPPL operations rather than reimplementing interpolation.
- Make charge weighting (`dt * Q`), selected-particle hashing/ranges, normalization, periodic background handling, and writable particle destinations explicit.
- Move Lorentz conversion, sign policies, per-unit accumulation, and final gather into `FieldComposer`.
- Use small trivially copyable device policies and named kernels. Copy views/scalars locally before capture.
- Keep whole-bunch and bin-restricted entry points on one implementation so they cannot drift apart.

Gate: fast gate, binned drift, and `Drift-4-multi-emit-open`. Review each new kernel under CUDA/HIP compilation rules even though local runtime is CPU-only.

### Stage 8 - Replace implicit bin loops with `IterationPlan`

Work:

- Add `WholeBunchPlan` and `BinningPlan`, both producing `SolveUnit` records consumed by the same execution loop.
- Reuse `ParticleBinning::AdaptBins` behavior rather than rewriting binning.
- Keep rebin, merge, sort, empty-bin handling, mean-momentum calculation, $\gamma$ calculation, longitudinal mesh stretch, and table-print cadence unchanged.
- Keep bin state solver-owned and persistent where the current algorithm requires history.
- Ensure the plan captures selections/indices, not owning particle or field objects.

Gate: `Drift-4-open-bins-fromfile`, `Drift-4-multi-emit-open`, and `AWAGun-1-emittedfromfile`; compare `nBins` over the full run, not only the final value.

### Stage 9 - Replace correction branches with `CorrectionPlan` and `SolvePass`

Work:

- Express primary, image-charge, and shifted-Green work as ordered `SolvePass` records.
- Move image scatter transforms and shifted-Green kernel/mirror options behind focused policies.
- Preserve correction mutual-exclusion rules, step budgets, cathode/mirror geometry, field signs, solve-call counts, and dump suppression.
- Preserve the existing monolithic-versus-binned image behavior and the current binned-only shifted-Green restriction. Do not use the pass model to change physical semantics.
- Preserve the existing device-resident, GPU-aware MPI field mirror and out-of-place scratch field.
- Validate correction/backend compatibility before any particle or field mutation.

Gate:

- `AWAGun-1-emittedfromfile` and `SwissFEL-booster-SC-emittedfromfile` for shifted Green; require the log/timing data to show that the shifted-Green pass actually executed.
- Repeat the SwissFEL run-to-completion smoke at two ranks to exercise the distributed mirror/redistribution combination.
- Existing regression tests do not cover the separate image-charge pass. Record that gap; do not claim image-pass parity from shifted-Green tests.

### Stage 10 - Make `Pic3DSolver` the complete 3D algorithm

Work:

- Replace the hollowed-out `BinnedFieldSolver` orchestration with `Pic3DSolver::execute`.
- Move 3D beam-frame rotation, mesh-frame policy, per-unit longitudinal stretch, field conversion, and final restoration out of `ParallelTracker` into a dedicated `FieldFramePolicy` used by `Pic3DSolver`.
- Keep the tracker in lab/reference coordinates when it calls the common interface. This is essential for a future 2.5D solver that performs its own Frenet transform.
- Preserve today's explicit primary-container 3D selection while keeping the context collection-shaped.
- Ensure every exit path restores coordinates/fields to the documented caller frame.
- Delete the temporary `LegacyPic3DAlgorithm` only after all modes use the new path.

Gate: complete solver smoke matrix at one rank plus the two-rank SwissFEL run-to-completion smoke. Compare final coordinates and fields indirectly through every available `.rt` metric.

### Stage 11 - Complete ownership cleanup and module organization

Work:

- Remove `PartBunch::setSolver`, `PartBunch::computeSelfFields`, concrete solver discovery/casts, and obsolete solver slots once no source path uses them.
- Move solver warm-up/initialization and its counter reset out of `PartBunch::pre_run` into `SelfFieldSystem`/`Pic3DSolver` initialization.
- Remove the old `FieldSolver`, `BinnedFieldSolver`, compatibility facades, and stale `LoadBalancer` coupling.
- Stop using obsolete `ippl::PicManager` solver/field/load-balancer slots once the new ownership path no longer depends on them; handle required unused hooks explicitly rather than leaving hidden discovery paths.
- Keep particle ownership and particle-statistics behavior in `PartBunch`.
- Move 3D solver files into `src/SpaceCharge/Pic3D/`; perform mechanical file moves separately from semantic changes.
- Keep `PartBunch/Binning` in place during functional extraction. A later mechanical move to `SpaceCharge/Pic3D/Binning` is optional and must not be combined with behavior changes.
- Verify no runtime algorithm holds `FieldSolverCmd*`, `BinningCmd*`, `EmissionSource*`, or `PartBunch*`.
- Run formatting only over touched files and run `git diff --check`.

Gate: complete one-rank solver matrix and the two-rank smoke. A source dependency scan must show that `ParallelTracker` sees only `SelfFieldSystem`/`SolveContext`, and common headers see no 3D-only types.

### Stage 12 - Perform the 2.5D extension audit without implementing it

Use `25dsolver/src/PartBunch/Solve2d5.h/.hpp` only as an integration reference.

Confirm that a future `Pic2d5Solver` can:

- implement `SelfFieldAlgorithm` without inheriting `Pic3DSolver` or `BinnedFieldSolver`;
- receive particle collections in lab/reference coordinates and select all required active containers;
- own its reference-path data, Frenet transform, slice mesh/layout/fields, line-density data, and 2D Poisson-solver array;
- support delayed initialization after the reference-path/orbit-threader output becomes available without adding a concrete-solver branch back to `ParallelTracker`;
- avoid any Cartesian `FieldLayout`, ORB, binning, or 3D correction type in the common interface;
- add `Pic2d5Config` to the factory/config variant without changing `ParallelTracker`;
- publish solver-specific timings beneath the stable common diagnostic lifecycle.

Document the future mapping:

| Current 2.5D prototype responsibility | Future owner |
|---|---|
| reference-path loading/interpolation | `Pic2d5Solver` reference-path component |
| lab/Frenet coordinate conversion | 2.5D frame policy |
| slice mesh, charge/field arrays, line density | 2.5D workspace |
| array of 2D open Poisson solvers | 2.5D backend component |
| scatter, slice solve, gather, lab-field restore | `Pic2d5Solver::execute` orchestration |

Do not copy the prototype's inheritance from `BinnedFieldSolver`, `PartBunch*` storage, low-level IPPL scatter/gather calls, or current single-rank layout assumption into the new common design.

Gate: architectural include/dependency audit and a short design note only. No dummy 2.5D implementation and no 2.5D test changes.

### Stage 13 - Final production verification and handoff

Work:

- Rebuild Release/SERIAL from a clean worktree-specific build directory with unit tests disabled.
- Run the complete smoke matrix and compare one-rank results to the Stage 0 baseline and checked-in references.
- Run the selected two-rank run-to-completion checks.
- Run `FodoCell-fromfile` at the final gate to retain parity with the repository's current CI smoke subset, while keeping it out of the per-stage solver loop.
- Obtain production-target compile checks for CUDA and HIP with unit tests disabled. Add SYCL/OpenMP compile checks when runners are available.
- Review every device lambda, host/device sync boundary, field-layout refresh, and MPI buffer path.
- Produce a compatibility report listing preserved inputs/outputs and the known image-charge regression coverage gap.
- Leave unit-test migration as a clearly separate follow-up before the branch is considered merge-ready under the repository's current required CI.

Gate: no production build failures, all selected regression comparisons pass, GPU production targets compile, no old orchestration path remains, and the 2.5D audit passes.

## Build and regression commands

Use a build directory dedicated to this worktree:

```bash
cmake -S /Users/aliemen/Desktop/opalx/refactor-take-two \
      -B /Users/aliemen/Desktop/opalx/build_refactor_take_two \
      -DBUILD_TYPE=Release \
      -DPLATFORMS=SERIAL \
      -DOPALX_ENABLE_UNIT_TESTS=OFF \
      -DCMAKE_C_COMPILER=/opt/homebrew/bin/clang \
      -DCMAKE_CXX_COMPILER=/opt/homebrew/bin/clang++
cmake --build /Users/aliemen/Desktop/opalx/build_refactor_take_two -j8
```

The smoke runner introduced in Stage 0 should invoke:

```bash
mpirun -np 1 /Users/aliemen/Desktop/opalx/build_refactor_take_two/src/opalx \
    <test>.in --info 2
```

from a copied test directory under the build tree, then apply the `.rt` comparison algorithm already present in `.github/workflows/cpu-serial.yml`. Never execute in the existing regression-test directories because that checkout already contains unrelated generated outputs that must be preserved.

## Regression matrix

| Test | Coverage | Use |
|---|---|---|
| `Drift-1-fromfile` | `TYPE=NONE`, tracker/no-self-field baseline | every stage |
| `Drift-3-open-fromfile` | whole-bunch OPEN, standard Green | every behavior-changing stage |
| `Drift-3-periodic-fromfile` | periodic FFT and background handling | backend/domain/full gates |
| `Drift-3-open-integrated-fromfile` | OPEN integrated Green | backend/full gates; not proof of shifted-Green correction parity because it has no binning |
| `Drift-4-open-bins-fromfile` | adaptive binning, acceleration, per-bin Lorentz accumulation | diagnostics/workspace/bin/full gates |
| `Drift-4-multi-emit-open` | open solver, binning, continuous plus delayed emission source | scatter/gather/bin/full gates and existing CI parity |
| `AWAGun-1-emittedfromfile` | adaptive binning, emission, shifted Green | bin/correction/full gates |
| `SwissFEL-booster-SC-emittedfromfile` | shifted Green, emission, frequent repartition request | correction/full gates; two-rank smoke compared to its Stage 0 two-rank baseline where deterministic |
| `FodoCell-fromfile` | current CI `TYPE=NONE` tracker/element coverage | final gate only |

The existing suite has no dedicated image-charge-pass case and no solver case for P3M, CG, or Barnes-Hut. Those are recorded verification gaps, not a reason to fold new regression inputs into this refactor. `AWAGun-1-emittedfromfile` requests repartition less frequently than its total step count, so it covers shifted Green but not ORB; SwissFEL is the redistribution test.

## GPU acceptance checklist

At every stage that adds or changes a kernel:

1. No `this` or host-only object capture in `KOKKOS_LAMBDA`.
2. All captured records are trivially copyable and all views refer to the correct execution memory space.
3. Host inspection uses an explicit mirror/sync.
4. Particle views are reacquired after migration/reallocation.
5. CUDA-visible callers have the required public access.
6. Backend-specific code is excluded at compile time, not merely by a runtime condition.
7. MPI send/receive buffers preserve the current device-aware path and lifetime until completion.
8. Kernel labels and timers are present and balanced.
9. Persistent workspace is reused across units/passes/steps.
10. CUDA and HIP compile the production `opalx` target before the stage is accepted.

The Mac can validate SERIAL runtime semantics only. A local pass is never evidence that the new code is GPU-safe.

## Deferred work and explicit non-goals

- Adapting, moving, or adding unit tests.
- Implementing or merging the 2.5D solver.
- Fixing the 2.5D prototype's current multi-rank limitation.
- Adding a new image-charge regression input.
- Changing solver mathematics, input syntax, physical defaults, or numerical tolerances.
- Importing unrelated uncommitted work from the original `opalx` worktree.
- Redesigning IPPL field/scatter/gather APIs.
- Moving Binning files while semantic extraction is still in progress.

## Completion criteria

The refactor is complete for this pass when:

- `ParallelTracker` only builds context and invokes `SelfFieldSystem`;
- `PartBunch` owns particles but not solver orchestration or fields;
- `Pic3DSolver` owns the complete current 3D PIC algorithm and its domain/workspace/backend state;
- binning and corrections are plans producing explicit units/passes;
- layout changes have one refresh path and no stale particle/device views;
- common headers contain no 3D/IPPL/ORB dependencies;
- all selected one-rank regressions match, selected two-rank smokes finish, and production CUDA/HIP targets compile;
- no 2.5D code has been added, but the audited extension path requires a new concrete algorithm rather than another refactor of the common boundary;
- deferred unit-test work and the image-charge coverage gap are recorded plainly.
