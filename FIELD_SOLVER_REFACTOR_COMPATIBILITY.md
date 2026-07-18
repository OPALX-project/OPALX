# Field Solver Refactor Compatibility Report

## Status

All implementation stages and the local Stage 13 verification are complete in the
`refactor-take-two` worktree. Release/SERIAL production builds, installation, the selected
one-rank regression matrix, retained Stage 0 comparisons, and selected two-rank runs pass.

CUDA and HIP production compiles cannot be performed on this Apple Silicon Mac. Static review
found no obvious source-level CUDA/HIP blocker, but the branch still requires production-only
compiles in Linux CUDA and HIP environments before it can satisfy the GPU part of the production
gate. Unit-test migration remains the explicitly deferred follow-up requested for this refactor.

## Verification configuration

- Source worktree: `/Users/aliemen/Desktop/opalx/refactor-take-two`
- Build tree: `/Users/aliemen/Desktop/opalx/build_refactor_take_two_final`
- Install tree: `/private/tmp/opalx-refactor-take-two-install`
- Configuration: `Release`, `SERIAL`, `OPALX_ENABLE_UNIT_TESTS=OFF`
- Compilers: `/opt/homebrew/bin/clang` and `/opt/homebrew/bin/clang++`
- IPPL: fetched from current `master`, commit `2180f35b0771b2c0a2694a20a9f7fd576ced242c`
- Production target: `/Users/aliemen/Desktop/opalx/build_refactor_take_two_final/src/opalx`

The build directory did not exist before the final configure. A complete build and a subsequent
incremental build after the final ownership cleanup both completed successfully. `cmake --install`
also completed successfully after configuring the temporary prefix. The installed tree contains
the common self-field API, the IPPL adapter API, the Pic3D API and template headers, including
`SpaceCharge/Pic3D/FieldMirror.hpp`, and the production executable.

## Preserved input and output contracts

No OPALX input syntax was changed. The existing parser and command objects still accept the
checked-in configurations used by the matrix, including:

- no self field (`TYPE=NONE`);
- periodic and open FFT solves;
- standard and integrated Green functions;
- adaptive binning and delayed emission;
- shifted-Green correction;
- the existing `ZEROFACE_R0Z` and `ZEROFACE_MAXSTEPS` image-charge options.

The regression runs preserve the existing `.stat`, H5, timing, element-position, and auxiliary
output names. All checked-in `.rt` assertions pass without changing tolerances or references. The
`nBins` column remains present and its complete sequence is unchanged in the retained comparisons.

Solver diagnostics retain the historical warm-up distinction: the low-level `Solve` timer has one
more call than the tracking-level field-solver count. Existing timing names for scatter, gather,
shifted Green, mirror communication, and redistribution remain observable.

## Final regression matrix

All paths below are relative to
`/Users/aliemen/Desktop/opalx/build_refactor_take_two_final/smoke`.

| Case | MPI ranks | Result | Final artifact |
|---|---:|---|---|
| `Drift-1-fromfile` | 1 | 7/7 checked-in assertions pass | `Drift-1-fromfile/stage13-verified-20260718T123059.977216Z` |
| `Drift-3-open-fromfile` | 1 | 8/8 checked-in assertions pass | `Drift-3-open-fromfile/stage13-verified-20260718T123100.970848Z` |
| `Drift-3-periodic-fromfile` | 1 | 8/8 checked-in assertions pass | `Drift-3-periodic-fromfile/stage13-verified-20260718T123100.857621Z` |
| `Drift-3-open-integrated-fromfile` | 1 | 8/8 checked-in assertions pass | `Drift-3-open-integrated-fromfile/stage13-verified-20260718T123101.123671Z` |
| `Drift-4-open-bins-fromfile` | 1 | 8/8 checked-in assertions pass | `Drift-4-open-bins-fromfile/stage13-verified-20260718T123120.321537Z` |
| `Drift-4-multi-emit-open` | 1 | 9/9 checked-in assertions pass | `Drift-4-multi-emit-open/stage13-verified-20260718T123114.643986Z` |
| `AWAGun-1-emittedfromfile` | 1 | 9/9 checked-in assertions pass | `AWAGun-1-emittedfromfile/stage13-verified-20260718T123140.052396Z` |
| `SwissFEL-booster-SC-emittedfromfile` | 1 | 9/9 checked-in assertions pass | `SwissFEL-booster-SC-emittedfromfile/stage13-verified-20260718T123137.639452Z` |
| `FodoCell-fromfile` | 1 | 7/7 checked-in assertions pass | `FodoCell-fromfile/stage13-verified-20260718T123115.465408Z` |
| `SwissFEL-booster-SC-emittedfromfile` | 2 | 9/9 Stage 0 two-rank assertions pass | `SwissFEL-booster-SC-emittedfromfile/stage13-verified-np2-20260718T123252.375919Z` |
| `Drift-3-open-fromfile` | 2 | run-to-completion pass | `Drift-3-open-fromfile/stage13-verified-open-np2-20260718T123249.834310Z` |
| image-configured `Drift-4-multi-emit-open` | 2 | parser/configuration run-to-completion pass | `Drift-4-multi-emit-open/stage13-verified-image-np2-20260718T123248.607618Z` |

The smoke runner copied every case into the build tree. The checked-in `RegressionTests` and
reference directories were not modified.

## Stage 0 and semantic parity evidence

The final executable was also run against the five retained one-rank Stage 0 `.stat` files. Every
applicable `.rt` comparison passes:

- `Drift-1-fromfile/stage13-verified-baseline-20260718T123156.958880Z`
- `Drift-3-open-fromfile/stage13-verified-baseline-20260718T123158.470386Z`
- `Drift-3-periodic-fromfile/stage13-verified-baseline-20260718T123200.346574Z`
- `Drift-4-open-bins-fromfile/stage13-verified-baseline-20260718T123219.459084Z`
- `Drift-4-multi-emit-open/stage13-verified-baseline-20260718T123217.737553Z`

Full `nBins` sequence comparisons are exact:

- AWAGun final versus Stage 11: 235 rows;
- SwissFEL two-rank final versus Stage 0: 344 rows;
- SwissFEL two-rank final versus Stage 11: 344 rows.

Correction and redistribution counters also retain the expected semantics:

| Evidence | Tracking calls | Low-level solves | Shifted Green | Mirror | Redistribution |
|---|---:|---:|---:|---:|---:|
| AWAGun, one rank | 5202 | 5203 | 2408 | 2408 | not exercised |
| SwissFEL, two ranks | 6824 | 6825 | 3412 | 3412 | 96 binary repartitions |

The SwissFEL run also records 96 `allReduce` and 96 `scatterR` operations, confirming that the
two-rank redistribution path executed rather than merely parsing its configuration.

## Architecture and GPU review

The final source scan finds no production references to the deleted `BinnedFieldSolver`, old
`FieldSolver`, `FieldContainer`, `LoadBalancer`, `ImageChargeScatterController`, or
`ippl::PicManager` orchestration paths. `ParallelTracker` sees only the common
`SolveContext`/`SelfFieldSystem` boundary for self-field execution, while `TrackRun` uses only
configuration, factory, and system assembly interfaces. The final cleanup also removes a
discarded per-step global-bounds calculation that previously performed two MPI reductions and a
barrier without using their results.

`PicDomainManager` now preserves layout-refresh and backend-rebind work across exceptions. It marks
the backend dirty before either an extent rebuild or ORB can mutate the layout. A later update
repeats every workspace field-layout refresh before migration and clears the state only after the
backend refresh succeeds.

Static review of every refactor-owned Kokkos kernel found:

- no `this`, parser, owning class, or `PartBunch` capture;
- public CUDA-visible enclosing methods;
- current particle and field views reacquired after migration or layout changes;
- explicit host mirrors/deep copies for host field inspection;
- persistent solver workspace and scratch allocation;
- backend destruction before the fields it borrows;
- send/receive buffers retained until MPI completion in the shifted-Green mirror path.

The existing device-resident communication design is preserved. FFT adapters set
`use_gpu_aware=true`, and `FieldMirror` communicates buffers in the field view's memory space.
CUDA and HIP runtime use therefore requires genuinely GPU-aware MPI; unconditional host staging
was deliberately not introduced.

This Mac has no `nvcc`, `hipcc`, `amdclang++`, CUDA toolkit, ROCm installation, or compatible GPU.
The repository's CUDA and HIP workflows document suitable Linux toolchain environments, but they
currently enable the deferred unit tests and cannot be used unchanged for this production-only
gate. Required external evidence is:

1. a clean `opalx` production compile with `PLATFORMS=CUDA` and unit tests disabled;
2. a clean `opalx` production compile with `PLATFORMS=HIP` and unit tests disabled;
3. a runtime shifted-Green/redistribution check on hardware with GPU-aware MPI.

## Recorded gaps and deferred follow-up

- The checked-in suite has no case that proves an active separate image-charge pass. The explicit
  two-rank smoke parses `ZEROFACE_R0Z=true` and `ZEROFACE_MAXSTEPS=7` and exits successfully, but
  its timing data does not show an image/mirror pass. It is configuration coverage only.
- Unit tests were intentionally not built or migrated. The legacy unit-test files that include
  deleted APIs must be adapted in a separate follow-up before unit-test-enabled CI can pass.
- No dedicated regression case covers CG, P3M, or Barnes-Hut solvers.
- SYCL and OpenMP production compile checks were not run in this final local SERIAL gate.
- CUDA/HIP compile and hardware runtime evidence remain external acceptance requirements.
- The checked-in CUDA/HIP workflows must either run a production-only compile with unit tests
  disabled or wait for the separate unit-test migration.

The future 2.5D integration boundary, ownership model, and prototype hazards are documented in
[`FIELD_SOLVER_2D5_EXTENSION_AUDIT.md`](FIELD_SOLVER_2D5_EXTENSION_AUDIT.md). No 2.5D code or tests
were added during this refactor.
