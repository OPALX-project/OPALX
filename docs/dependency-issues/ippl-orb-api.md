# OPALX / IPPL ORB API mismatch

Recorded on 2026-10-16 while building the `feature/parser-only` branch. This note
preserves evidence for a possible upstream issue. No issue has been submitted.

Suggested issue title: **Default IPPL master lacks the ORB overload used by
CartesianDomainUpdater; replacement overload has different axis-mask semantics**.

## Revisions and environment

| Component | Observed configuration |
| --- | --- |
| OPALX base revision | `1df5352d42cdce661010dc33b136bd2d30eb280c` |
| OPALX local changes | Uncommitted parser-only feature; the failing CartesianDomainUpdater caller is unchanged |
| IPPL checkout (`master`) | `a7e505730a76dd1edc7d0c8b96fb4223d9f84cfb` |
| IPPL ORB file source | `40a0e43257db0ecbe5c9e6f154feeea724224cf4`, branch `559-orb-does-not-allow-arbitrary-decomposition-axis` |
| Build | CMake, `Debug`, `PLATFORMS=SERIAL`, unit tests enabled |
| Compiler | GNU C++ 13.3.0, Ubuntu 24.04 toolchain |
| Runtime | Linux 5.15.167.4 WSL2; Open MPI / OpenRTE 4.1.6; Kokkos 5.2.0 |

OPALX's `cmake/Dependencies.cmake` defaults `IPPL_GIT_TAG` to the moving branch
`master`. These observations apply to the exact revisions above, not necessarily
to a later version of that branch.

## Build failure

The unchanged call in
[`CartesianDomainUpdater.cpp`](../../src/SpaceCharge/CartesianPIC3D/CartesianDomainUpdater.cpp)
is:

```cpp
orb_m.binaryRepartition(primary().R, layout, false, orbForbiddenAxes(layout))
```

The fetched IPPL checkout declares only the three-argument overload. Building
`opalx_exe` fails with a diagnostic of this form:

```text
CartesianDomainUpdater.cpp:291:37: error: no matching function for call to
  ... binaryRepartition(..., ..., bool, std::array<bool, 3>)
candidate expects 3 arguments, 4 provided
```

Expected behavior: the configured IPPL dependency provides an API compatible
with OPALX's caller, so the standard serial build compiles.

The failure was observed with local parser-only changes present. An independent
clean build of the unmodified OPALX base revision was not performed. The failing
call and the default dependency selection are outside the parser-only changes.

## Reproduction recipe for an upstream investigation

In a separate checkout of the OPALX base revision, configure a fresh build using
the same toolchain and installed prerequisites:

```bash
cmake -S . -B build-orb-repro \
  -DCMAKE_BUILD_TYPE=Debug \
  -DPLATFORMS=SERIAL \
  -DOPALX_ENABLE_UNIT_TESTS=ON \
  -DIPPL_GIT_TAG=a7e505730a76dd1edc7d0c8b96fb4223d9f84cfb
cmake --build build-orb-repro --target opalx_exe -j 4
```

The explicit IPPL commit avoids a later `master` hiding the reported mismatch.
This fresh-checkout recipe is provided for reproduction; it has not itself been
executed as part of this investigation.

## Local build workaround actually used

The final local dependency checkout remains IPPL `master` at `a7e5057...`.
Only these two files under `build_serial/_deps/ippl-src` were replaced with their
contents from IPPL commit `40a0e43...`:

- `src/Decomposition/OrthogonalRecursiveBisection.h`
- `src/Decomposition/OrthogonalRecursiveBisection.hpp`

An exact snapshot of the resulting diff against `a7e5057...` is preserved in
[`ippl-orb-local.patch`](ippl-orb-local.patch). It is an investigation artifact,
not a proposed production fix.

The build cache currently uses:

```text
IPPL_GIT_TAG:STRING=master
FETCHCONTENT_UPDATES_DISCONNECTED:BOOL=ON
```

The disconnected-update setting prevents automatic dependency fetches during
this local work. It does not make the two replacement files part of OPALX's
source configuration. They live in an ignored build directory, will not travel
with a feature commit, and can be lost when that dependency checkout is recreated
or changed. A fresh build still needs a compatible dependency configuration.

Selecting the entire IPPL `559-...` branch was also attempted. It exposed
include-order/circular-include compilation errors, including missing
`ippl::parallel_for` / `RangePolicy` declarations. That branch was not retained
as the local dependency checkout. Copying its two ORB files onto `master` is the
configuration on which the tests below passed.

## Axis-mask contract still needs resolution

**The local workaround is only validated for the tested single-rank serial use.**

OPALX's helper says the mask identifies axes ORB must **not** cut, and returns:

```cpp
const std::array<bool, 3> parallel = layout.isParallel();
return {!parallel[0], !parallel[1], !parallel[2]};
```

The copied IPPL overload names that parameter `allowedAxes`, documents `true`
as permission to cut an axis, and skips an axis when its mask entry is false.
These are opposite mask conventions. Merely adding the overload resolves the
compiler error but does not establish correct multi-rank repartitioning.

The API and mask contract should be agreed between OPALX and IPPL before this
workaround is treated as a general fix. A durable resolution should pin or
otherwise select a compatible IPPL revision and exercise distributed-axis
restrictions with multiple MPI ranks.

## Verification completed with the local workaround

All runs below used the serial Debug build and one MPI rank, with normal
simulation mode for the regression cases:

| Check | Result |
| --- | --- |
| Full CTest suite | 127/127 entries passed: 126 existing entries plus the parser-only suite; about 110 seconds |
| Parser-only integration suite | 27/27 cases passed, included in the CTest result above |
| `Drift-1-fromfile` | 7/7 reference assertions passed; 100 recorded steps |
| `Drift-4-multi-emit-open` | 9/9 reference assertions passed; 92 recorded steps |
| `FodoCell-fromfile` | 7/7 reference assertions passed; 2,001 recorded steps |

Regression inputs were from `OPALX-project/regression-tests-x`, revision
`222a048` (the three cases selected by OPALX's serial CI). The existing
`tools/run_regression_smoke.py` runner compared outputs to the supplied `.rt`
assertions and reference statistics without modifying the input cases.
The multi-emission case first exceeded a 180-second allowance; its rerun with a
600-second allowance completed and passed. The timeout was not a failed
reference comparison.

These results do not validate multi-rank ORB semantics, GPU builds, or an
unpatched fresh dependency checkout. In particular, a passing single-rank
`TestCartesianDomainUpdaterORB` entry does not establish multi-rank coverage;
that test has a path conditional on the MPI size being at least two.

Local evidence retained during this session:

- `/tmp/opalx-full-ctest.log`
- `/tmp/opalx-regression-Drift-1-fromfile.log`
- `/tmp/opalx-regression-Drift-4-multi-emit-open-long.log`
- `/tmp/opalx-regression-FodoCell-fromfile.log`
- `build_serial/Testing/Temporary/LastTest.log`
- `build_serial/smoke/<case>/parser-only-regression*/` (inputs, outputs, logs, and run manifests)

The `/tmp` and build-tree artifacts are local and temporary; retain copies if
attaching them to an upstream issue. The report and patch snapshot are the
source-tree records intended to survive build-directory cleanup.
