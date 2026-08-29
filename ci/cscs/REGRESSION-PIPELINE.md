# OPALX Regression Test Pipeline

This document describes the CSCS CI/CD regression test pipeline for OPALX.
It covers the architecture, file layout, every step in the pipeline, key
workarounds, and how to resume work or debug failures.

## Overview

The regression pipeline runs **after the build and unit-test stages** on
every PR and on master. It uses the pre-built **release** artifact, clones
the NightlyBuildX scripts and regression-tests-x repo, and runs the full
regression test suite via the NightlyBuildX `run_tests` driver.

Three architectures are supported:

| CI target | Architecture string | Runner | UENV |
|-----------|--------------------|--------|------|
| OpenMP (Eiger Zen2) | `cpu-zen2-openmp` | `.baremetal-runner-eiger-zen2` | `prgenv-gnu/25.6:v2` |
| GH200 (Daint) | `gpu-cuda-gh200` | `.baremetal-runner-daint-gh200` | CUDA 13.1.1 squashfs |
| MI300 (Beverin) | `gpu-rocm-mi300` | `.baremetal-runner-beverin-mi300` | `prgenv-gnu/25.07-6.3.3:v12` |

## File Layout

```
ci/cscs/
├── common.yml                        # Stages + PR number extractor
├── cscs-openmp.yml                   # Entry point: includes build + run + regression
├── cscs-gh200.yml                    # Entry point: includes build + run + regression
├── cscs-mi300.yml                    # Entry point: includes build + run + regression
├── regression-common.yml             # Shared regression job template (.opalx-regression-common)
├── regression-test-mpirun-wrapper    # mpirun → srun translation script
├── config/
│   ├── openmp.conf                   # NightlyBuildX config for OpenMP
│   ├── cuda.conf                     # NightlyBuildX config for CUDA/GH200
│   └── rocm.conf                     # NightlyBuildX config for ROCm/MI300
├── openmp/
│   ├── build_openmp.yml              # Build stage (pre-existing)
│   ├── run_openmp.yml                # Unit test stage (pre-existing)
│   └── regression_openmp.yml         # Regression job definition
├── cuda/
│   ├── build_sm90.yml                # Build stage (pre-existing)
│   ├── run_sm90.yml                  # Unit test stage (pre-existing)
│   └── regression_sm90.yml           # Regression job definition
└── rocm/
    ├── build_rocm-6.3.yml            # Build stage (pre-existing)
    ├── run_rocm-6.3.yml              # Unit test stage (pre-existing)
    └── regression_rocm-6.3.yml       # Regression job definition
```

## Stage Ordering

Defined in `common.yml`:

```yaml
stages:
  - opalx_build
  - opalx_test
  - opalx_regression
```

Regression jobs depend only on the release build (via `needs:`), so they
can start as soon as the build finishes, even if unit tests are still
running.

## Entry Point Includes

Each architecture entry point includes three files:

```yaml
# cscs-openmp.yml
include:
  - local: "ci/cscs/openmp/build_openmp.yml"
  - local: "ci/cscs/openmp/run_openmp.yml"
  - local: "ci/cscs/openmp/regression_openmp.yml"
```

The CSCS admin console points at these entry point files. The names
`cscs-openmp`, `cscs-gh200`, and `cscs-mi300` are the pipeline names used
with `cscs-ci run`.

## Architecture-Specific Regression Files

Each `regression_<arch>.yml` defines:

- `BUILD_ARCH` — the architecture string used for directories and display
- `REGTEST_CONFIG` — the NightlyBuildX config file under `ci/cscs/config/`
- `CTEST_SITE` — diagnostic site name
- `SRUN_ARGS` — srun arguments (uenv, view, repo, GPU flags)
- `WRAPPER` — optional MPI wrapper binary path (CUDA only)
- The runner to use (`.baremetal-runner-*`)
- `SLURM_CPUS_PER_TASK` per architecture
- A single job (e.g. `opalx-regression-openmp-release`) that:
  - Extends the common template + arch runner
  - `needs: ["opalx-build-<arch>-release"]`
  - Sets `BUILD_DIR: "build-$CI_COMMIT_SHORT_SHA-release"`

## The mpirun → srun Wrapper

**File:** `ci/cscs/regression-test-mpirun-wrapper`

The regression test `.local` scripts (committed in regression-tests-x) call
`mpirun -np N <opalx> <args>`. On CSCS systems we must use `srun` with uenv
flags instead. This wrapper intercepts `mpirun`, extracts `-np N`, and
translates it to:

```bash
srun ${SRUN_ARGS} ${WRAPPER} env MALLOC_CHECK_=0 GLIBC_TUNABLES=glibc.malloc.tcache_count=0 -n N <args>
```

The wrapper is copied to `ci-bin/mpirun` and prepended to `PATH` so that
all `.local` scripts use it transparently.

## Step-by-Step Pipeline Execution

The shared template `.opalx-regression-common` in `regression-common.yml`
executes the following steps:

### 1. PR Number Extraction

Uses the `.pr-number-extractor` from `common.yml` to set `CDASH_LABEL`
(e.g., `PR-508`).

### 2. Build Path Setup

```bash
export BUILD_PATH=$(pwd)/$BUILD_DIR
# e.g. BUILD_PATH=/capstor/.../build-89e66c54-release
```

### 3. Publish Branch Normalization

CSCS CI creates merge pipelines with refs like `__CSCSCI__pr508`. These are
normalized:

- `__CSCSCI__pr527` → `pr-527`
- `master` → `master`
- anything else → empty (skip push)

### 4. PR Title Fetch from GitHub

If `OPALX_PUBLISH_BRANCH` matches `pr-<number>`, the pipeline fetches the PR
title from the GitHub API (no token needed, repos are public):

```bash
curl -fsSL "https://api.github.com/repos/OPALX-project/OPALX/pulls/${PR_NUM}" \
  | python3 -c "import json,sys; print(json.load(sys.stdin)['title'])"
```

The result is appended: `pr-527 Fix remaining OPALX and external dependency warnings`.
This full string is passed to `run_tests --opalx-branch="..."`. NightlyBuildX
sanitizes it for directory names but writes the full string to
`branch-name.txt` for display in the overview HTML.

### 5. Verify Release Build Artifact

```bash
ls -l $BUILD_PATH/bin/
export LD_LIBRARY_PATH=$BUILD_PATH/lib:$LD_LIBRARY_PATH
export MALLOC_CHECK_=0
```

### 6. Fix ctest Paths (Build → Test Directory)

The build artifacts contain `CTestTestfile.cmake` and `DartConfiguration.tcl`
with hardcoded paths from the build job's working directory (e.g.,
`f7t/1/.../build-89e66c54-release`). The test job runs in a different
directory (e.g., `f7t/6/...`). Without fixing, ctest uses stale paths and
test executables can't find `libippl.so`.

```bash
CACHE_DIR=$(grep "CMAKE_CACHEFILE_DIR:INTERNAL=" $BUILD_PATH/CMakeCache.txt | cut -d'=' -f2)
find $BUILD_PATH -type f -name "CTestTestfile.cmake" -exec sed -i "s|$CACHE_DIR|$BUILD_PATH|g" {} \;
find $BUILD_PATH -type f -name "DartConfiguration.tcl" -exec sed -i "s|$CACHE_DIR|$BUILD_PATH|g" {} \;
```

This is the same fix used by the regular unit-test pipeline
(`run_openmp.yml`, `run_sm90.yml`, `run_rocm-6.3.yml`).

### 7. Install NightlyBuildX Tools (matplotlib, pandoc, xsltproc, perl)

A persistent conda environment is created at
`$SCRATCH/nightlybuildx-tools-$(uname -m)` using micromamba. The check
verifies not just that binaries exist, but that Python can actually import
`site` (to catch corrupted environments from uenv updates):

```bash
if [ ! -x "$TOOLS_ENV/bin/python3" ] || ! "$TOOLS_ENV/bin/python3" -c "import site" 2>/dev/null || ...; then
    rm -rf "$TOOLS_ENV"
    # download micromamba, create environment with python matplotlib pandoc libxslt perl
fi
```

If the environment is broken, it is deleted and recreated from scratch.

### 8. Clone NightlyBuildX

```bash
git clone --depth 1 --branch alps-regression-testing https://github.com/OPALX-project/NightlyBuildX.git nightlybuildx
```

The `alps-regression-testing` branch is a work-in-progress branch that adds
ALPS/CSCS-friendly options to NightlyBuildX:

- `--no-salloc` — skip the `salloc` wrapper when already inside a SLURM
  allocation.
- Standard-folder aware `--opalx-build-dir` — detects `bin/opalx` as well as
  `src/opalx`.
- Robust OPALX revision detection via `OPALX_SRC_DIR`.
- `run-reg-tests.py` exits non-zero when tests fail.

Once verified, this branch will be proposed back to NightlyBuildX via a PR.

### 9. Pass `--no-salloc` to `run_tests`

The NightlyBuildX `run_tests` script tries to use `salloc` when running
locally, but inside CSCS CI the allocation already exists. The `--no-salloc`
flag disables the salloc wrapper.

### 10. Provide `OPALX_SRC_DIR` for Revision Detection

Because the regression runner uses a pre-built binary from a GitLab artifact,
the OPALX source directory is not in the fixed relative location that
NightlyBuildX expects. Setting `OPALX_SRC_DIR=$CI_PROJECT_DIR` lets it read the
git revision from the actual checkout.

### 11. Install mpirun → srun Wrapper

```bash
cp $CI_PROJECT_DIR/ci/cscs/regression-test-mpirun-wrapper ci-bin/mpirun
chmod +x ci-bin/mpirun
export PATH=$(pwd)/ci-bin:$PATH
```

### 12. Run Regression and Unit Tests

```bash
OPALX_SRC_DIR=$CI_PROJECT_DIR \
bash nightlybuildx/scripts/run_tests \
  --config=$CI_PROJECT_DIR/ci/cscs/config/$REGTEST_CONFIG \
  --opalx-exe-path=$BUILD_PATH/bin \
  --opalx-build-dir=$BUILD_PATH \
  --opalx-branch="${OPALX_PUBLISH_BRANCH:-$CI_COMMIT_REF_NAME}" \
  --no-salloc \
  --reg-tests \
  --unit-tests \
  --no-gpl \
  --publish-dir=$CI_PROJECT_DIR/regression-results-$BUILD_ARCH
```

This runs both unit tests (`ctest -L unit`) and regression tests
(`run-reg-tests.py`) using the pre-built opalx binary. The `--no-gpl` flag
selects matplotlib instead of gnuplot for comparison plots.

### 13. Collect Artifacts

```yaml
artifacts:
  paths:
    - regression-results-$BUILD_ARCH
  expire_in: 1 week
```

### 14. Push Results to opal-live-doc

If `PSI_GIT_SSHKEY` is set and `OPALX_PUBLISH_BRANCH` is non-empty, the
pipeline pushes results to `git@gitea.psi.ch:AMAS/opal-live-doc.git`:

1. Reconstruct the SSH private key from the `PSI_GIT_SSHKEY` environment
   variable (which may be base64 or PEM with whitespace removed).
2. Clone `opal-live-doc` with sparse checkout of `docs/opalx-regression-test/`.
3. Copy regression results into the repo.
4. Re-render the overview HTML using `run_tests --do-not-compile-run`.
5. Commit and push to `master`, with up to 5 retry attempts to handle
   concurrent pushes from different architectures.

If `PSI_GIT_SSHKEY` is not set, the push is skipped (useful for testing).

## NightlyBuildX Config Files

Located at `ci/cscs/config/`:

```bash
# openmp.conf
branch="master"
regtests_branch="master"
architecture="cpu-zen2-openmp"
do_unittests='yes'
```

```bash
# cuda.conf
branch="master"
regtests_branch="master"
architecture="gpu-cuda-gh200"
do_unittests='yes'
```

```bash
# rocm.conf
branch="master"
regtests_branch="master"
architecture="gpu-rocm-mi300"
do_unittests='yes'
```

## Key Workarounds and Their Reasons

### mpirun → srun wrapper

The `.local` scripts in regression-tests-x call `mpirun -np N`. On CSCS we
must use `srun` with uenv flags. The wrapper translates one to the other
without modifying the regression-tests-x repo.

### ctest path fixing

Build artifacts contain hardcoded paths from the build runner's working
directory. The test runner uses a different path. Without `sed` rewriting,
test executables launched via `srun` can't find shared libraries.

### `--no-salloc`

NightlyBuildX `run_tests` tries to allocate via `salloc`, but inside CSCS
CI the allocation already exists. The `--no-salloc` flag disables this path.

### `run-reg-tests.py` exit code

NightlyBuildX `run-reg-tests.py` exits non-zero when
`totalNrPassed != totalNrTests`, so CI fails when regression tests fail.

### `OPALX_SRC_DIR`

The pre-built binary is not located next to the OPALX source tree, so the
revision detection heuristic fails. `OPALX_SRC_DIR=$CI_PROJECT_DIR` points it
at the actual GitLab checkout.

### NightlyBuildX tools environment validation

A cached conda environment can become broken when the uenv is updated. The
check runs `python3 -c "import site"` to verify the environment is
functional, not just that the binary exists.

## External Repositories

| Repository | Branch | Purpose |
|-----------|--------|---------|
| `https://github.com/OPALX-project/NightlyBuildX.git` | `alps-regression-testing` | Scripts: `run_tests`, `run-reg-tests.py`, `OpalRegressionTests/`, HTML templates |
| `https://github.com/OPALX-project/regression-tests-x.git` | `master` | Regression test input files, `.local` scripts, reference data |
| `git@gitea.psi.ch:AMAS/opal-live-doc.git` | `master` | Published results website (GitLab Pages at `amas.pages.psi.ch`) |

## Debugging Tips

1. **Check job logs for the echo markers** — every step starts with
   `=== <step name> ===` to quickly locate where failures occur.

2. **`libippl.so: cannot open shared object file`** — the ctest path fix
   (step 6) is missing or failed. Check that `CACHE_DIR` was correctly
   extracted from `CMakeCache.txt`.

3. **`Fatal Python error: init_import_site` / `unknown encoding`** — the
   cached NightlyBuildX tools environment is corrupted. Delete
   `$SCRATCH/nightlybuildx-tools-$(uname -m)` or let the CI recreate it
   (the check should catch this automatically).

4. **`cmake` binary not found in tests** — this was a build/test uenv mount
   point mismatch. Fixed in `demos/alpine/CMakeLists.txt` by replacing
   `${CMAKE_COMMAND} -P` with `bash -c` for the `LandauDampingInit` fixture.

5. **YAML parse errors in GitLab** — ensure `echo` lines containing `: `
   are single-quoted, and multi-line commands use `>-` block scalars.

6. **Push failures to opal-live-doc** — the pipeline retries up to 5 times
   with increasing backoff. If all fail, check the SSH key in
   `PSI_GIT_SSHKEY` and that `gitea.psi.ch` is reachable.

## Remaining TODOs

- [ ] Switch regression jobs from all PRs to **master-only** once fully verified
- [ ] Open a PR from `alps-regression-testing` back to NightlyBuildX main
- [ ] Increase `SLURM_TIMELIMIT` from 30 min if timeouts occur
- [ ] Review concurrency/resource usage of regression jobs if needed
- [ ] Ensure `PSI_GIT_SSHKEY` is set as a CI/CD variable for the push step
