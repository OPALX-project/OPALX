# OPALX Regression Test Pipeline

This document describes the CSCS CI/CD regression test pipeline for OPALX: the file layout,
every step of the job, how a run is published to the results site, and how to debug failures. 

## Overview

The regression job runs **after the build stage** on every PR and on master. It uses the
pre-built **release** artifact and runs the OPALX unit tests (`ctest -L unit`) and every test of
regression-tests-x with `tools/regtest.sh` from the results repository
[AMAS/opalx-regression](https://gitea.psi.ch/AMAS/opalx-regression). Runs of master and of PRs
are pushed to that repository, whose Pages build shows them at
https://amas.pages.psi.ch/opalx-regression/results/. The repository's README describes the run
format and the tools in detail.

Three architectures are supported:

| CI target | Architecture name on the site | Runner |
|-----------|-------------------------------|--------|
| OpenMP (Eiger Zen2) | `cpu-zen2-openmp` | `.baremetal-runner-eiger-zen2` |
| GH200 (Daint) | `gpu-cuda-gh200` | `.baremetal-runner-daint-gh200` |
| MI300 (Beverin) | `gpu-rocm-mi300` | `.baremetal-runner-beverin-mi300` |

The uenv of each target is set in its `build_*.yml` and `regression_*.yml`.

## File Layout

```
ci/cscs/
├── common.yml                        # Stages + PR number extractor
├── cscs-openmp.yml                   # Entry point: includes build + run + regression
├── cscs-gh200.yml                    # Entry point: includes build + run + regression
├── cscs-mi300.yml                    # Entry point: includes build + run + regression
├── regression-common.yml             # Shared regression job template (.opalx-regression-common)
├── regression-test-mpirun-wrapper    # Launcher: "-np N" -> srun with the job's uenv
├── openmp/
│   ├── build_openmp.yml              # Build stage
│   ├── run_openmp.yml                # Unit test stage (CDash)
│   └── regression_openmp.yml         # Regression job definition
├── cuda/
│   ├── build_sm90.yml
│   ├── run_sm90.yml
│   └── regression_sm90.yml
└── rocm/
    ├── build_rocm-6.3.yml
    ├── run_rocm-6.3.yml
    └── regression_rocm-6.3.yml
```

## Stages and Entry Points

`common.yml` defines the stages `opalx_build`, `opalx_test` and `opalx_regression`. Regression
jobs depend only on the release build (via `needs:`), so they start as soon as the build
finishes, even if unit test stages are still running.

Each entry point (`cscs-openmp.yml`, `cscs-gh200.yml`, `cscs-mi300.yml`) includes the build, run
and regression files of its architecture. The CSCS admin console points at these files.

## Architecture-Specific Regression Files

Each `regression_<arch>.yml` defines:

- `BUILD_ARCH`: the architecture name on the results site. The site compares runs of the same
  name over time, so change it only together with the site (see "Names" in the results
  repository README).
- `CTEST_SITE`: diagnostic site name.
- `SRUN_ARGS`: srun arguments (uenv, view, repo, GPU flags) used by the launcher wrapper.
- `WRAPPER`: optional MPI wrapper binary path (CUDA only).
- The runner (`.baremetal-runner-*`) and `SLURM_CPUS_PER_TASK`.
- A single job (e.g. `opalx-regression-openmp-release`) that extends the common template and
  the runner, `needs: ["opalx-build-<arch>-release"]` and sets
  `BUILD_DIR: "build-$CI_COMMIT_SHORT_SHA-release"`.

## The Launcher Wrapper

**File:** `ci/cscs/regression-test-mpirun-wrapper`

`regtest.sh` starts every regression test with the command given by `--launcher`, here
`regression-test-mpirun-wrapper -np {ranks}`. `{ranks}` is 1, as in the `.local` scripts of
regression-tests-x. The wrapper turns `-np N` into

```bash
srun ${SRUN_ARGS} -n N ${WRAPPER} env MALLOC_CHECK_=0 GLIBC_TUNABLES=glibc.malloc.tcache_count=0 <opalx> <test>.in --info 2
```

so each test runs inside the job's allocation with the uenv. The unit tests use the srun
command the build configured into ctest (`MPIEXEC_EXECUTABLE`, `MPIEXEC_PREFLAGS`).

## Step-by-Step Job Execution

The shared template `.opalx-regression-common` in `regression-common.yml` runs these steps.
Every step prints a `=== <step> ===` marker in the job log.

1. **PR number extraction.** `.pr-number-extractor` from `common.yml` sets `CDASH_LABEL`.
2. **Build path.** `BUILD_PATH=$(pwd)/$BUILD_DIR`.
3. **Publish branch.** CSCS CI runs PRs as merge pipelines with refs like `__CSCSCI__pr508`.
   These become `pr-508`, `master` stays `master`, and any other branch gets an empty
   `OPALX_PUBLISH_BRANCH` (tested, not published).
4. **PR title.** For `pr-<N>` the title is fetched from the GitHub API (no token, the repo is
   public) and appended: `pr-527 Fix remaining OPALX and external dependency warnings`.
5. **Run branch name.** `RUN_BRANCH` is the publish branch (or the ref name) with every
   character outside `[A-Za-z0-9._-]` replaced by `_`, the rule NightlyBuildX used, so new
   runs join the runs imported from opal-live-doc:
   `pr-508 CI/CD and regression testing fixes at CSCS` becomes
   `pr-508_CI_CD_and_regression_testing_fixes_at_CSCS`. The site shows underscores as spaces.
6. **Build artifact.** Lists `$BUILD_PATH/bin/` and adds `$BUILD_PATH/lib` to
   `LD_LIBRARY_PATH`.
7. **ctest paths.** `CTestTestfile.cmake` and `DartConfiguration.tcl` contain the build job's
   working directory (e.g. `f7t/1/.../build-89e66c54-release`); `sed` rewrites it to
   `$BUILD_PATH`, as the unit test stages (`run_*.yml`) do.
8. **PSI Gitea key.** `PSI_GIT_SSHKEY` (base64 or PEM with whitespace removed) is written to
   `~/.ssh/psi_gitea_key` and checked with `ssh-keygen`. The job fails without it, because
   the tools come from the private results repository. The log shows the public key, its
   fingerprint and Gitea's greeting, which names the key and says whether it is a user or a
   deploy key.
9. **Tools.** A sparse, blob-filtered clone of `tools/` of `RESULTS_REPO` at
   `RESULTS_TOOLS_REF` (default `main`). If the key cannot read the repository, the job stops
   here. A `git push --dry-run` of a new branch then checks write access; without it the log
   warns that publishing will fail.
10. **Python and ctest.** The tools need Python 3.9 or newer, and ctest 3.21 or newer for its
    JUnit report; SLES 15's system python3 is 3.6. A micromamba environment with `python` and
    `cmake` from conda-forge is kept at `$SCRATCH/opalx-regtest-tools-$(uname -m)` and rebuilt
    when it is missing or Python cannot import `site`. Only `PATH` changes: the environment's
    binaries find their libraries by RPATH, so srun, git and ssh keep the system libraries.
11. **regression-tests-x.** Shallow clone of `master`.
12. **Tests.** One call:

    ```bash
    bash opalx-regression/tools/regtest.sh --arch "$BUILD_ARCH" --build "$BUILD_PATH" \
      --source "$CI_PROJECT_DIR" --tests regression-tests-x --branch "$RUN_BRANCH" \
      --launcher "$CI_PROJECT_DIR/ci/cscs/regression-test-mpirun-wrapper -np {ranks}" \
      --timeout 1200 --out "$RESULTS_DIR"
    ```

    It runs `ctest -L unit`, then every test directory with `<name>.in`, a reference `.stat`
    and no `disabled` file, each in a copy under `TMPDIR` (set to `$CI_PROJECT_DIR/regtest-tmp`
    on scratch). It applies the `.rt` checks to the new and the reference `.stat` files and
    writes `runs/<RUN_BRANCH>/<BUILD_ARCH>/<run_id>/` with `run.json` and `series.json.gz`.
    `--source` is needed because `CMakeCache.txt` names the build job's checkout. A test that
    runs longer than 1200 s (the NightlyBuildX limit) counts as crashed. The exit code (0
    passed, 3 tests failed, other values errors) is kept in `REGTEST_EXIT`.
13. **Publish.** For master and PRs, `publish_run.sh <run dir> $RESULTS_REPO` clones only that
    run directory, commits it and pushes to `main`, rebasing and retrying up to 5 times if
    another job pushed first. The site's `index.json` is written by the Pages build, so
    concurrent jobs cannot conflict. Feature branches are not published; their run is in the
    job artifacts.
14. **Job status.** The job exits with `REGTEST_EXIT`, so a run with failing tests is published
    and the job still fails.

Artifacts (`when: always`, one week): `regression-results-$BUILD_ARCH`, the packed run.

## Key Workarounds and Their Reasons

- **ctest path fixing.** Without it, test executables launched via srun cannot find
  `libippl.so` (step 7).
- **`--source`.** `regtest.sh` takes branch and commit from git; the build artifact's
  `CMakeCache.txt` points at a checkout that does not exist in this job.
- **Launcher wrapper.** Starts each test with srun and the uenv without changing
  regression-tests-x.
- **Captured exit code.** Failing runs are still published and the job is still failed.
- **Tools environment validation.** A cached conda environment can break when the uenv or
  the scratch file system changes; importing `site` catches that, not only a missing binary.

## External Repositories

| Repository | Branch | Purpose |
|-----------|--------|---------|
| `https://github.com/OPALX-project/regression-tests-x.git` | `master` | Test inputs, reference `.stat` files, `.rt` checks |
| `git@gitea.psi.ch:AMAS/opalx-regression.git` | `main` | `tools/` (`regtest.sh`, `run_suite.py`, `pack_run.py`, `publish_run.sh`) and the published runs under `data/runs/`; Pages site at `amas.pages.psi.ch/opalx-regression` |

## Debugging Tips

1. **Find the step** by its `=== <step> ===` marker.
2. **`libippl.so: cannot open shared object file`**: the ctest path fix failed. Check that
   `CACHE_DIR` was read from `CMakeCache.txt`.
3. **`the key in PSI_GIT_SSHKEY cannot read ...`**: the key has no access to
   AMAS/opalx-regression. Gitea's greeting just above names the key. Add the public key that
   step 8 printed to the repository as a deploy key with write access, or set a new key.
4. **`the key can read but not push`**: a deploy key without write access, or a user without
   write permission on the repository.
5. **`ctest wrote no unit.xml`**: ctest is older than 3.21, or not the one from the tools
   environment. The `ctest --version` line of step 10 shows which one ran.
6. **`Fatal Python error: init_import_site` / `unknown encoding`**: the tools environment is
   corrupted. Delete `$SCRATCH/opalx-regtest-tools-$(uname -m)`; the next job recreates it.
7. **A test shows as crashed** on the site: its log tail is on the run page. Exit code -9
   means it hit the 1200 s timeout.
8. **YAML parse errors in GitLab**: quote `echo` lines containing `: ` and use `>-` block
   scalars for multi-line commands.

## Remaining TODOs

- [ ] Switch regression jobs from all PRs to **master-only** once fully verified
- [ ] Increase `SLURM_TIMELIMIT` from 45 min if timeouts occur
- [ ] Review concurrency/resource usage of regression jobs if needed
- [ ] Check that the OPALX commit of a PR run (HEAD of the `__CSCSCI__` merge pipeline
      checkout) links to a commit that GitHub knows
- [ ] Delete the old NightlyBuildX environment `$SCRATCH/nightlybuildx-tools-*` once no
      pipeline uses it
