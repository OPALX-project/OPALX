# OPALX in Docker

A CPU (OpenMP) build of OPALX in a container, for trying OPALX without
installing compilers, MPI or any of its dependencies. All you need is Docker.

The image is built for x86-64 and ARM64, and Docker picks the one that matches
your machine. It runs natively on Linux, on Windows (WSL 2 or Docker Desktop),
on Intel and Apple Silicon Macs, and on ARM Linux machines. GPU builds are not
provided yet. [Step by step on a Mac](#step-by-step-on-a-mac) walks through
installing Docker, building the image and running the tests.

## Run a simulation

From the directory that holds your input file:

```bash
docker run --rm --user "$(id -u):$(id -g)" -v "$PWD:/data" \
    ghcr.io/opalx-project/opalx:cpu case.in --info 2
```

Output files are written next to the input file. Any other files the input
refers to (distributions, field maps) must be in the same directory or below
it, because only that directory is visible inside the container.

| Option | Effect |
|---|---|
| `--user "$(id -u):$(id -g)"` | output files belong to you instead of root |
| `-e OMP_NUM_THREADS=8` | OpenMP threads per MPI rank (default: all CPUs) |
| `-e OPALX_NP=4 --shm-size=8g` | 4 MPI ranks; MPI needs more shared memory than Docker's 64 MB default |

Particles that are emitted or generated randomly are drawn per thread, so a
different thread count gives a different (statistically equivalent) sample.

`--help` works as an argument, and `--version-full --info 5` shows the OPALX
version and git revision. Any command on the image's `PATH` runs instead of
OPALX, e.g. `bash` for a shell (add `-it`).

On Windows PowerShell, leave out `--user` and use `-v "${PWD}:/data"`.

### Tags

| Tag | Contents |
|---|---|
| `cpu`, `latest` | current `master` |
| `1.0.0-cpu` | release `OPALX-1.0.0` |
| `sha-<commit>-cpu` | a specific commit |

## Build the image yourself

From the repository root:

```bash
docker build -f docker/Dockerfile -t opalx:cpu --build-arg GIT_SHA=$(git rev-parse HEAD) .
```

The image is built for your machine's architecture; add `--platform linux/arm64`
or `linux/amd64` for the other one (this needs QEMU emulation and is much
slower). The build downloads IPPL, Kokkos, Heffte, HDF5, H5hut and LAPACK, so
it needs network access. Build arguments:

| Argument | Default | |
|---|---|---|
| `BUILD_TYPE` | `Release` | CMake build type |
| `JOBS` | number of CPUs | parallel compile jobs; lower it if the build runs out of memory |
| `GIT_SHA` | `unknown` | git revision reported by OPALX |
| `UBUNTU_VERSION` | `24.04` | base image |

To check an image with the smoke regression tests used in CI:

```bash
git clone --depth=1 https://github.com/OPALX-project/regression-tests-x.git ../regression-tests-x
docker/smoke-test.sh opalx:cpu ../regression-tests-x
```

The tests are cloned next to the repository rather than into it, so that they
are not sent to the next image build. Other tests can be named after the
checkout path, e.g.
`docker/smoke-test.sh opalx:cpu ../regression-tests-x RFCavity`. Some of them
currently miss their tolerances with any build, including the serial build
CI uses, so only the default CI subset is a reliable check of an image. ARM64
results differ from x86-64 ones in the last digits, because ARM64 compilers
fuse multiply-adds.

## Development environment

The `dev` stage has the toolchain only (GCC, gfortran, Open MPI, CMake). The
source stays on your machine and is mounted into the container:

```bash
docker build -f docker/Dockerfile --target dev -t opalx-dev .
docker/dev-shell.sh
```

Inside the container, the repository is at `/work`:

```bash
cmake -S /work -B /work/build_openmp -DCMAKE_BUILD_TYPE=Release -DPLATFORMS=OPENMP -DOPALX_ENABLE_UNIT_TESTS=ON
cmake --build /work/build_openmp -j "$(nproc)"
ctest --test-dir /work/build_openmp --output-on-failure
```

Always run CMake inside the container: the build directory records container
paths. Set `OPALX_WORK` to mount a different directory at `/work`, for
example a parent directory that also holds build and run directories.

## Step by step on a Mac

These steps set up Docker on a Mac with Apple Silicon or an Intel chip, build
the image and check it with the smoke tests. Run the commands in Terminal.

1. Install [Docker Desktop](https://docs.docker.com/desktop/setup/install/mac-install/),
   with the Apple Silicon or Intel installer to match your Mac, and start it.
   `docker version` then shows both a Client and a Server section.

2. Install git and Python, unless `git --version` and `python3 --version`
   already work:

   ```bash
   xcode-select --install
   ```

3. Get the OPALX source and build the image, which takes several minutes:

   ```bash
   git clone https://github.com/OPALX-project/OPALX.git
   cd OPALX
   docker build -f docker/Dockerfile -t opalx:cpu --build-arg GIT_SHA=$(git rev-parse HEAD) .
   ```

   If the build stops with `Killed` or exit code 137, it ran out of memory.
   Docker Desktop gets half of the Mac's memory by default. Raise the limit
   under *Settings → Resources → Advanced*, or add `--build-arg JOBS=4` to
   compile fewer files at a time.

   Instead of building the image, you can pull the published one with
   `docker pull ghcr.io/opalx-project/opalx:cpu` and use that name in place of
   `opalx:cpu` below. The source is still needed for the test scripts.

4. Check that the image matches your Mac. This prints `aarch64` on Apple
   Silicon and `x86_64` on Intel:

   ```bash
   docker run --rm --entrypoint uname opalx:cpu -m
   ```

   An `x86_64` image on Apple Silicon runs under emulation and is much
   slower. That happens, for example, when `DOCKER_DEFAULT_PLATFORM` is set
   to `linux/amd64`.

5. Get the regression tests and run the smoke tests. The last line of the
   output should read `All 3 tests passed`:

   ```bash
   git clone --depth=1 https://github.com/OPALX-project/regression-tests-x.git ../regression-tests-x
   docker/smoke-test.sh opalx:cpu ../regression-tests-x
   ```

6. To run the unit tests as well, build the dev image and compile OPALX in it
   as described in [Development environment](#development-environment).
