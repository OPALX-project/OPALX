# OPALX

## Manual and Documentation

[OPALX documentation](https://opalx-project.github.io/Manual/)

## Dependencies

In order to compile OPALX, make sure you have the following dependencies installed on your system:

```
cmake/3.25.2
openmpi/4.1.5_slurm
gcc/12.3.0             
gnutls/3.5.19
cuda/12.8.1
```

Other dependencies are fetched and installed in the opalx installation.

```
ippl/3.2.0
hdf5/1.10.8_slurm  
H5hut/2.0.0rc6_slurm  
gtest/1.13.0-1
```

## Building OPALX

```bash
git clone https://github.com/OPALX-project/OPALX.git
cd OPALX/tools
./gen_OPALXrevision # Note that this step is (now) optional
```

### Setting up cmake

#### cmake command for CPU build

Building OPALX without multi-threading (only MPI):
```bash
mkdir build_serial && cd build_serial
cmake .. \
    -DBUILD_TYPE=Debug \
    -DPLATFORMS=SERIAL
```

and for multi-threading with OpenMP:

```bash
mkdir build_openmp && cd build_openmp
cmake .. \
    -DBUILD_TYPE=Debug \
    -DPLATFORMS=OPENMP
```

In order to enable the compilation of unit tests, set `-DOPALX_ENABLE_UNIT_TESTS=ON` in the cmake command. The resulting executables will appear in `unit_tests` directory.


#### cmake command for GPU build
```bash
mkdir build_cuda && cd build_cuda
```

For example, for A100 with Amper80 Architecture (Gwendolen), and the debug mode, the cmake should be something like:

```bash
cmake .. \
    -DBUILD_TYPE=Debug \
    -DPLATFORMS=CUDA \
    -DARCH=AMPERE80
```

For the release mode, use `Release` instead of `Debug` as the argument for `-DBUILD_TYPE`. For other GPUs use the correct flag for their corresponding architecture. For example, for P100 or GTX 1080 with Pascal61 architecture on Merlin login node, use `-DARCH=PASCAL61` instead of `-DARCH=AMPERE80`. 

#### Notes:

- Use -DBUILD_TYPE=Release for optimized builds.
- ARCH is required for CUDA builds so OPALX can configure Kokkos properly.
- All IPPL/Kokkos flags (FFT, solvers, tests, ALPINE, `Kokkos_ARCH_*`, etc.) are now set automatically.

#### Further Options
| Flag | Default | Description |
|------|---------|-------------|
| `OPALX_EMBED_BUILD_METADATA` | OFF | Embeds user, machine, and date metadata into `BuildInfo.h`; leave OFF to reduce rebuild churn. |
| `OPALX_USE_INSTALLED_HDF5` | OFF | Disables the use of system-built HDF5 dependency. |
| `OPALX_USE_INSTALLED_H5HUT` | OFF |  Disables the use of system-built H5hut dependency. |
| `OPALX_USE_INSTALLED_GTEST` | OFF |  Disables the use of system-built GoogleTest dependency. |
| `IPPL_ENABLE_ALPINE` | OFF | Disables Alpine features in IPPL by default; set IPPL_ENABLE_ALPINE to ON to enable. |
| `IPPL_ENABLE_TEST` | OFF | Disables IPPL tests; corresponds to IPPL_ENABLE_TESTS OFF in IPPL default features. |
| `OPALX_ENABLE_UNIT_TESTS` | OFF | Disables building unit tests using GoogleTest. |
| `OPALX_ENABLE_EXAMPLES` | OFF | Disables building the Example module. |
| `OPALX_ENABLE_TESTS` | OFF | Disables building integration tests in the test/ directory. |
| `OPALX_ENABLE_COVERAGE` | OFF | Disables code coverage instrumentation. |
| `OPALX_ENABLE_HIP_PROFILER` | OFF | Enables the HIP Systems Profiler on HIP builds and adds the corresponding compile definition. |
| `OPALX_ENABLE_NSYS_PROFILER` | OFF | Disables Nvidia Nsight Systems Profiler; requires CUDA platform and adds compile definition -DOPALX_ENABLE_NSYS_PROFILER. |
| `OPALX_ENABLE_SANITIZER` | OFF | Disables sanitizer tools (e.g., AddressSanitizer). |
| `OPALX_USE_ALTERNATIVE_VARIANT` | OFF | Disables modified variant implementation; required for CUDA 12.2 + GCC 12.3.0 compatibility. |
| `OPALX_USE_STANDARD_FOLDERS` | OFF | Places generated binaries and libraries under the build tree `bin/` and `lib/` folders. |
| `OPALX_SKIP_FAILING_TESTS` | OFF | Skips building and running tests that are currently marked as failing. |
| `OPALX_ENABLE_SCRIPTS` | OFF | Disables generation of job script templates for benchmarks/tests. |
| `OPALX_FIELD_DEBUG` | OFF | Disables FieldSolver field-dump debug code; emits field dumps during simulation when enabled. |
| `OPALX_USE_KOKKOS_MATH_CONSTANTS` | ON | Sources `Physics.h` mathematical constants such as `pi`, `e`, and `log10e` from `Kokkos::numbers`; turn OFF to use literal fallback values. |

Enable flags with `-D<FLAG>=ON` during CMake configuration.

### Compilation

Finally, compile OPALX with 
```bash
make
```
using single thread, and
```bash
make -j 4
```
using `4` threads for example.


## Continuous Integration

[![CPU Serial](https://github.com/OPALX-project/OPALX/actions/workflows/cpu-serial.yml/badge.svg)](https://github.com/OPALX-project/OPALX/actions/workflows/cpu-serial.yml)
[![CPU OpenMP](https://github.com/OPALX-project/OPALX/actions/workflows/cpu-openmp.yml/badge.svg)](https://github.com/OPALX-project/OPALX/actions/workflows/cpu-openmp.yml)
[![GPU CUDA](https://github.com/OPALX-project/OPALX/actions/workflows/gpu-cuda.yml/badge.svg)](https://github.com/OPALX-project/OPALX/actions/workflows/gpu-cuda.yml)
[![GPU HIP](https://github.com/OPALX-project/OPALX/actions/workflows/gpu-hip.yml/badge.svg)](https://github.com/OPALX-project/OPALX/actions/workflows/gpu-hip.yml)
[![GPU SYCL](https://github.com/OPALX-project/OPALX/actions/workflows/gpu-sycl.yml/badge.svg)](https://github.com/OPALX-project/OPALX/actions/workflows/gpu-sycl.yml)

OPALX uses GitHub Actions for required and non-required pull request checks. The compile CI workflows run on non-draft pull requests to `master` when the `compile-ci` label is selected. Required checks have to pass before a PR can be merged; non-required checks provide additional feedback but are not part of the required branch protection status checks.

### Formatting Checks

- `Check clang-format`: checks the formatting of changed source and unit test files in the pull request using clang-format 21.1.8. This check is **not required**.

See [IPPL's formatting documentation](https://github.com/IPPL-framework/ippl/blob/master/doc/extras/CodeFormattingSetup.md) for instructions on setting up pre-commit hooks.

### Compilation

- `CPU Serial`: compiles debug OPALX and the unit tests with `PLATFORMS=SERIAL`. This check is **required**.
- `CPU OpenMP`: compiles debug OPALX and the unit tests with `PLATFORMS=OPENMP`. This check is **required**.
- `GPU CUDA`: compiles debug OPALX and the unit tests with CUDA for NVIDIA `AMPERE80`. This check is **required**.
- `GPU HIP`: compiles debug OPALX and the unit tests with HIP for AMD `AMD_GFX90A`. This check is **required**.
- `GPU SYCL`: compiles release OPALX and the unit tests with SYCL for Intel `INTEL_PVC`. This check is **not required**, because the CI runner can compile this configuration but cannot test it.

### Running Tests

- Unit tests run with the `CPU Serial` build.
- Unit tests run with the `CPU OpenMP` build.
- A subset of fast regression tests from [OPALX-project/regression-tests-x](https://github.com/OPALX-project/regression-tests-x) runs with the `CPU Serial` binary and checks results against the tests' `.rt` metrics. The current regression subset is:
  - `Drift-1-fromfile`
  - `Drift-4-multi-emit-open`
  - `FodoCell-fromfile`

## Job scripts
To execute opalx on merlin's gpus (compile for PASCAL61), the job script should looks like
```bash
#!/bin/bash
#SBATCH --error=merlin.error
#SBATCH --output=merlin.out
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cluster=gmerlin6
#SBATCH --partition=gpu-short
#SBATCH --account=merlin
##SBATCH --exclusive
#SBATCH --gpus=2
#SBATCH --nodelist=merlin-g-001

##unlink core
ulimit -c unlimited

srun ./opalx DriftTest-1.in  --info 10 --kokkos-map-device-id-by=mpi_rank
```

and for Gwendolen (compile for AMPERE80)
```bash
#!/bin/bash
#SBATCH --error=gwendolen.error
#SBATCH --output=gwendolen.out
#SBATCH --time=00:02:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --clusters=gmerlin6
#SBATCH --partition=gwendolen # Mandatory, as gwendolen is not the default partition
#SBATCH --account=gwendolen   # Mandatory, as gwendolen is not the default account
##SBATCH --exclusive
#SBATCH --gpus=2

##unlink core
ulimit -c unlimited

srun ./opalx DriftTest-1.in  --info 10 --kokkos-map-device-id-by=mpi_rank
```

The documentation has been moved to the [OPAL Landing Page](https://amas.pages.psi.ch/opal/).
