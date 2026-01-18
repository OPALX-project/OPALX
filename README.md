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
boost/1.82.0_slurm   
gsl/2.7                
gtest/1.13.0-1
```

## Building OPALX

```
git clone https://github.com/OPALX-project/OPALX.git
cd OPALX
./gen_OPALrevision
```

### Setting up cmake

#### cmake command for CPU build

Building OPALX without multi-threading (only MPI):
```
mkdir build_serial && cd build_serial
cmake .. \
    -DBUILD_TYPE=Debug \
    -DPLATFORMS=SERIAL
```

and for multi-threading with OpenMP:

```
mkdir build_openmp && cd build_openmp
cmake .. \
    -DBUILD_TYPE=Debug \
    -DPLATFORMS=OPENMP
```

In order to enable the compilation of unit tests, set `-DOPALX_ENABLE_UNIT_TESTS=ON` in the cmake command. The resulting executables will appear in `unit_tests` directory.


#### cmake command for GPU build
```
mkdir build_cuda && cd build_cuda
```

For example, for A100 with Amper80 Architecture (Gwendolen), and the debug mode, the cmake should be something like:

```
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
`-DIPPL_ENABLE_ALPINE` default `OFF`

`-DIPPL_ENABLE_TEST` default `OFF`

`-DOPALX_ENABLE_UNIT_TESTS` default `OFF`

`-DOPALX_ENABLE_EXAMPLES` default `OFF`

`-DOPALX_ENABLE_TESTS` default `OFF`

`-DOPALX_ENABLE_COVERAGE` default `OFF`

`-DOPALX_ENABLE_NSYS_PROFILER` default `OFF`

`-DOPALX_ENABLE_SANITIZER` default `OFF`

`-DOPALX_USE_ALTERNATIVE_VARIANT` default `OFF`

`-DOPALX_ENABLE_SCRIPTS` default `OFF`

### Compilation

Finally, compile OPALX with 
```
make
```
using single thread, and
```
make -j 4
```
using `4` threads for example.


## Job scripts
To execute opalx on merlin's gpus (compile for PASCAL61), the job script should looks like
```
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
```
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

The documentation has been moved to the [Wiki](https://gitlab.psi.ch/OPAL/src/wikis/home).
