# OPALX

## Dependencies

```
cmake/3.25.2
openmpi/4.1.5_slurm
fftw/3.3.10_merlin6    
gsl/2.7                
H5hut/2.0.0rc6_slurm
gcc/12.3.0             
boost/1.82.0_slurm     
gtest/1.13.0-1         
hdf5/1.10.8_slurm     
gnutls/3.5.19
cuda/12.8.1
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
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_CXX_STANDARD=20 \
    -DIPPL_ENABLE_FFT=ON \
    -DIPPL_ENABLE_SOLVERS=ON \
    -DIPPL_ENABLE_ALPINE=OFF \
    -DIPPL_ENABLE_TESTS=OFF  \
    -DIPPL_PLATFORMS=serial
```

and for multi-threading with OpenMP:

```
mkdir build_openmp && cd build_openmp
cmake .. \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_CXX_STANDARD=20 \
    -DIPPL_ENABLE_FFT=ON \
    -DIPPL_ENABLE_SOLVERS=ON \
    -DIPPL_ENABLE_ALPINE=OFF \
    -DIPPL_ENABLE_TESTS=OFF  \
    -DIPPL_PLATFORMS=openmp
```

#### cmake command for GPU build
```
mkdir build_cuda && cd build_cuda
```

For example, for A100 with Amper80 Architecture (Gwendolen), and the debug mode, the cmake should be something like:

```
cmake .. \
    -DCMAKE_BUILD_TYPE=Debug \
    -DIPPL_PLATFORMS=CUDA \
    -DKokkos_ARCH_AMPERE80=ON \
    -DCMAKE_CXX_STANDARD=20 \
    -DIPPL_ENABLE_FFT=ON \
    -DIPPL_ENABLE_SOLVERS=ON \
    -DIPPL_ENABLE_ALPINE=OFF \
    -DIPPL_ENABLE_TESTS=OFF
```

For the release mode, use `Release` instead of `Debug` as the argument for `-DCMAKE_BUILD_TYPE`. For other GPUs use the correct flag for their corresponding architecture. For example, for P100 or GTX 1080 with Pascal61 architecture on Merlin login node, use `-DKokkos_ARCH_PASCAL61=ON` instead of `-DKokkos_ARCH_AMPERE80=ON`. 

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
