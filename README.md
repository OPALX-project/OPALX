# OPALX

# BUILDING OPALX on Merlin

## Modules needed OPENMP build

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

## Clone repo and build opalx with OPENMP 

```
$ git clone https://github.com/OPALX-project/OPALX.git
$ cd OPALX
$ ./gen_OPALrevision
```

### To compile for OPENMP:
```
$ mkdir build_openmp && cd build_openmp
$ cmake .. \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_CXX_STANDARD=20 \
    -DIPPL_ENABLE_FFT=ON \
    -DIPPL_ENABLE_SOLVERS=ON \
    -DIPPL_ENABLE_ALPINE=OFF \
    -DIPPL_ENABLE_TESTS=OFF  \
    -DIPPL_PLATFORMS=openmp
```

### To compile for GPU, for example Amper80 on Gwendolen
```
$ mkdir build_cuda && cd build_cuda
```

in debug mode:

```
$ cmake .. \
    -DCMAKE_BUILD_TYPE=Debug \
    -DIPPL_PLATFORMS=CUDA \
    -DKokkos_ARCH_AMPERE80=ON \
    -DCMAKE_CXX_STANDARD=20 \
    -DIPPL_ENABLE_FFT=ON \
    -DIPPL_ENABLE_SOLVERS=ON \
    -DIPPL_ENABLE_ALPINE=OFF \
    -DIPPL_ENABLE_TESTS=OFF
```

and release (optimized) mode:
```
$ cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DIPPL_PLATFORMS=CUDA \
    -DKokkos_ARCH_AMPERE80=ON \
    -DCMAKE_CXX_STANDARD=20 \
    -DIPPL_ENABLE_FFT=ON \
    -DIPPL_ENABLE_SOLVERS=ON \
    -DIPPL_ENABLE_ALPINE=OFF \
    -DIPPL_ENABLE_TESTS=OFF
```

### To compile for other GPU architecture, like Pascal on the Merlin's login node
```
$ mkdir build_cuda_login && cd build_cuda_login
```

in debug mode:

```
$ cmake .. \
    -DCMAKE_BUILD_TYPE=Debug \
    -DIPPL_PLATFORMS=CUDA \
    -DKokkos_ARCH_PASCAL61=ON \
    -DCMAKE_CXX_STANDARD=20 \
    -DIPPL_ENABLE_FFT=ON \
    -DIPPL_ENABLE_SOLVERS=ON \
    -DIPPL_ENABLE_ALPINE=OFF \
    -DIPPL_ENABLE_TESTS=OFF
```

and release (optimized) mode:
```
$ cmake .. \
    -DCMAKE_BUILD_TYPE=Release \
    -DIPPL_PLATFORMS=CUDA \
    -DKokkos_ARCH_PASCAL61=ON \
    -DCMAKE_CXX_STANDARD=20 \
    -DIPPL_ENABLE_FFT=ON \
    -DIPPL_ENABLE_SOLVERS=ON \
    -DIPPL_ENABLE_ALPINE=OFF \
    -DIPPL_ENABLE_TESTS=OFF
```

### Submitting jobs on Gwendolen and Merlin GPUs
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
