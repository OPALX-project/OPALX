#!/bin/bash
#SBATCH --job-name=adapt_bins_test_gpu
#SBATCH --error=output/bins_%j.err
#SBATCH --output=output/bins_%j.out
#SBATCH --time=00:20:00
#SBATCH --nodes=1                   # Request node
#SBATCH --ntasks-per-node=4        # ranks per node
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=8           # "threads" per task (for e.g. multithreading in Kokkod:parallel_for?)
#SBATCH --cluster=gmerlin6 # gmerlin6
#SBATCH --partition=gwendolen # Mandatory, as gwendolen is not the default partition
#SBATCH --account=gwendolen   # Mandatory, as gwendolen is not the default account
##SBATCH --exclusive
##SBATCH --nodelist=merlin-c-001   # Modify node list if needed for non-GPU nodes
#SBATCH --gpus=4

# for gpu: use "--gpus=1", "--cluster=gmerlin6" and "--partition=gpu-short" instead of "--cluster=merlin6", "--partition=hourly"


export OMP_PROC_BIND=spread # I guess?
export OMP_PLACES=threads
export OMP_NUM_THREADS=$(nproc) # set this to how many threads the program should be run with!

echo "Number of threads: $(nproc)"

#source /data/user/liemen_a/.bashrc # for load_ippl...
#load_ippl
#module clear && module use unstable && module load cmake/3.20.5 gcc/12.3.0 gtest/1.13.0-1 openmpi/4.1.5_slurm cuda/12.1.1
#module clear && module use unstable && module load gcc/12.3.0 gtest/1.13.0-1 openmpi/4.1.5_slurm && module use Libraries && module load ucx/1.14.1_slurm fftw/3.3.10_merlin6 boost gsl hdf5 H5hut cuda/12.1.1 cmake/3.25.2

cd /data/user/liemen_a/ippl/build_ippl_cuda/
# cmake ../ippl/ -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_STANDARD=20 -DIPPL_PLATFORMS="CUDA;OPENMP" -DKokkos_ARCH_AMPERE80=ON -DUSE_ALTERNATIVE_VARIANT=ON -DENABLE_TESTS=ON -DENABLE_SOLVERS=ON -DENABLE_FFT=ON -DONLY_BINNING=ON
make -j $(nproc)
echo "Finished compiling. Now running the program..."

#cd /data/user/liemen_a/build_ippl_cuda/test/binning/test/pic3d/
#srun ./Binning_pic3d 8 8 8 10000000 1 --info 10 

cd /data/user/liemen_a/ippl/build_ippl_cuda/test/binning/test/alpine/
#srun ./BinningLandauDamping 16 16 16 1000000 10 FFT 0.01 LeapFrog Flattop --overallocate 2.0 --info 10
srun ./BinningLandauDamping 32 32 32 1000000 5 FFT 0.01 LeapFrog Landau --overallocate 2.0 --info 10  # --kokkos-map-device-id-by=mpi_rank
#srun ./BinningLandauDamping 16 16 16 1000000 1 FFT 0.01 LeapFrog Landau --overallocate 2.0 --info 10



# srun ./Binning_pic3d 8 8 8 1000000 1 --info 10 --kokkos-disable-cuda # use this to run ONLY on CPU


# srun ./Binning_pic3d 8 8 8 1000000 1 --info 10 --kokkos-map-device-id-by=mpi_rank
# srun --cpus-per-task=$(nproc) ./Binning_pic3d 8 8 8 1000000 1 --info 10



