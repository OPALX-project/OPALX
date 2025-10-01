#!/bin/bash
#SBATCH --job-name=opalx-1-128
#SBATCH --cluster=merlin7
#SBATCH --partition=hourly      # Using 'hourly' will grant higher priority
#SBATCH --nodes=1               # No. of nodes
#SBATCH --ntasks=128
#SBATCH --cpus-per-task=1      # No. of OMP threads
#SBATCH --time=01:00:00         # Define max time job will run (e.g. here 5 mins)
#SBATCH --hint=multithread    # Without hyperthreading
#SBATCH --exclusive         

#SBATCH --output=1-128.out  # Name of output file
#SBATCH --error=1-128.err    # Name of error file

module purge
module use unstable
module load cmake/3.25.2
module load gcc/10.3.0       
module load openmpi/4.0.5
module load fftw/3.3.10
module load hdf5/1.10.7
module load H5hut/2.0.0rc6
module load boost/1.76.0
module load gtest/1.11.0
module load gnutls/3.7.10
module load gsl/2.7
module load cuda/12.9.1
module load Python/3.9.10
module load pmix/5.0.8
srun --mpi=pmix /data/user/binder_j/bin/opalx-elements/build2/src/opalx template.in --info 10
