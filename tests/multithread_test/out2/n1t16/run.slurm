#!/bin/bash
#SBATCH --job-name=opalx-1-16
#SBATCH --cluster=merlin7
#SBATCH --partition=hourly      # Using 'hourly' will grant higher priority
#SBATCH --nodes=1              # No. of nodes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16      # No. of OMP threads
#SBATCH --time=01:00:00         # Define max time job will run (e.g. here 5 mins)
#SBATCH --hint=multithread    # Without hyperthreading
#SBATCH --exclusive         

#SBATCH --output=1-16.out  # Name of output file
#SBATCH --error=1-16.err    # Name of error file

export OMP_PROC_BIND=spread
export OMP_PLACES=threads

module purge
module use Spack unstable
module load gcc/13.2.0 openmpi/5.0.7-dnpr-A100-gpu
module load boost/1.82.0-lgrt fftw/3.3.10.6-zv2b-omp gnutls/3.8.9-mcdr googletest/1.14.0-msmu gsl/2.7.1-hxwy h5hut/2.0.0rc7-zy7s openblas/0.3.29-zkwb cmake/3.31.6-oe7u
srun --mpi=pmix /data/user/binder_j/bin/opalx/build_openmp/src/opalx template.in --info 10
