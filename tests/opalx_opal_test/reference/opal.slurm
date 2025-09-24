#!/bin/bash
#SBATCH --job-name=opal
#SBATCH --cluster=merlin7
#SBATCH --partition=hourly      # Using 'hourly' will grant higher priority
#SBATCH --nodes=1               # No. of nodes
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1      # No. of OMP threads
#SBATCH --time=01:00:00         # Define max time job will run (e.g. here 5 mins)
#SBATCH --hint=multithread    # Without hyperthreading
#SBATCH --exclusive         

#SBATCH --output=opal.out  # Name of output file
#SBATCH --error=opal.err  # Name of error file

module load gcc/10.3.0 openmpi/4.0.5 OPAL/2021.1.0

srun --mpi=pmix opal ref-100steps.in --info 10
