#!/bin/bash
#SBATCH --partition=hourly      # Using 'hourly' will grant higher priority
#SBATCH --nodes=1               # No. of nodes
#SBATCH --ntasks-per-node=1     # No. of MPI ranks per node. Merlin CPU nodes have 44 cores
#SBATCH --cpus-per-task=16      # No. of OMP threads
#SBATCH --time=01:00:00         # Define max time job will run (e.g. here 5 mins)
#SBATCH --hint=nomultithread    # Without hyperthreading
##SBATCH --exclusive            # The allocations will be exclusive if turned on (remove extra hashtag to turn on)

#SBATCH --output=multi_run1-1e7-2.in.out  # Name of output file
#SBATCH --error=multi_run1-1e7-2.in.err    # Name of error file

export OMP_NUM_THREADS=16
export OMP_PROC_BIND=spread
export OMP_PLACES=threads


srun --cpus-per-task=16 /data/user/binder_j/opalx-elements/build_rel/src/opalx multi_run1-1e7-2.in --info 10
