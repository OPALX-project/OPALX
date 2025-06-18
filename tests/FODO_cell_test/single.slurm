#!/bin/bash
#SBATCH --partition=hourly      # Using 'hourly' will grant higher priority
#SBATCH --nodes=1               # No. of nodes
#SBATCH --ntasks-per-node=1     # No. of MPI ranks per node. Merlin CPU nodes have 44 cores
#SBATCH --cpus-per-task=22      # No. of OMP threads
#SBATCH --time=00:45:00         # Define max time job will run (e.g. here 5 mins)
#SBATCH --hint=nomultithread    # Without hyperthreading
##SBATCH --exclusive            # The allocations will be exclusive if turned on (remove extra hashtag to turn on)

#SBATCH --output=single.out  # Name of output file
#SBATCH --error=single.err    # Name of error file

export OMP_NUM_THREADS=22
export OMP_PROC_BIND=spread
export OMP_PLACES=threads


srun --cpus-per-task=22  /data/user/binder_j/opalx-elements/build/src/opalx opalx/multi.in --info 10
