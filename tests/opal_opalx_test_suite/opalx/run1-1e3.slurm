#!/bin/bash
#SBATCH --partition=hourly      # Using 'hourly' will grant higher priority
#SBATCH --nodes=1               # No. of nodes
#SBATCH --ntasks-per-node=1     # No. of MPI ranks per node. Merlin CPU nodes have 44 cores
#SBATCH --cpus-per-task=1      # No. of OMP threads
#SBATCH --time=00:30:00         # Define max time job will run (e.g. here 5 mins)
#SBATCH --hint=nomultithread    # Without hyperthreading
##SBATCH --exclusive            # The allocations will be exclusive if turned on (remove extra hashtag to turn on)

#SBATCH --output=run1-1e3.in.out  # Name of output file
#SBATCH --error=run1-1e3.in.err    # Name of error file

export OMP\_NUM\_THREADS=1
export OMP\_PROC\_BIND=spread
export OMP\_PLACES=threads


srun --cpus-per-task=1 /data/user/binder_j/opalx-test/build/src/opalx run1-1e3.in --info 10
