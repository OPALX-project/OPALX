#!/bin/bash
#SBATCH --partition=hourly      # Using 'hourly' will grant higher priority
#SBATCH --nodes=1               # No. of nodes
#SBATCH --ntasks-per-node=1     # No. of MPI ranks per node. Merlin CPU nodes have 22 cores
#SBATCH --cpus-per-task=22      # No. of OMP threads
#SBATCH --time=00:30:00         # Define max time job will run (e.g. here 5 mins)
#SBATCH --hint=nomultithread    # Without hyperthreading
##SBATCH --exclusive            # The allocations will be exclusive if turned on (remove extra hashtag to turn on)

#SBATCH --output=multi.out  # Name of output file
#SBATCH --error=multi.err    # Name of error file

export OMP\_NUM\_THREADS=22
export OMP\_PROC\_BIND=spread
export OMP\_PLACES=threads


srun --cpus-per-task=22  /data/user/binder_j/opalx-elements/build/src/opalx multi.in --info 10
