#!/bin/bash
#SBATCH --job-name=adapt_bins_test_gpu
#SBATCH --error=output/bins_%j.err
#SBATCH --output=output/bins_%j.out
#SBATCH --time=00:02:00
#SBATCH --nodes=1                   # Request node
#SBATCH --ntasks-per-node=4         # cores per node
#SBATCH --mem-per-cpu=4G
#SBATCH --cpus-per-task=1           # "threads" per task (for e.g. multithreading in Kokkod:parallel_for?)
#SBATCH --cluster=gmerlin6 # gmerlin6
#SBATCH --partition=gwendolen # Mandatory, as gwendolen is not the default partition
#SBATCH --account=gwendolen   # Mandatory, as gwendolen is not the default account
##SBATCH --exclusive
##SBATCH --nodelist=merlin-c-001   # Modify node list if needed for non-GPU nodes
#SBATCH --gpus=1

# for gpu: use "--gpus=1", "--cluster=gmerlin6" and "--partition=gpu-short" instead of "--cluster=merlin6", "--partition=hourly"


export OMP_PROC_BIND=spread # I guess?
export OMP_PLACES=threads

cd /data/user/liemen_a/build_ippl/test/binning/
## ./Binning_pic3d 32 32 32 1000 10 --info 10
mpirun ./Binning_pic3d 32 32 32 100000 1 --info 10



