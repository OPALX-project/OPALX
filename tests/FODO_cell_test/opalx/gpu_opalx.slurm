#!/bin/bash
#SBATCH --job-name=gpu_fodo_opalx
#SBATCH --cluster=gmerlin7
#SBATCH --partition=a100-hourly
#SBATCH --gpus=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-gpu=2 
#SBATCH --time=01:00:00        

#SBATCH --output=gpu_fodo_opalx.out  # Name of output file
#SBATCH --error=gpu_fodo_opalx.err    # Name of error file

module purge
module use Spack unstable
module load gcc/13.2.0 openmpi/5.0.7-dnpr-A100-gpu
module load boost/1.82.0-lgrt fftw/3.3.10.6-zv2b-omp gnutls/3.8.9-mcdr googletest/1.14.0-msmu gsl/2.7.1-hxwy h5hut/2.0.0rc7-zy7s openblas/0.3.29-zkwb cmake/3.31.6-oe7u

srun /data/user/binder_j/bin/opalx/build_gpu/src/opalx test.in --info 10
