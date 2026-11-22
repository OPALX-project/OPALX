#!/bin/bash -l
#SBATCH --job-name=bb-track12
#SBATCH --account=project_465002453
#SBATCH --partition=small-g
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
set -euo pipefail
module load LUMI/25.03 partition/G cpeAMD rocm buildtools
export OMP_NUM_THREADS=1 MPICH_GPU_SUPPORT_ENABLED=1 OMP_PROC_BIND=false
export MPICH_GPU_IPC_ENABLED=0
run=/flash/project_465002453/bb-validation-e24dc9218-20261120
test -f "$run/smokes-validated-locally"
cd "$run/track12"
exe=/flash/project_465002453/opalx-beambeam-recovery-20261120/build-hip/src/opalx
sha256sum "$exe" track12_timed.in input/*.emittedfromfile input/primary_fixed.fromfile
/usr/bin/time -p srun -n1 --gpus-per-task=1 "$run/select-gpu.sh" "$exe" track12_timed.in --info 4 > opalx.log 2>&1
test -s timing.dat
touch completed
