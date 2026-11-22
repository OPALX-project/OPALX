#!/bin/bash -l
#SBATCH --job-name=bb-cutoff-track12
#SBATCH --account=project_465002453
#SBATCH --partition=small-g
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
set -euo pipefail
run=/flash/project_465002453/bb-cutoff-track12-20261122
src=/flash/project_465002453/opalx-beambeam-recovery-20261120
cd "$run/track12"
trap 'result=$?; printf "%s\n" "$result" > exit-status; exit "$result"' EXIT
test -f "$run/build-and-unit.completed"
test -f "$run/unit-mpi2/completed"
module load LUMI/25.03 partition/G cpeAMD rocm buildtools
export OMP_NUM_THREADS=1 MPICH_GPU_SUPPORT_ENABLED=1 OMP_PROC_BIND=false
export MPICH_GPU_IPC_ENABLED=0
(cd "$src" && sha256sum -c "$run/source.sha256")
sha256sum -c "$run/build/binaries.sha256"
sha256sum -c "$run/inputs.sha256"
{
    date -u
    hostname
    scontrol show job "$SLURM_JOB_ID"
    sha256sum "$run/bin/opalx" track12_timed.in input/*
} > execution-fingerprint.txt
srun -n1 --gpus-per-task=1 "$run/select-gpu.sh" "$run/bin/device-info" > gpu-binding.log 2>&1
/usr/bin/time -p srun -n1 --gpus-per-task=1 "$run/select-gpu.sh" "$run/bin/opalx" track12_timed.in --info 4 > opalx.log 2>&1
test -s timing.dat
touch completed
