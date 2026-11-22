#!/bin/bash -l
#SBATCH --job-name=bb-exp4-cutoff
#SBATCH --account=project_465002453
#SBATCH --partition=small-g
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
set -euo pipefail
campaign=/flash/project_465002453/bb-experiment4-cutoff-20261122
validated=/flash/project_465002453/bb-cutoff-track12-20261122
case_name=${1:?Specify cutoff or old}
case "$case_name" in
    cutoff)
        binary="$validated/bin/opalx"
        expected=f44bf28796a8134cc96d8a3afbb5e8a62a90adf1be72c152dd28123266da73c3 ;;
    old)
        binary="$validated/pre-cutoff/opalx"
        expected=31e97e63d5208ef210c4188b717b48ccd739fba37874f38735e08754372d1307 ;;
    *) exit 2 ;;
esac
cd "$campaign/$case_name"
trap 'result=$?; printf "%s\n" "$result" > exit-status; exit "$result"' EXIT
test ! -e opalx.log
test -f "$validated/unit-mpi2/completed"
module load LUMI/25.03 partition/G cpeAMD rocm buildtools
export OMP_NUM_THREADS=1 MPICH_GPU_SUPPORT_ENABLED=1 OMP_PROC_BIND=false
export MPICH_GPU_IPC_ENABLED=0
printf '%s  %s\n' "$expected" "$binary" | sha256sum -c -
sha256sum -c inputs.sha256
{
    date -u
    hostname
    scontrol show job "$SLURM_JOB_ID"
    sha256sum "$binary" *.in input/*
} > execution-fingerprint.txt
srun -n1 --gpus-per-task=1 "$validated/select-gpu.sh" "$validated/bin/device-info" > gpu-binding.log 2>&1
/usr/bin/time -p srun -n1 --gpus-per-task=1 "$validated/select-gpu.sh" "$binary" gamma_gamma_pairs-20mm-dynamic-xy-staged.in --info 4 > opalx.log 2>&1
test -s timing.dat
touch completed
