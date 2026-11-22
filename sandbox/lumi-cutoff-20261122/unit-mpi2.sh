#!/bin/bash -l
#SBATCH --job-name=bb-cutoff-mpi2
#SBATCH --account=project_465002453
#SBATCH --partition=dev-g
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --gpus-per-node=2
#SBATCH --cpus-per-task=4
#SBATCH --time=00:10:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
set -euo pipefail
run=/flash/project_465002453/bb-cutoff-track12-20261122
src=/flash/project_465002453/opalx-beambeam-recovery-20261120
cd "$run/unit-mpi2"
trap 'result=$?; printf "%s\n" "$result" > exit-status; exit "$result"' EXIT
test -f "$run/build-and-unit.completed"
module load LUMI/25.03 partition/G cpeAMD rocm buildtools
export OMP_NUM_THREADS=1 MPICH_GPU_SUPPORT_ENABLED=1 OMP_PROC_BIND=false
export MPICH_GPU_IPC_ENABLED=0
(cd "$src" && sha256sum -c "$run/source.sha256")
sha256sum -c "$run/build/binaries.sha256"
srun -n2 --gpus-per-task=1 "$run/select-gpu.sh" "$run/bin/device-info" > gpu-binding.log 2>&1
srun -n2 --gpus-per-task=1 "$run/select-gpu.sh" "$src/build-hip/unit_tests/SpaceCharge/TestBeamBeamFieldServices" > unit-n2.log 2>&1
touch completed
