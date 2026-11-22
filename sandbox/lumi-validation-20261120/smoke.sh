#!/bin/bash -l
#SBATCH --job-name=bb-smokes
#SBATCH --account=project_465002453
#SBATCH --partition=dev-g
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --gpus-per-node=2
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
set -euo pipefail
module load LUMI/25.03 partition/G cpeAMD rocm buildtools
export OMP_NUM_THREADS=1 MPICH_GPU_SUPPORT_ENABLED=1 OMP_PROC_BIND=false
# Cray MPICH's intra-node GPU IPC import fails with per-task ROCr visibility.
# Retain GPU-aware MPI but use its non-IPC transport for this validation.
export MPICH_GPU_IPC_ENABLED=0
run=/flash/project_465002453/bb-validation-e24dc9218-20261120
src=/flash/project_465002453/opalx-beambeam-recovery-20261120
exe=$src/build-hip/src/opalx
test -f "$run/build-and-unit.completed"
cd "$run"
{
  git -C "$src" rev-parse HEAD
  sha256sum "$exe"
  scontrol show job "$SLURM_JOB_ID"
} > execution-fingerprint.txt
amdclang++ -x hip --offload-arch=gfx90a "$run/device-info.cpp" -o "$run/device-info"
srun -n2 --gpus-per-task=1 "$run/select-gpu.sh" "$run/device-info"
for test in Distribution/TestEmittedFromFile PartBunch/TestBeamBeam SpaceCharge/TestBeamBeamFieldServices Structure/TestBeamBeamDiagnosticsWriter; do
  srun -n2 --gpus-per-task=1 "$run/select-gpu.sh" "$src/build-hip/unit_tests/$test"
done
for sigma in 3 2 1 0; do
  cd "$run/experiment-1/cases/${sigma}sigma"
  /usr/bin/time -p srun -n1 --gpus-per-task=1 "$run/select-gpu.sh" "$exe" "rigid_fields_${sigma}sigma.in" --info 4 > opalx.log 2>&1
done
cd "$run/experiment-2"
/usr/bin/time -p srun -n1 --gpus-per-task=1 "$run/select-gpu.sh" "$exe" track12_timed.in --info 4 > opalx.log 2>&1
cd "$run/experiment-3"
/usr/bin/time -p srun -n1 --gpus-per-task=1 "$run/select-gpu.sh" "$exe" gamma_gamma_pairs_timed.in --info 4 > opalx.log 2>&1
cd "$run/experiment-4"
/usr/bin/time -p srun -n1 --gpus-per-task=1 "$run/select-gpu.sh" "$exe" gamma_gamma_pairs-20mm-dynamic-xy-smoke.in --info 4 > opalx.log 2>&1
touch "$run/smokes-completed"
