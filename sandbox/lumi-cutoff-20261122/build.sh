#!/bin/bash -l
#SBATCH --job-name=bb-cutoff-build
#SBATCH --account=project_465002453
#SBATCH --partition=dev-g
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --time=00:30:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
set -euo pipefail
run=/flash/project_465002453/bb-cutoff-track12-20261122
src=/flash/project_465002453/opalx-beambeam-recovery-20261120
cd "$run/build"
trap 'result=$?; printf "%s\n" "$result" > exit-status; exit "$result"' EXIT
module load LUMI/25.03 partition/G cpeAMD rocm buildtools
export OMP_NUM_THREADS=1 MPICH_GPU_SUPPORT_ENABLED=1 OMP_PROC_BIND=false
export MPICH_GPU_IPC_ENABLED=0
{
    hostname
    date -u
    module -t list
    git -C "$src" rev-parse HEAD
    git -C "$src" status --short
    sha256sum "$run/cutoff.patch" "$run/pre-cutoff/opalx"
    sed -n '/^CMAKE_BUILD_TYPE:/p; /^IPPL_GIT_TAG:/p; /^Kokkos_ARCH_AMD_GFX90A:/p; /^CMAKE_CXX_COMPILER:/p' "$src/build-hip/CMakeCache.txt"
} > provenance.txt 2>&1
(cd "$src" && sha256sum -c "$run/source.sha256")
cmake --build "$src/build-hip" --target opalx_exe TestBeamBeamFieldServices -j20 > build.log 2>&1
cp -p "$src/build-hip/src/opalx" "$run/bin/opalx"
amdclang++ -x hip --offload-arch=gfx90a "$run/device-info.cpp" -o "$run/bin/device-info"
sha256sum "$run/bin/opalx" "$src/build-hip/unit_tests/SpaceCharge/TestBeamBeamFieldServices" > binaries.sha256
srun -n1 --gpus-per-task=1 "$run/select-gpu.sh" "$run/bin/device-info" > gpu-binding.log 2>&1
srun -n1 --gpus-per-task=1 "$run/select-gpu.sh" "$src/build-hip/unit_tests/SpaceCharge/TestBeamBeamFieldServices" > unit-n1.log 2>&1
(cd "$src" && sha256sum -c "$run/source.sha256")
touch "$run/build-and-unit.completed"
