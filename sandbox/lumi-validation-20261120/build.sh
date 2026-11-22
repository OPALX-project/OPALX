#!/bin/bash -l
#SBATCH --job-name=bb-build-smoke
#SBATCH --account=project_465002453
#SBATCH --partition=dev-g
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --time=01:00:00
set -euo pipefail
cd /flash/project_465002453/opalx-beambeam-recovery-20261120
module load LUMI/25.03 partition/G cpeAMD rocm buildtools
export OMP_NUM_THREADS=1
export MPICH_GPU_SUPPORT_ENABLED=1
export OMP_PROC_BIND=false
link_flags=$(sed -n 's/^CMAKE_EXE_LINKER_FLAGS:STRING=//p' /flash/project_465002453/opalx/build/CMakeCache.txt)
cmake -S . -B build-hip -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER=amdclang++ -DCMAKE_C_COMPILER=cc -DCMAKE_Fortran_COMPILER=ftn \
  -DOPALX_ENABLE_UNIT_TESTS=ON -DPLATFORMS=HIP -DIPPL_PLATFORMS=HIP \
  -DARCH=AMD_GFX90A -DKokkos_ARCH_AMD_GFX90A=ON -DKokkos_ARCH_ZEN3=ON \
  -DKokkos_ENABLE_HIP=ON -DCMAKE_HIP_ARCHITECTURES=gfx90a \
  -DCMAKE_HIP_FLAGS=--offload-arch=gfx90a -DHeffte_ENABLE_ROCM=ON \
  -DHeffte_ENABLE_GPU_AWARE_MPI=OFF -DBUILD_SHARED_LIBS=OFF \
  -DCMAKE_EXE_LINKER_FLAGS="$link_flags"
cmake --build build-hip -j20
git rev-parse HEAD
sha256sum build-hip/src/opalx
srun -n1 build-hip/unit_tests/Distribution/TestEmittedFromFile
srun -n1 build-hip/unit_tests/SpaceCharge/TestBeamBeamFieldServices
srun -n1 build-hip/unit_tests/Structure/TestBeamBeamDiagnosticsWriter
touch /flash/project_465002453/bb-validation-e24dc9218-20261120/build-and-unit.completed
