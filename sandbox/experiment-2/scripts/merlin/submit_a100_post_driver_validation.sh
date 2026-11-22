#!/usr/bin/env bash
# Build-validation chain for Merlin6 after a CUDA driver update:
# unit tests (1 A100) -> Track-12 smoke (1 A100) -> production (4 A100s).

set -euo pipefail

stage="${1:-submit}"
if [[ $# -gt 0 ]]; then
    shift
fi

run_root=""
build_dir=""
cluster="gmerlin6"
partition="gwendolen"
account="gwendolen"
omp_threads=4

while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root) run_root="$2"; shift 2 ;;
        --build-dir) build_dir="$2"; shift 2 ;;
        --cluster) cluster="$2"; shift 2 ;;
        --partition) partition="$2"; shift 2 ;;
        --account) account="$2"; shift 2 ;;
        --omp-threads) omp_threads="$2"; shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 2 ;;
    esac
done

if [[ -z "$run_root" || -z "$build_dir" ]]; then
    echo "Usage: $0 [submit|unit-tests|smoke|production] --run-root PATH --build-dir PATH" >&2
    exit 2
fi

opalx="$build_dir/src/opalx"
unit_tests="$run_root/unit-tests"
common=(
    --run-root "$run_root"
    --build-dir "$build_dir"
    --cluster "$cluster"
    --partition "$partition"
    --account "$account"
    --omp-threads "$omp_threads"
)

if [[ "$stage" == "submit" ]]; then
    production="$run_root/production"
    smoke="$run_root/smoke"
    slurm="$run_root/slurm"
    test -x "$opalx"
    test -f "$production/track12_timed.in"
    test -d "$production/input"
    mkdir -p "$smoke" "$slurm" "$unit_tests"
    ln -sfn "$production/input" "$smoke/input"
    sed -E 's/REAL n_steps = [0-9]+;/REAL n_steps = 3;/' \
        "$production/track12_timed.in" > "$smoke/track12_timed.in"

    script_path="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
    submission_id="$(date -u +%Y%m%dT%H%M%SZ)_$$"
    submitted_script="$slurm/submit_post_driver_validation_${submission_id}.sh"
    cp "$script_path" "$submitted_script"
    chmod u+x "$submitted_script"

    unit_raw="$(sbatch --parsable --clusters="$cluster" --account="$account" \
        --job-name=opalx-unit-a100 --partition="$partition" --nodes=1 --ntasks=4 \
        --cpus-per-task=10 --gpus-per-node=1 --time=02:00:00 \
        --output="$slurm/unit_tests_%j.out" \
        "$submitted_script" unit-tests "${common[@]}")"
    unit_job="${unit_raw%%;*}"

    smoke_raw="$(sbatch --parsable --clusters="$cluster" --account="$account" \
        --dependency="afterok:$unit_job" --job-name=track12-smoke --partition="$partition" \
        --nodes=1 --ntasks=1 --cpus-per-task="$omp_threads" --gpus-per-node=1 \
        --time=00:15:00 --output="$slurm/smoke_%j.out" \
        "$submitted_script" smoke "${common[@]}")"
    smoke_job="${smoke_raw%%;*}"

    production_raw="$(sbatch --parsable --clusters="$cluster" --account="$account" \
        --dependency="afterok:$smoke_job" --job-name=track12-production --partition="$partition" \
        --nodes=1 --ntasks=4 --cpus-per-task="$omp_threads" --gpus-per-node=4 \
        --time=01:00:00 --output="$slurm/production_%j.out" \
        "$submitted_script" production "${common[@]}")"
    production_job="${production_raw%%;*}"

    printf 'Jobs: unit-tests=%s smoke=%s production=%s\n' \
        "$unit_job" "$smoke_job" "$production_job"
    printf 'Submitted script: %s\n' "$submitted_script"
    exit 0
fi

export OMP_NUM_THREADS="$omp_threads"
# Every stage runs within one DGX node.  Keep UCX on shared-memory and CUDA-IPC
# transports so it does not select the external mlx5 fabric for intra-node MPI.
export UCX_TLS="${UCX_TLS:-sm,self,cuda_copy,cuda_ipc}"

if [[ "$stage" == "unit-tests" ]]; then
    cd "$unit_tests"
    {
        printf 'utc_start=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        printf 'host=%s\n' "$(hostname)"
        printf 'slurm_job_id=%s\n' "${SLURM_JOB_ID:-}"
        printf 'build_dir=%s\n' "$build_dir"
        nvidia-smi --query-gpu=index,name,driver_version,uuid --format=csv,noheader
    } > run_manifest.txt
    ctest --test-dir "$build_dir" --output-on-failure -j1 2>&1 | tee ctest.log
    date -u +%Y-%m-%dT%H:%M:%SZ > completed
    exit 0
fi

case "$stage" in
    smoke) case_dir="$run_root/smoke"; mpi_ranks=1 ;;
    production) case_dir="$run_root/production"; mpi_ranks=4 ;;
    *) echo "Unknown stage: $stage" >&2; exit 2 ;;
esac

cd "$case_dir"
{
    printf 'case=%s\n' "$stage"
    printf 'utc_start=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'host=%s\n' "$(hostname)"
    printf 'slurm_job_id=%s\n' "${SLURM_JOB_ID:-}"
    printf 'mpi_ranks=%s\n' "$mpi_ranks"
    printf 'omp_threads=%s\n' "$omp_threads"
    printf 'opalx=%s\n' "$opalx"
    sha256sum "$opalx" track12_timed.in input/*
    nvidia-smi --query-gpu=index,name,driver_version,uuid --format=csv,noheader
} > run_manifest.txt

/usr/bin/time -p -o runtime_compute.txt \
    mpiexec -n "$mpi_ranks" "$opalx" track12_timed.in --info 4 \
    --kokkos-map-device-id-by=mpi_rank > opalx.log 2>&1
test -s timing.dat
date -u +%Y-%m-%dT%H:%M:%SZ > completed
