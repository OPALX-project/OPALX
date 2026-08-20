#!/usr/bin/env bash
# Submit the fixed 400k pair-4 field probe to four Merlin6 A100 GPUs.

set -euo pipefail

STAGE="submit"
if [[ $# -gt 0 && "$1" =~ ^(submit|run)$ ]]; then
    STAGE="$1"
    shift
fi

RUN_ROOT=""
OPALX=""
CLUSTER="gmerlin6"
PARTITION="gwendolen"
ACCOUNT="gwendolen"
OMP_THREADS=4
CASE_NAME="production_rect_1024x128_fixed400k"

usage() {
    echo "Usage: $0 [submit|run] --run-root PATH --opalx PATH [options]"
    echo "Options:"
    echo "  --cluster NAME    Slurm cluster (default: gmerlin6)"
    echo "  --partition NAME  Slurm partition (default: gwendolen)"
    echo "  --account NAME    Slurm account (default: gwendolen)"
    echo "  --omp-threads N   CPU threads per MPI rank (default: 4)"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root) RUN_ROOT="$2"; shift 2 ;;
        --opalx) OPALX="$2"; shift 2 ;;
        --cluster) CLUSTER="$2"; shift 2 ;;
        --partition) PARTITION="$2"; shift 2 ;;
        --account) ACCOUNT="$2"; shift 2 ;;
        --omp-threads) OMP_THREADS="$2"; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

if [[ -z "$RUN_ROOT" || -z "$OPALX" ]]; then
    usage >&2
    exit 2
fi
if [[ ! -x "$OPALX" ]]; then
    echo "OPALX executable is not available: $OPALX" >&2
    exit 1
fi

case_dir="$RUN_ROOT/$CASE_NAME"
if [[ ! -f "$case_dir/pair4_field_probe.in" || \
      ! -f "$RUN_ROOT/input/primary_fixed.fromfile" ]]; then
    echo "Prepared field-probe case is missing below $RUN_ROOT" >&2
    exit 1
fi

COMMON_ARGS=(
    --run-root "$RUN_ROOT"
    --opalx "$OPALX"
    --cluster "$CLUSTER"
    --partition "$PARTITION"
    --account "$ACCOUNT"
    --omp-threads "$OMP_THREADS"
)

if [[ "$STAGE" == "submit" ]]; then
    mkdir -p "$RUN_ROOT/slurm"
    script_path="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
    job_raw="$(sbatch --parsable --clusters="$CLUSTER" --account="$ACCOUNT" \
        --job-name=track12-pair4-field --partition="$PARTITION" \
        --nodes=1 --ntasks=4 --cpus-per-task="$OMP_THREADS" --gpus-per-node=4 \
        --time=00:30:00 --output="$RUN_ROOT/slurm/field_%j.out" \
        "$script_path" run "${COMMON_ARGS[@]}")"
    job_id="${job_raw%%;*}"
    printf 'Run root: %s\n' "$RUN_ROOT"
    printf 'Job: %s\n' "$job_id"
    printf 'Monitor: squeue -M %s -j %s\n' "$CLUSTER" "$job_id"
    exit 0
fi

if [[ "$STAGE" != "run" ]]; then
    echo "Unknown stage: $STAGE" >&2
    exit 2
fi

cd "$case_dir"
export OMP_NUM_THREADS="$OMP_THREADS"
{
    printf 'utc_start=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'host=%s\n' "$(hostname)"
    printf 'slurm_job_id=%s\n' "${SLURM_JOB_ID:-}"
    printf 'mpi_ranks=%s\n' "${SLURM_NTASKS:-4}"
    printf 'omp_threads=%s\n' "$OMP_THREADS"
    printf 'opalx=%s\n' "$OPALX"
    sha256sum "$OPALX" pair4_field_probe.in input/pair4_probe.fromfile \
        ../input/primary_fixed.fromfile
    nvidia-smi --query-gpu=index,name,uuid --format=csv,noheader
} > run_manifest.txt

/usr/bin/time -p -o runtime_compute.txt \
    mpiexec -n "${SLURM_NTASKS:-4}" \
    "$OPALX" pair4_field_probe.in --kokkos-map-device-id-by=mpi_rank \
    > opalx.log 2>&1

test -f pair4_field_probe_c1.h5
date -u +%Y-%m-%dT%H:%M:%SZ > completed
