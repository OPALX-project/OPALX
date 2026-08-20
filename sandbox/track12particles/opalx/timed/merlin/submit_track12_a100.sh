#!/usr/bin/env bash
# Submit the exact-timed track12 CAIN comparison to four Merlin6 A100 GPUs.

set -euo pipefail

STAGE="submit"
if [[ $# -gt 0 && "$1" =~ ^(submit|run)$ ]]; then
    STAGE="$1"
    shift
fi

RUN_ROOT=""
OPALX=""
CASE_NAME=""
CLUSTER="gmerlin6"
PARTITION="gwendolen"
ACCOUNT="gwendolen"
OMP_THREADS=4

usage() {
    echo "Usage: $0 [submit|run] --run-root PATH --opalx PATH [options]"
    echo "Options:"
    echo "  --case NAME       Internal run stage: smoke or production"
    echo "  --cluster NAME    Slurm cluster (default: gmerlin6)"
    echo "  --partition NAME  Slurm partition (default: gwendolen)"
    echo "  --account NAME    Slurm account (default: gwendolen)"
    echo "  --omp-threads N   CPU threads per MPI rank (default: 4)"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root) RUN_ROOT="$2"; shift 2 ;;
        --opalx) OPALX="$2"; shift 2 ;;
        --case) CASE_NAME="$2"; shift 2 ;;
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

COMMON_ARGS=(
    --run-root "$RUN_ROOT"
    --opalx "$OPALX"
    --cluster "$CLUSTER"
    --partition "$PARTITION"
    --account "$ACCOUNT"
    --omp-threads "$OMP_THREADS"
)

if [[ "$STAGE" == "submit" ]]; then
    production="$RUN_ROOT/production"
    smoke="$RUN_ROOT/smoke"
    slurm="$RUN_ROOT/../slurm"
    if [[ ! -f "$production/track12_timed.in" || ! -d "$production/input" ]]; then
        echo "Prepared production case is missing below $production" >&2
        exit 1
    fi
    mkdir -p "$smoke" "$slurm"
    ln -sfn "$production/input" "$smoke/input"
    sed 's/REAL n_steps = 1501;/REAL n_steps = 3;/' \
        "$production/track12_timed.in" > "$smoke/track12_timed.in"

    script_path="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
    smoke_raw="$(sbatch --parsable --clusters="$CLUSTER" --account="$ACCOUNT" \
        --job-name=track12-smoke --partition="$PARTITION" --nodes=1 --ntasks=4 \
        --cpus-per-task="$OMP_THREADS" --gpus-per-node=4 --time=00:15:00 \
        --output="$slurm/smoke_%j.out" \
        "$script_path" run "${COMMON_ARGS[@]}" --case smoke)"
    smoke_job="${smoke_raw%%;*}"
    production_raw="$(sbatch --parsable --clusters="$CLUSTER" --account="$ACCOUNT" \
        --dependency="afterok:$smoke_job" \
        --job-name=track12-production --partition="$PARTITION" --nodes=1 --ntasks=4 \
        --cpus-per-task="$OMP_THREADS" --gpus-per-node=4 --time=04:00:00 \
        --output="$slurm/production_%j.out" \
        "$script_path" run "${COMMON_ARGS[@]}" --case production)"
    production_job="${production_raw%%;*}"

    printf 'Run root: %s\n' "$RUN_ROOT"
    printf 'Jobs: smoke=%s production=%s\n' "$smoke_job" "$production_job"
    printf 'Monitor: squeue -M %s -j %s,%s\n' \
        "$CLUSTER" "$smoke_job" "$production_job"
    exit 0
fi

if [[ "$STAGE" != "run" || ("$CASE_NAME" != "smoke" && "$CASE_NAME" != "production") ]]; then
    echo "run stage requires --case smoke or --case production" >&2
    exit 2
fi

case_dir="$RUN_ROOT/$CASE_NAME"
cd "$case_dir"
export OMP_NUM_THREADS="$OMP_THREADS"
{
    printf 'case=%s\n' "$CASE_NAME"
    printf 'utc_start=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    printf 'host=%s\n' "$(hostname)"
    printf 'slurm_job_id=%s\n' "${SLURM_JOB_ID:-}"
    printf 'mpi_ranks=%s\n' "${SLURM_NTASKS:-4}"
    printf 'omp_threads=%s\n' "$OMP_THREADS"
    printf 'opalx=%s\n' "$OPALX"
    sha256sum "$OPALX" track12_timed.in input/*
    nvidia-smi --query-gpu=index,name,uuid --format=csv,noheader
} > run_manifest.txt

/usr/bin/time -p -o runtime_compute.txt \
    mpiexec -n "${SLURM_NTASKS:-4}" \
    "$OPALX" track12_timed.in --kokkos-map-device-id-by=mpi_rank \
    > opalx.log 2>&1
date -u +%Y-%m-%dT%H:%M:%SZ > completed
