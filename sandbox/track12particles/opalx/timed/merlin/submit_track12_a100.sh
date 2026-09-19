#!/usr/bin/env bash
# Submit the exact-timed track12 CAIN comparison to Merlin6 A100 GPUs.

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
MPI_RANKS=4
PRODUCTION_TIME="04:00:00"

usage() {
    echo "Usage: $0 [submit|run] --run-root PATH --opalx PATH [options]"
    echo "Options:"
    echo "  --case NAME       Internal run stage: smoke or production"
    echo "  --cluster NAME    Slurm cluster (default: gmerlin6)"
    echo "  --partition NAME  Slurm partition (default: gwendolen)"
    echo "  --account NAME    Slurm account (default: gwendolen)"
    echo "  --omp-threads N   CPU threads per MPI rank (default: 4)"
    echo "  --mpi-ranks N     MPI ranks and A100 GPUs (default: 4)"
    echo "  --production-time TIME  Production Slurm limit (default: 04:00:00)"
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
        --mpi-ranks) MPI_RANKS="$2"; shift 2 ;;
        --production-time) PRODUCTION_TIME="$2"; shift 2 ;;
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
    --mpi-ranks "$MPI_RANKS"
    --production-time "$PRODUCTION_TIME"
)

if [[ "$STAGE" == "submit" ]]; then
    production="$RUN_ROOT/production"
    smoke="$RUN_ROOT/smoke"
    slurm="$RUN_ROOT/slurm"
    if [[ ! -f "$production/track12_timed.in" || ! -d "$production/input" ]]; then
        echo "Prepared production case is missing below $production" >&2
        exit 1
    fi
    mkdir -p "$smoke" "$slurm"
    ln -sfn "$production/input" "$smoke/input"
    sed 's/REAL n_steps = 1501;/REAL n_steps = 3;/' \
        "$production/track12_timed.in" > "$smoke/track12_timed.in"

    script_path="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
    submission_id="$(date -u +%Y%m%dT%H%M%SZ)_$$"
    submitted_script="$slurm/submit_track12_${submission_id}.sh"
    cp "$script_path" "$submitted_script"
    chmod u+x "$submitted_script"
    smoke_raw="$(sbatch --parsable --clusters="$CLUSTER" --account="$ACCOUNT" \
        --job-name=track12-smoke --partition="$PARTITION" --nodes=1 --ntasks="$MPI_RANKS" \
        --cpus-per-task="$OMP_THREADS" --gpus-per-node="$MPI_RANKS" --time=00:15:00 \
        --output="$slurm/smoke_%j.out" \
        "$submitted_script" run "${COMMON_ARGS[@]}" --case smoke)"
    smoke_job="${smoke_raw%%;*}"
    production_raw="$(sbatch --parsable --clusters="$CLUSTER" --account="$ACCOUNT" \
        --dependency="afterok:$smoke_job" \
        --job-name=track12-production --partition="$PARTITION" --nodes=1 --ntasks="$MPI_RANKS" \
        --cpus-per-task="$OMP_THREADS" --gpus-per-node="$MPI_RANKS" \
        --time="$PRODUCTION_TIME" \
        --output="$slurm/production_%j.out" \
        "$submitted_script" run "${COMMON_ARGS[@]}" --case production)"
    production_job="${production_raw%%;*}"

    printf 'Run root: %s\n' "$RUN_ROOT"
    printf 'Jobs: smoke=%s production=%s\n' "$smoke_job" "$production_job"
    printf 'Submitted script: %s\n' "$submitted_script"
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
    printf 'mpi_ranks=%s\n' "${SLURM_NTASKS:-$MPI_RANKS}"
    printf 'omp_threads=%s\n' "$OMP_THREADS"
    printf 'requested_production_time=%s\n' "$PRODUCTION_TIME"
    printf 'opalx=%s\n' "$OPALX"
    sha256sum "$OPALX" track12_timed.in input/*
    nvidia-smi --query-gpu=index,name,uuid --format=csv,noheader
} > run_manifest.txt

# IPPL's timing.dat writer inherits the --info level.  Its default level zero
# creates the file but suppresses every row, so request level four explicitly.
/usr/bin/time -p -o runtime_compute.txt \
    mpiexec -n "${SLURM_NTASKS:-$MPI_RANKS}" \
    "$OPALX" track12_timed.in --info 4 --kokkos-map-device-id-by=mpi_rank \
    > opalx.log 2>&1
test -s timing.dat
date -u +%Y-%m-%dT%H:%M:%SZ > completed
