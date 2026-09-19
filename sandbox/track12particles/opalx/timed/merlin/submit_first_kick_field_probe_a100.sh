#!/usr/bin/env bash
# Submit the fixed-source pair-4 field probe to Merlin6 A100 GPUs.

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
MPI_RANKS=4
CASES="production_rect_1024x128_fixed400k"

usage() {
    echo "Usage: $0 [submit|run] --run-root PATH --opalx PATH [options]"
    echo "Options:"
    echo "  --cluster NAME    Slurm cluster (default: gmerlin6)"
    echo "  --partition NAME  Slurm partition (default: gwendolen)"
    echo "  --account NAME    Slurm account (default: gwendolen)"
    echo "  --omp-threads N   CPU threads per MPI rank (default: 4)"
    echo "  --mpi-ranks N     MPI ranks and A100 GPUs (default: 4)"
    echo "  --cases LIST      Comma-separated prepared case names"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root) RUN_ROOT="$2"; shift 2 ;;
        --opalx) OPALX="$2"; shift 2 ;;
        --cluster) CLUSTER="$2"; shift 2 ;;
        --partition) PARTITION="$2"; shift 2 ;;
        --account) ACCOUNT="$2"; shift 2 ;;
        --omp-threads) OMP_THREADS="$2"; shift 2 ;;
        --mpi-ranks) MPI_RANKS="$2"; shift 2 ;;
        --cases) CASES="$2"; shift 2 ;;
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

IFS=',' read -r -a CASE_NAMES <<< "$CASES"
for case_name in "${CASE_NAMES[@]}"; do
    case_dir="$RUN_ROOT/$case_name"
    if [[ ! -f "$case_dir/pair4_field_probe.in" || \
          ! -f "$RUN_ROOT/input/primary_fixed.fromfile" ]]; then
        echo "Prepared field-probe case is missing: $case_dir" >&2
        exit 1
    fi
done

COMMON_ARGS=(
    --run-root "$RUN_ROOT"
    --opalx "$OPALX"
    --cluster "$CLUSTER"
    --partition "$PARTITION"
    --account "$ACCOUNT"
    --omp-threads "$OMP_THREADS"
    --mpi-ranks "$MPI_RANKS"
    --cases "$CASES"
)

if [[ "$STAGE" == "submit" ]]; then
    mkdir -p "$RUN_ROOT/slurm"
    script_path="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/$(basename "${BASH_SOURCE[0]}")"
    submission_id="$(date -u +%Y%m%dT%H%M%SZ)_$$"
    submitted_script="$RUN_ROOT/slurm/submit_field_${submission_id}.sh"
    cp "$script_path" "$submitted_script"
    chmod u+x "$submitted_script"
    job_raw="$(sbatch --parsable --clusters="$CLUSTER" --account="$ACCOUNT" \
        --job-name=track12-pair4-field --partition="$PARTITION" \
        --nodes=1 --ntasks="$MPI_RANKS" --cpus-per-task="$OMP_THREADS" \
        --gpus-per-node="$MPI_RANKS" \
        --time=00:30:00 --output="$RUN_ROOT/slurm/field_%j.out" \
        "$submitted_script" run "${COMMON_ARGS[@]}")"
    job_id="${job_raw%%;*}"
    printf 'Run root: %s\n' "$RUN_ROOT"
    printf 'Job: %s\n' "$job_id"
    printf 'Submitted script: %s\n' "$submitted_script"
    printf 'Monitor: squeue -M %s -j %s\n' "$CLUSTER" "$job_id"
    exit 0
fi

if [[ "$STAGE" != "run" ]]; then
    echo "Unknown stage: $STAGE" >&2
    exit 2
fi

export OMP_NUM_THREADS="$OMP_THREADS"
for case_name in "${CASE_NAMES[@]}"; do
    case_dir="$RUN_ROOT/$case_name"
    cd "$case_dir"
    {
        printf 'case=%s\n' "$case_name"
        printf 'utc_start=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
        printf 'host=%s\n' "$(hostname)"
        printf 'slurm_job_id=%s\n' "${SLURM_JOB_ID:-}"
        printf 'mpi_ranks=%s\n' "${SLURM_NTASKS:-$MPI_RANKS}"
        printf 'omp_threads=%s\n' "$OMP_THREADS"
        printf 'opalx=%s\n' "$OPALX"
        sha256sum "$OPALX" pair4_field_probe.in input/pair4_probe.fromfile \
            ../input/primary_fixed.fromfile
        nvidia-smi --query-gpu=index,name,uuid --format=csv,noheader
    } > run_manifest.txt

    # IPPL's timing.dat writer inherits the --info level.  Level zero creates
    # an empty file, so request level one and verify the result below.
    /usr/bin/time -p -o runtime_compute.txt \
        mpiexec -n "${SLURM_NTASKS:-$MPI_RANKS}" \
        "$OPALX" pair4_field_probe.in --info 1 --kokkos-map-device-id-by=mpi_rank \
        > opalx.log 2>&1

    test -f pair4_field_probe_c1.h5
    test -s timing.dat
    date -u +%Y-%m-%dT%H:%M:%SZ > completed
done
