#!/usr/bin/env bash
# Submit and run the BeamBeam manufactured-field convergence study on Merlin6.

set -euo pipefail

STAGE="submit"
if [[ $# -gt 0 && "$1" =~ ^(submit|single|rank|finalize)$ ]]; then
    STAGE="$1"
    shift
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT_PATH="${SCRIPT_DIR}/$(basename "${BASH_SOURCE[0]}")"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../../.." && pwd)"
RUN_ROOT=""
OPALX="/psi/home/adelmann/opalx-dev/build-opalx-a100-5a354101f-pinned/src/opalx"
PYTHON="python3"
PARTITION="gwendolen"
CLUSTER="gmerlin6"
ACCOUNT="gwendolen"
OMP_THREADS=4
RANKS=""
RANK_ONLY=false

usage() {
    echo "Usage: $0 [submit] [options]"
    echo "       $0 single|rank|finalize [internal options]"
    echo "Options:"
    echo "  --run-root PATH    Remote study output directory"
    echo "  --repo-root PATH   OPALX source tree containing this workflow"
    echo "  --opalx PATH       CUDA OPALX executable"
    echo "  --python PATH      Python used only to prepare input decks"
    echo "  --partition NAME   Slurm A100 partition (default: gwendolen)"
    echo "  --cluster NAME     Slurm federation cluster (default: gmerlin6)"
    echo "  --account NAME     Slurm account (default: gwendolen)"
    echo "  --omp-threads N    CPU threads per OPALX rank (default: 4)"
    echo "  --rank-only        Run only the fixed-primary 1/2/4-rank decomposition scan"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --run-root) RUN_ROOT="$2"; shift 2 ;;
        --repo-root) REPO_ROOT="$2"; shift 2 ;;
        --opalx) OPALX="$2"; shift 2 ;;
        --python) PYTHON="$2"; shift 2 ;;
        --partition) PARTITION="$2"; shift 2 ;;
        --cluster) CLUSTER="$2"; shift 2 ;;
        --account) ACCOUNT="$2"; shift 2 ;;
        --omp-threads) OMP_THREADS="$2"; shift 2 ;;
        --ranks) RANKS="$2"; shift 2 ;;
        --rank-only) RANK_ONLY=true; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
done

if [[ ! -x "$OPALX" ]]; then
    echo "OPALX executable is not available: $OPALX" >&2
    exit 1
fi
if [[ ! -x "$PYTHON" ]] && ! command -v "$PYTHON" >/dev/null 2>&1; then
    echo "Python is not available: $PYTHON" >&2
    exit 1
fi

REVISION="$(git -C "$REPO_ROOT" rev-parse HEAD)"
if [[ -z "$RUN_ROOT" ]]; then
    TIMESTAMP="$(date -u +%Y%m%dT%H%M%SZ)"
    RUN_ROOT="/psi/home/adelmann/opalx-dev/a100-convergence-${REVISION:0:10}-${TIMESTAMP}"
fi

COMMON_ARGS=(
    --run-root "$RUN_ROOT"
    --repo-root "$REPO_ROOT"
    --opalx "$OPALX"
    --python "$PYTHON"
    --partition "$PARTITION"
    --cluster "$CLUSTER"
    --account "$ACCOUNT"
    --omp-threads "$OMP_THREADS"
)
if [[ "$RANK_ONLY" == true ]]; then
    COMMON_ARGS+=(--rank-only)
fi

run_case() {
    local label="$1"
    shift
    local case_dir="$RUN_ROOT/cases/$label/3sigma"
    local input_file="rigid_fields_3sigma.in"
    if [[ -f "$case_dir/completed" ]]; then
        echo "cached $label"
        return 0
    fi
    echo "running $label"
    (
        cd "$case_dir"
        export OMP_NUM_THREADS="$OMP_THREADS"
        /usr/bin/time -p -o runtime_compute.txt \
            "$@" "$OPALX" "$input_file" \
            --kokkos-map-device-id-by=mpi_rank >opalx.log 2>&1
        touch completed
    )
}

finalize_study() {
    local missing=0
    local -a label_files
    local completion_marker="study-completed"
    if [[ "$RANK_ONLY" == true ]]; then
        label_files=("$RUN_ROOT"/mpi_rank{1,2,4}_labels.txt)
        completion_marker="rank-study-completed"
    else
        label_files=("$RUN_ROOT"/rank{1,2,4}_labels.txt)
    fi
    while IFS= read -r expected; do
        [[ -n "$expected" ]] || continue
        if [[ ! -f "$RUN_ROOT/cases/$expected/3sigma/completed" ]]; then
            echo "missing completion marker for $expected" >&2
            missing=1
        fi
    done < <(cat "${label_files[@]}")
    if [[ "$missing" -ne 0 ]]; then
        return 1
    fi
    date -u +%Y-%m-%dT%H:%M:%SZ > "$RUN_ROOT/$completion_marker"
    echo "study complete: $RUN_ROOT"
}

if [[ "$STAGE" == "submit" ]]; then
    mkdir -p "$RUN_ROOT/slurm"
    "$PYTHON" "$REPO_ROOT/sandbox/reduced-order-model/python/prepare_a100_convergence.py" \
        --run-root "$RUN_ROOT" --opalx "$OPALX"

    single_resources=(--nodes=1 --ntasks=4 --gpus-per-node=4 --time=02:00:00)
    if [[ "$RANK_ONLY" == true ]]; then
        single_resources=(--nodes=1 --ntasks=1 --gpus-per-node=1 --time=00:30:00)
    fi
    single_raw="$(sbatch --parsable --clusters="$CLUSTER" --account="$ACCOUNT" \
        --job-name=bbconv1 --partition="$PARTITION" "${single_resources[@]}" \
        --cpus-per-task="$OMP_THREADS" \
        --output="$RUN_ROOT/slurm/single_%j.out" \
        "$SCRIPT_PATH" single "${COMMON_ARGS[@]}")"
    single_job="${single_raw%%;*}"
    rank2_raw="$(sbatch --parsable --clusters="$CLUSTER" --account="$ACCOUNT" \
        --dependency="afterok:$single_job" \
        --job-name=bbconv2 --partition="$PARTITION" --nodes=1 --ntasks=2 \
        --cpus-per-task="$OMP_THREADS" --gpus-per-node=2 --time=00:30:00 \
        --output="$RUN_ROOT/slurm/rank2_%j.out" \
        "$SCRIPT_PATH" rank "${COMMON_ARGS[@]}" --ranks 2)"
    rank2_job="${rank2_raw%%;*}"
    rank4_raw="$(sbatch --parsable --clusters="$CLUSTER" --account="$ACCOUNT" \
        --dependency="afterok:$rank2_job" \
        --job-name=bbconv4 --partition="$PARTITION" --nodes=1 --ntasks=4 \
        --cpus-per-task="$OMP_THREADS" --gpus-per-node=4 --time=00:30:00 \
        --output="$RUN_ROOT/slurm/rank4_%j.out" \
        "$SCRIPT_PATH" rank "${COMMON_ARGS[@]}" --ranks 4)"
    rank4_job="${rank4_raw%%;*}"

    printf 'Study root: %s\n' "$RUN_ROOT"
    printf 'Jobs: single-rank=%s, two-rank=%s, four-rank=%s\n' \
        "$single_job" "$rank2_job" "$rank4_job"
    printf 'Monitor: squeue -M %s -j %s,%s,%s\n' \
        "$CLUSTER" "$single_job" "$rank2_job" "$rank4_job"
    exit 0
fi

if [[ ! -f "$RUN_ROOT/study_manifest.csv" ]]; then
    echo "Study has not been prepared: $RUN_ROOT" >&2
    exit 1
fi

if [[ "$STAGE" == "finalize" ]]; then
    finalize_study
    exit 0
fi

if [[ "$STAGE" == "single" ]]; then
    baseline_label="$(head -n 1 "$RUN_ROOT/rank_baseline_label.txt")"
    run_case "$baseline_label" srun --exclusive --nodes=1 --ntasks=1 \
        --cpus-per-task="$OMP_THREADS" --gpus-per-task=1 \
        --gpu-bind=single:1 </dev/null
    if [[ "$RANK_ONLY" == true ]]; then
        exit 0
    fi
    pids=()
    labels=()
    status=0
    while IFS= read -r label; do
        [[ -n "$label" ]] || continue
        [[ "$label" != "$baseline_label" ]] || continue
        run_case "$label" srun --exclusive --nodes=1 --ntasks=1 \
            --cpus-per-task="$OMP_THREADS" --gpus-per-task=1 \
            --gpu-bind=single:1 </dev/null &
        pids+=("$!")
        labels+=("$label")
        if [[ ${#pids[@]} -eq 4 ]]; then
            for index in "${!pids[@]}"; do
                if ! wait "${pids[$index]}"; then
                    echo "failed ${labels[$index]}" >&2
                    status=1
                fi
            done
            pids=()
            labels=()
        fi
    done < "$RUN_ROOT/rank1_labels.txt"
    for index in "${!pids[@]}"; do
        if ! wait "${pids[$index]}"; then
            echo "failed ${labels[$index]}" >&2
            status=1
        fi
    done
    exit "$status"
fi

if [[ "$STAGE" == "rank" ]]; then
    if [[ "$RANKS" != "2" && "$RANKS" != "4" ]]; then
        echo "rank stage requires --ranks 2 or --ranks 4" >&2
        exit 2
    fi
    label="$(head -n 1 "$RUN_ROOT/mpi_rank${RANKS}_labels.txt")"
    run_case "$label" mpiexec -n "$RANKS"
    if [[ "$RANKS" == "4" ]]; then
        finalize_study
    fi
fi
