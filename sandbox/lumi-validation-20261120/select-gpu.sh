#!/bin/bash
set -euo pipefail
# srun --gpus-per-task supplies a per-task ROCr mask. A second rank-based
# HIP mask would hide rank 1's sole exposed device. Retain the historical
# rank selector only for allocations without a Slurm ROCr mask.
if [[ -n ${ROCR_VISIBLE_DEVICES:-} ]]; then
    unset HIP_VISIBLE_DEVICES
else
    export HIP_VISIBLE_DEVICES=${SLURM_LOCALID:-0}
fi
printf 'GPU binding: rank=%s ROCR_VISIBLE_DEVICES=%s HIP_VISIBLE_DEVICES=%s\n' \
    "${SLURM_LOCALID:-0}" "${ROCR_VISIBLE_DEVICES:-unset}" "${HIP_VISIBLE_DEVICES:-unset}" >&2
exec "$@"
