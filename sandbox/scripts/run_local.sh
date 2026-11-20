#!/usr/bin/env bash
set -euo pipefail
experiment="$1"
shift
scripts_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python_bin="${OPALX_PYTHON:-${HOME}/.venv-h6/bin/python}"
export PYTHONDONTWRITEBYTECODE=1
exec "$python_bin" "$scripts_dir/run_experiment.py" --experiment "$experiment" "$@"
