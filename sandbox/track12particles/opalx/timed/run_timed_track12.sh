#!/usr/bin/env bash
set -euo pipefail

run_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "${run_dir}/../../../../" && pwd)"
python_bin="${OPALX_PYTHON:-${HOME}/.venv-h6/bin/python}"
opalx_exe="${OPALX_EXE:-${repo_dir}/build_openmp/src/opalx}"
primary_macroparticles="${PRIMARY_MACROPARTICLES:-400000}"
nx="${NX:-1024}"
ny="${NY:-128}"
nz="${NZ:-128}"

"${python_bin}" "${run_dir}/prepare_timed_track12.py" \
    --primary-macroparticles "${primary_macroparticles}" \
    --nx "${nx}" \
    --ny "${ny}" \
    --nz "${nz}"

if [[ ! -x "${opalx_exe}" ]]; then
    echo "OPALX executable not found: ${opalx_exe}" >&2
    exit 1
fi

(
    cd "${run_dir}"
    "${opalx_exe}" track12_timed.in 2>&1 | tee track12_timed.out
)

"${python_bin}" "${run_dir}/compare_timed_track12.py"
