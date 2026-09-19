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
if [[ $# -gt 1 ]]; then
    echo "Usage: $0 [new-output-directory]" >&2
    exit 2
fi
output_dir="${1:-}"
if [[ -z "$output_dir" ]]; then
    mkdir -p "${run_dir}/runs"
    output_dir="$(mktemp -d "${run_dir}/runs/track12.XXXXXX")"
fi

"${python_bin}" "${run_dir}/prepare_timed_track12.py" \
    --output-dir "${output_dir}" \
    --primary-macroparticles "${primary_macroparticles}" \
    --nx "${nx}" \
    --ny "${ny}" \
    --nz "${nz}"

if [[ ! -x "${opalx_exe}" ]]; then
    echo "OPALX executable not found: ${opalx_exe}" >&2
    exit 1
fi

(
    cd "${output_dir}"
    # IPPL suppresses timing.dat contents at its default --info level zero.
    "${opalx_exe}" track12_timed.in --info 4 2>&1 | tee track12_timed.out
    test -s timing.dat
)

"${python_bin}" "${run_dir}/compare_timed_track12.py" --run-dir "${output_dir}"
