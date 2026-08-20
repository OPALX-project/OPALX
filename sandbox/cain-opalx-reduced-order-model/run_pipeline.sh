#!/usr/bin/env bash
set -euo pipefail

workflow_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_dir="$(cd "${workflow_dir}/../.." && pwd)"
python_bin="${OPALX_PYTHON:-${HOME}/.venv-h6/bin/python}"
build_dir="${OPALX_BUILD_DIR:-${repo_dir}/build_openmp}"

"${python_bin}" "${workflow_dir}/python/convert_cain_pairs.py"
"${python_bin}" "${workflow_dir}/python/validate_pair_input.py"

emitted_test="${build_dir}/unit_tests/Distribution/TestEmittedFromFile"
if [[ -x "${emitted_test}" ]]; then
    "${emitted_test}"
else
    echo "Focused C++ test not found at ${emitted_test}; build TestEmittedFromFile first." >&2
fi

if [[ "${RUN_OPALX:-0}" == "1" ]]; then
    opalx_exe="${OPALX_EXE:-${build_dir}/src/opalx}"
    if [[ ! -x "${opalx_exe}" ]]; then
        echo "OPALX executable not found at ${opalx_exe}" >&2
        exit 1
    fi
    (
        cd "${workflow_dir}/opalx"
        "${opalx_exe}" gamma_gamma_pairs_timed.in
    )
    "${python_bin}" "${workflow_dir}/python/validate_opalx_output.py"
fi
