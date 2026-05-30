#!/usr/bin/env zsh
set -euo pipefail

script_dir="${0:A:h}"
repo_dir="${script_dir:h:h}"

opalx_exe="${OPALX_EXE:-${repo_dir}/build/src/opalx}"
old_opal_profile="${OLD_OPAL_PROFILE:-/Users/adelmann/OPAL-2022.1/etc/profile.d/opal.sh}"
run_root="${1:-/tmp/opalx-rbend-compare}"

if [[ ! -x "${opalx_exe}" ]]; then
    echo "OPALX executable not found or not executable: ${opalx_exe}" >&2
    exit 1
fi

if [[ ! -r "${old_opal_profile}" ]]; then
    echo "Old OPAL profile not readable: ${old_opal_profile}" >&2
    exit 1
fi

python_exe="${PYTHON:-}"
if [[ -z "${python_exe}" ]]; then
    for candidate in \
        "/Users/adelmann/.venv-h6/bin/python" \
        "python3"; do
        if command -v "${candidate}" >/dev/null 2>&1; then
            python_exe="${candidate}"
            break
        fi
    done
fi

if [[ -z "${python_exe}" ]]; then
    echo "No Python interpreter found. Set PYTHON=/path/to/python." >&2
    exit 1
fi

opalx_dir="${run_root}/opalx"
opal_dir="${run_root}/opal-2022"
compare_dir="${run_root}/comparison"

mkdir -p "${opalx_dir}" "${opal_dir}" "${compare_dir}"

cp "${script_dir}/test-rbend-1.in" \
   "${script_dir}/test-rbend-1_distribution.txt" \
   "${opalx_dir}/"
cp "${script_dir}/test-rbend-1_opal.in" \
   "${script_dir}/test-rbend-1_opal_distribution.txt" \
   "${opal_dir}/"

(
    cd "${opalx_dir}"
    "${opalx_exe}" test-rbend-1.in
)

(
    set +u
    source "${old_opal_profile}"
    set -u
    cd "${opal_dir}"
    opal test-rbend-1_opal.in
)

MPLBACKEND=Agg \
MPLCONFIGDIR=/tmp/opalx-matplotlib \
XDG_CACHE_HOME=/tmp/opalx-cache \
PYTHONFAULTHANDLER=1 \
"${python_exe}" "${script_dir}/compare_opal_opalx.py" \
    --opalx-dir "${opalx_dir}" \
    --opal-dir "${opal_dir}" \
    --out-dir "${compare_dir}" \
    --plots

echo
echo "Comparison written to ${compare_dir}"
