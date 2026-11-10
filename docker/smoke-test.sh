#!/usr/bin/env bash
#
# Runs smoke regression tests inside an OPALX runtime image and checks the
# results against the reference .stat files, like the CPU Serial CI job.
#
# Usage: docker/smoke-test.sh <image> <regression-tests-x checkout> [test ...]
#
# Without test names, the CI subset is run. Run directories are kept in a
# temporary directory that is printed at the end.

set -euo pipefail

if [ $# -lt 2 ]; then
    echo "usage: $0 <image> <regression-tests-x checkout> [test ...]" >&2
    exit 2
fi

image=$1
tests_root=$(cd "$2/RegressionTests" && pwd)
shift 2
[ $# -gt 0 ] || set -- Drift-1-fromfile Drift-4-multi-emit-open FodoCell-fromfile

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
runs=$(mktemp -d)
failed=()

for test in "$@"; do
    echo "=== $test"
    run_dir="$runs/$test"
    mkdir -p "$run_dir"
    # Everything but the reference results, including folders such as fieldmaps/
    find "$tests_root/$test" -mindepth 1 -maxdepth 1 ! -name reference -exec cp -r {} "$run_dir/" \;

    # The references come from a serial build. Emission draws random numbers
    # per OpenMP thread, so more threads give a different particle sample.
    if docker run --rm --user "$(id -u):$(id -g)" -e OMP_NUM_THREADS=1 \
           -v "$run_dir:/data" "$image" "$test.in" --info 2 \
       && python3 "$here/compare_stat.py" "$tests_root/$test" "$run_dir" "$test"; then
        echo "=== $test: PASS"
    else
        echo "=== $test: FAIL"
        failed+=("$test")
    fi
done

echo "Run directories: $runs"
if [ ${#failed[@]} -gt 0 ]; then
    echo "Failed: ${failed[*]}"
    exit 1
fi
echo "All $# tests passed"
