#!/usr/bin/env bash
#
# Opens a shell (or runs a command) in the OPALX development image.
#
# The repository is mounted at /work. Set OPALX_WORK to mount a different
# directory instead, e.g. a parent directory that also holds build and run
# directories. Files are created with your user and group IDs.
#
#   docker build -f docker/Dockerfile --target dev -t opalx-dev .
#   docker/dev-shell.sh
#   docker/dev-shell.sh cmake -S /work -B /work/build_openmp -DPLATFORMS=OPENMP

set -euo pipefail

repo=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
work=${OPALX_WORK:-$repo}
image=${OPALX_DEV_IMAGE:-opalx-dev}

tty_flags=-i
[ -t 0 ] && tty_flags=-it

# --shm-size: Docker's 64 MB default is too small for MPI shared-memory transport
exec docker run --rm $tty_flags \
    --user "$(id -u):$(id -g)" -e HOME=/tmp \
    --shm-size=8g \
    -v "$work:/work" -w /work \
    "$image" "${@:-bash}"
