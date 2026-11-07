#!/bin/sh
#
# Entry point of the OPALX runtime image.
#
#   docker run ... opalx:cpu case.in --info 2   runs opalx on case.in under mpirun
#   docker run ... opalx:cpu bash               runs any command found on PATH
#
# OPALX_NP sets the number of MPI ranks (default 1). Binding is disabled so
# that the OpenMP threads of a rank are not all pinned to a single core.

set -eu

case "${1-}" in
    "" | -* | */*) ;;
    *) command -v "$1" >/dev/null 2>&1 && exec "$@" ;;
esac

exec mpirun -np "${OPALX_NP:-1}" --bind-to none opalx "$@"
