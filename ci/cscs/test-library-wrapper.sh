#!/usr/bin/env bash
set -euo pipefail

if [[ $# -eq 0 ]]; then
  echo "Usage: $0 command [arguments...]" >&2
  exit 2
fi

# Run after srun activates the uenv view, which can prepend its own Kokkos libraries.
# The build job copies this script into the artifact root so it remains relocatable.
build_path=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
export LD_LIBRARY_PATH="$build_path/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
exec "$@"
