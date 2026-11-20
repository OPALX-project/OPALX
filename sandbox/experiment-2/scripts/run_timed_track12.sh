#!/usr/bin/env bash
# Compatibility entry point; use ../run.sh for the documented interface.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec bash "$here/../run.sh" "$@"
