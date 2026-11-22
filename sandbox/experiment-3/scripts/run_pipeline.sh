#!/usr/bin/env bash
# Compatibility entry point; scientific execution no longer runs regression tests.
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec bash "$here/../run.sh" "$@"
