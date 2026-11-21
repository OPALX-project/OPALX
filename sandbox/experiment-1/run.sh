#!/usr/bin/env bash
set -euo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
exec bash "$here/../scripts/run_local.sh" 1 "$@"
