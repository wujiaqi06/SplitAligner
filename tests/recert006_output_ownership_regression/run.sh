#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
exec python3 "$ROOT/tests/recert006_output_ownership_regression/ownership_regression.py" "$ROOT"
