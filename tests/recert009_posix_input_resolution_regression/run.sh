#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"

perl "$ROOT/tests/recert009_posix_input_resolution_regression/input_resolution_property_test.pl"
python3 "$ROOT/tests/recert009_posix_input_resolution_regression/controller_regression.py" "$ROOT"
bash "$ROOT/tests/recert009_posix_input_resolution_regression/symlink_retarget_controller.sh"

echo "[PASS] RECERT-009 POSIX existing-input resolution regression runner"
