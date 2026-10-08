#!/usr/bin/env bash
# Representation- and input-variant regression battery (see RESULTS.md).
set -euo pipefail
TEST_DIR="$(cd "$(dirname "$0")" && pwd)"
python3 "$TEST_DIR/variant_battery.py" "$@"
