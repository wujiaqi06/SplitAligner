#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -P -- "$(dirname "$0")" && pwd -P)"
PYTHONDONTWRITEBYTECODE=1 python3 "$SCRIPT_DIR/archive_reset_regression.py"
