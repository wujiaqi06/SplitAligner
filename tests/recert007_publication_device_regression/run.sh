#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"

perl "$ROOT/tests/recert007_publication_device_regression/publication_device_property_test.pl"
python3 "$ROOT/tests/recert007_publication_device_regression/device_regression.py" "$ROOT"

echo "[PASS] RECERT-007 publication-device regression runner"
