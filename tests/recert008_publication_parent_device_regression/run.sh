#!/usr/bin/env bash

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
perl "$ROOT/tests/recert008_publication_parent_device_regression/publication_parent_property_test.pl"
python3 "$ROOT/tests/recert008_publication_parent_device_regression/controller_regression.py" "$ROOT"
echo "[PASS] RECERT-008 publication-parent regression runner"
