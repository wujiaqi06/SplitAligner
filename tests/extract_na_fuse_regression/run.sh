#!/usr/bin/env bash

set -euo pipefail

TEST_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$TEST_DIR/../.." && pwd)"
OUTPUT_FILE="$TEST_DIR/input.matrix_with_fuse.na_fuse.txt"

perl "$REPO_ROOT/scripts/extract_na_fuse.pl" \
  -i "$TEST_DIR/input.matrix_with_fuse.txt" \
  -o "$OUTPUT_FILE"

diff -u "$TEST_DIR/expected.na_fuse.txt" "$OUTPUT_FILE"

echo "[PASS] extract_na_fuse regression test"
