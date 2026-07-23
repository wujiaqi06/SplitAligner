#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
TEST_DIR="$ROOT_DIR/tests/confirm_na_structure_regression"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

perl "$ROOT_DIR/scripts/confirm_na_structure.pl" \
  --fix "$TEST_DIR/fix.matrix_with_fuse.txt" \
  --free "$TEST_DIR/free.matrix_with_fuse.txt" \
  --fix_state "$TEST_DIR/fix.primitive_state.tsv" \
  --free_state "$TEST_DIR/free.primitive_state.tsv" \
  -o "$TMP_DIR/out" >/dev/null 2>&1

diff -u "$TEST_DIR/expected.fix.na_classified.txt" "$TMP_DIR/out.fix.na_classified.txt"
diff -u "$TEST_DIR/expected.free.na_classified.txt" "$TMP_DIR/out.free.na_classified.txt"

python3 - <<'PY' "$TMP_DIR/out.fix.na_classified.txt" "$TMP_DIR/out.free.na_classified.txt"
import sys
from pathlib import Path

fix_path = Path(sys.argv[1])
free_path = Path(sys.argv[2])

def load_row(path):
    with path.open() as fh:
        header = next(fh).rstrip("\n").split("\t")
        row = next(fh).rstrip("\n").split("\t")
    return dict(zip(header[1:], row[1:]))

fix = load_row(fix_path)
free = load_row(free_path)

old_free = {
    "B1": "NA_struct",
    "B2": "NA_topo",
    "B3": "NA_topo",
    "B4": "NA_fuse",
    "B5": "NA_topo",
    "B6": "NA_topo",
    "B7": "NA_topo",
    "B8": "NA_topo",
    "B9": "NA_topo",
}

print("Regression test passed.")
print("New free output:", free)
print("Old buggy free output would have been:", old_free)
print("Fixed output:", fix)
PY
