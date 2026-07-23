#!/usr/bin/env bash

set -euo pipefail

TEST_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$TEST_DIR/../.." && pwd)"
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT
OUT_PREFIX="$TMP_DIR/final"
SUPPORT_TREE="$TMP_DIR/species_tree.support_b.nwk"
cp "$TEST_DIR/species_tree.forSplit.nwk" "$TMP_DIR/species_tree.forSplit.nwk"
cp "$TEST_DIR/species_tree.branch_map.txt" "$TMP_DIR/species_tree.branch_map.txt"

(
  cd "$TMP_DIR"
  perl "$REPO_ROOT/scripts/confirm_na_structure.pl" \
    --fix "$TEST_DIR/fix.matrix_with_fuse.txt" \
    --free "$TEST_DIR/free.matrix_with_fuse.txt" \
    --fix_state "$TEST_DIR/input.primitive_state.tsv" \
    --free_state "$TEST_DIR/input.primitive_state.tsv" \
    --species_tree "$TMP_DIR/species_tree.forSplit.nwk" \
    -o "$OUT_PREFIX"
)

grep -q $'^B5\tinternal\t1\t1\t1\t100.0000000000\t0.0000000000$' "$OUT_PREFIX.support_b.txt"
grep -q '100.0000000000' "$SUPPORT_TREE"

count="$(grep -o '100.0000000000' "$SUPPORT_TREE" | wc -l | tr -d ' ')"
[ "$count" = "2" ]

if grep -q ')0\.0000000000' "$SUPPORT_TREE"; then
    echo "[FAIL] Duplicate rooted support label was written as 0 in support_b.nwk" >&2
    exit 1
fi

echo "[PASS] support duplicate root regression test"
