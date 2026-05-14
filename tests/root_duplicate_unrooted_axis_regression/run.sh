#!/usr/bin/env bash

set -euo pipefail

TEST_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$TEST_DIR/../.." && pwd)"
WORK_ROOT="$TEST_DIR/tmp_hashseed_runs"

rm -rf "$WORK_ROOT"
mkdir -p "$WORK_ROOT"

for seed in 1 2 3 4 5 6 7 8; do
    RUN_DIR="$WORK_ROOT/seed_$seed"
    mkdir -p "$RUN_DIR"

    (
        cd "$RUN_DIR"
        PERL_HASH_SEED="$seed" PERL_PERTURB_KEYS=1 \
        perl "$REPO_ROOT/SplitAligner.pl" \
            --mode matrix \
            --species "$TEST_DIR/species.nwk" \
            --gene "$TEST_DIR/genes.nwk" \
            --label test
    )
done

REF_DIR="$WORK_ROOT/seed_1"
for seed in 2 3 4 5 6 7 8; do
    RUN_DIR="$WORK_ROOT/seed_$seed"
    cmp "$REF_DIR/species_tree.splits.txt" "$RUN_DIR/species_tree.splits.txt"
    cmp "$REF_DIR/test.matrix_no_fuse.txt" "$RUN_DIR/test.matrix_no_fuse.txt"
    cmp "$REF_DIR/test.matrix_with_fuse.txt" "$RUN_DIR/test.matrix_with_fuse.txt"
done

ROOT_COLLAPSED_LINE="$(grep $'^A||B..C..D\t' "$REF_DIR/species_tree.splits.txt" || true)"
[ -n "$ROOT_COLLAPSED_LINE" ]
case "$ROOT_COLLAPSED_LINE" in
    $'A||B..C..D\tB1') ;;
    *)
        echo "[FAIL] Expected root-collapsed split A||B..C..D to map to B1, got: $ROOT_COLLAPSED_LINE" >&2
        exit 1
        ;;
esac

grep -q $'^B1\t.*duplicate_unrooted_split_winner_over=B6$' "$REF_DIR/species_tree.branch_map.txt"
grep -q $'^B6\t.*duplicate_unrooted_split_loser_of=B1$' "$REF_DIR/species_tree.branch_map.txt"

echo "[PASS] root duplicate unrooted axis regression test"
