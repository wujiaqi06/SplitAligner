#!/usr/bin/env bash

# Example guide for SplitAligner
#
# Example 1: examples/302mammal
# - Small toy example for quick smoke testing.
# - Includes a subset of free/fix gene trees and can be run quickly.
# - Good for checking that the workflow and outputs are generated as expected.
#
# Example 2: examples/preprint_302mammal
# - Full dataset used for the preprint-scale 302-mammal analysis.
# - Includes 2275 free-topology gene trees and 2275 fixed-topology gene trees.
# - Intended for reproducing the main analysis outputs, including Support(b).
# 
# Example 3: examples/benchmark
# - Packaged benchmark alignment example derived from the benchmark bundle.
# - Uses a fixed-topology pruning series and is intended to test SplitAligner
#   behavior on benchmark trees.
# - Demonstrates fix-only finalization via `finalize_fix`.
#
# Usage:
#   bash examples/run.sh toy
#   bash examples/run.sh preprint
#   bash examples/run.sh benchmark
#
# Notes:
# - Run this script from the repository root.
# - SplitAligner.pl must be available on PATH, or this script will use the copy
#   in the repository root.

set -euo pipefail

MODE="${1:-toy}"
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export PATH="$PATH:$REPO_ROOT"

run_toy() {
    cd "$REPO_ROOT/examples/302mammal"
    perl "$REPO_ROOT/SplitAligner.pl" --mode matrix --species input/speciesTree302.nwk --gene input/free_tree.examples.nwk --label free
    perl "$REPO_ROOT/SplitAligner.pl" --mode matrix --species input/speciesTree302.nwk --gene input/fix_tree.examples.nwk --label fix
    perl "$REPO_ROOT/SplitAligner.pl" --mode finalize --free free.matrix_with_fuse.txt --fix fix.matrix_with_fuse.txt --final_label final --species_tree species_tree.forSplit.nwk
}

run_preprint() {
    cd "$REPO_ROOT/examples/preprint_302mammal"
    perl "$REPO_ROOT/SplitAligner.pl" --mode matrix --species input/speciesTree302.nwk --gene input/free.2275genes.nwk --label free
    perl "$REPO_ROOT/SplitAligner.pl" --mode matrix --species input/speciesTree302.nwk --gene input/fix.2275genes.nwk --label fix
    perl "$REPO_ROOT/SplitAligner.pl" --mode finalize --free free.matrix_with_fuse.txt --fix fix.matrix_with_fuse.txt --final_label final --species_tree species_tree.forSplit.nwk
}

run_benchmark() {
    cd "$REPO_ROOT/examples/benchmark"
    perl "$REPO_ROOT/SplitAligner.pl" --mode matrix --species input/benchmark.species_tree.nwk --gene input/benchmark.gene_trees.nwk --label benchmark
    perl "$REPO_ROOT/SplitAligner.pl" --mode finalize_fix --fix benchmark.matrix_with_fuse.txt --final_label benchmark_fix
}

case "$MODE" in
    toy)
        run_toy
        ;;
    preprint)
        run_preprint
        ;;
    benchmark)
        run_benchmark
        ;;
    *)
        echo "Usage: bash examples/run.sh [toy|preprint|benchmark]" >&2
        exit 1
        ;;
esac
