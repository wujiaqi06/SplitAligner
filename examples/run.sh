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
# Usage:
#   bash examples/run.sh toy
#   bash examples/run.sh preprint
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
    perl "$REPO_ROOT/SplitAligner.pl" --mode matrix --species speciesTree302.nwk --gene free.2275genes.nwk --label free
    perl "$REPO_ROOT/SplitAligner.pl" --mode matrix --species speciesTree302.nwk --gene fix.2275genes.nwk --label fix
    perl "$REPO_ROOT/SplitAligner.pl" --mode finalize --free free.matrix_with_fuse.txt --fix fix.matrix_with_fuse.txt --final_label final --species_tree species_tree.forSplit.nwk
}

case "$MODE" in
    toy)
        run_toy
        ;;
    preprint)
        run_preprint
        ;;
    *)
        echo "Usage: bash examples/run.sh [toy|preprint]" >&2
        exit 1
        ;;
esac
