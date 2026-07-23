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
# - Immutable inputs and expected outputs remain outside the generated run/ dir.
# - Each rerun uses ownership-aware --force inside its dedicated run/ workspace.

set -euo pipefail

MODE="${1:-toy}"
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
export PATH="$PATH:$REPO_ROOT"

example_workspace_error() {
    printf '[ERROR][Example workspace] %s\n' "$1" >&2
    return 1
}

enter_example_workspace() {
    local requested_example_root="$1"
    local physical_example_root
    local run_entry
    local physical_run
    local physical_fixture
    local current_physical
    local fixture_name

    if ! physical_example_root="$(cd -P -- "$requested_example_root" 2>/dev/null && pwd -P)"; then
        example_workspace_error "cannot resolve example root: $requested_example_root"
        return 1
    fi

    run_entry="$physical_example_root/run"

    if [[ -L "$run_entry" ]]; then
        example_workspace_error "run entry must not be a symbolic link: $run_entry"
        return 1
    fi
    if [[ -e "$run_entry" && ! -d "$run_entry" ]]; then
        example_workspace_error "run entry exists but is not a directory: $run_entry"
        return 1
    fi
    if [[ ! -e "$run_entry" ]]; then
        if ! mkdir -- "$run_entry"; then
            example_workspace_error "cannot create run directory: $run_entry"
            return 1
        fi
    fi

    if [[ -L "$run_entry" ]]; then
        example_workspace_error "run entry became a symbolic link: $run_entry"
        return 1
    fi
    if [[ ! -d "$run_entry" ]]; then
        example_workspace_error "run entry is not a directory after creation: $run_entry"
        return 1
    fi
    if ! physical_run="$(cd -P -- "$run_entry" 2>/dev/null && pwd -P)"; then
        example_workspace_error "cannot resolve run directory: $run_entry"
        return 1
    fi
    if [[ "$physical_run" != "$physical_example_root/run" ]]; then
        example_workspace_error "run directory is not the direct physical child of the example root: $run_entry"
        return 1
    fi

    for fixture_name in input expected; do
        if [[ -e "$physical_example_root/$fixture_name" ]]; then
            if ! physical_fixture="$(cd -P -- "$physical_example_root/$fixture_name" 2>/dev/null && pwd -P)"; then
                example_workspace_error "cannot resolve immutable fixture directory: $physical_example_root/$fixture_name"
                return 1
            fi
            if [[ "$physical_run" == "$physical_fixture" || "$run_entry" -ef "$physical_example_root/$fixture_name" ]]; then
                example_workspace_error "run directory overlaps immutable $fixture_name directory: $run_entry"
                return 1
            fi
        fi
    done

    if ! cd -P -- "$run_entry"; then
        example_workspace_error "cannot enter run directory: $run_entry"
        return 1
    fi
    if ! current_physical="$(pwd -P)" || [[ "$current_physical" != "$physical_run" ]]; then
        example_workspace_error "physical working directory changed before execution: $run_entry"
        return 1
    fi
}

run_toy() {
    local example_root="$REPO_ROOT/examples/302mammal"
    enter_example_workspace "$example_root"
    perl "$REPO_ROOT/SplitAligner.pl" --mode matrix --species ../input/speciesTree302.nwk --gene ../input/free_tree.examples.nwk --label free --force
    perl "$REPO_ROOT/SplitAligner.pl" --mode matrix --species ../input/speciesTree302.nwk --gene ../input/fix_tree.examples.nwk --label fix --force
    perl "$REPO_ROOT/SplitAligner.pl" --mode finalize --free free.matrix_with_fuse.txt --fix fix.matrix_with_fuse.txt --final_label final --species_tree species_tree.forSplit.nwk --force
}

run_preprint() {
    local example_root="$REPO_ROOT/examples/preprint_302mammal"
    enter_example_workspace "$example_root"
    perl "$REPO_ROOT/SplitAligner.pl" --mode matrix --species ../input/speciesTree302.nwk --gene ../input/free.2275genes.nwk --label free --force
    perl "$REPO_ROOT/SplitAligner.pl" --mode matrix --species ../input/speciesTree302.nwk --gene ../input/fix.2275genes.nwk --label fix --force
    perl "$REPO_ROOT/SplitAligner.pl" --mode finalize --free free.matrix_with_fuse.txt --fix fix.matrix_with_fuse.txt --final_label final --species_tree species_tree.forSplit.nwk --force
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
