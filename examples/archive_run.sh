#!/usr/bin/env bash

# Non-destructively archive a bundled example run workspace.
#
# Usage:
#   bash examples/archive_run.sh toy
#   bash examples/archive_run.sh preprint

set -euo pipefail

MODE="${1:-}"
SCRIPT_DIR="$(cd -P -- "$(dirname "$0")" && pwd -P)"
REPO_ROOT="$(cd -P -- "$SCRIPT_DIR/.." && pwd -P)"
LOCK_OWNED=0
LOCK_ENTRY=""

archive_error() {
    printf '[ERROR][Example archive] %s\n' "$1" >&2
    return 1
}

cleanup_lock() {
    if [[ "$LOCK_OWNED" == "1" && -n "$LOCK_ENTRY" ]]; then
        rmdir -- "$LOCK_ENTRY" 2>/dev/null || true
    fi
}

path_signature() {
    perl -e '
        use strict;
        use warnings;
        my ($path) = @ARGV;
        my @st = lstat($path);
        die "cannot lstat $path: $!\n" unless @st;
        die "path is a symbolic link: $path\n" if -l _;
        die "path is not a directory: $path\n" unless -d _;
        print "$st[0]:$st[1]\n";
    ' "$1"
}

archive_example_workspace() {
    local requested_example_root="$1"
    local physical_example_root
    local run_entry
    local physical_run
    local current_physical
    local initial_signature
    local current_signature
    local parent_signature
    local timestamp
    local archive_base
    local archive_name
    local archive_entry
    local physical_archive
    local archived_signature
    local counter=0

    if ! physical_example_root="$(cd -P -- "$requested_example_root" 2>/dev/null && pwd -P)"; then
        archive_error "cannot resolve example root: $requested_example_root"
        return 1
    fi

    run_entry="$physical_example_root/run"

    if [[ -L "$run_entry" ]]; then
        archive_error "run entry must not be a symbolic link: $run_entry"
        return 1
    fi
    if [[ ! -e "$run_entry" ]]; then
        printf '[INFO][Example archive] no run workspace exists: %s\n' "$run_entry"
        return 0
    fi
    if [[ ! -d "$run_entry" ]]; then
        archive_error "run entry exists but is not a directory: $run_entry"
        return 1
    fi
    if ! physical_run="$(cd -P -- "$run_entry" 2>/dev/null && pwd -P)"; then
        archive_error "cannot resolve run directory: $run_entry"
        return 1
    fi
    if [[ "$physical_run" != "$physical_example_root/run" ]]; then
        archive_error "run directory is not the direct physical child of the example root: $run_entry"
        return 1
    fi
    if ! initial_signature="$(path_signature "$run_entry" 2>/dev/null)"; then
        archive_error "cannot record run directory identity: $run_entry"
        return 1
    fi
    if ! parent_signature="$(path_signature "$physical_example_root" 2>/dev/null)"; then
        archive_error "cannot record example root identity: $physical_example_root"
        return 1
    fi
    if [[ "${initial_signature%%:*}" != "${parent_signature%%:*}" ]]; then
        archive_error "run directory is on a different device or is a mounted workspace; archive it manually without recursive deletion: $run_entry"
        return 1
    fi

    LOCK_ENTRY="$physical_example_root/.splitaligner-archive-run.lock"
    if ! mkdir -- "$LOCK_ENTRY" 2>/dev/null; then
        archive_error "another archive operation may be active, or the example root is not writable: $LOCK_ENTRY"
        return 1
    fi
    LOCK_OWNED=1
    trap cleanup_lock EXIT
    trap 'cleanup_lock; exit 1' HUP INT TERM

    if [[ -L "$run_entry" || ! -d "$run_entry" ]]; then
        archive_error "run entry changed before archive: $run_entry"
        return 1
    fi
    if ! cd -P -- "$run_entry"; then
        archive_error "cannot enter run directory: $run_entry"
        return 1
    fi
    if ! current_physical="$(pwd -P)" || [[ "$current_physical" != "$physical_run" ]]; then
        archive_error "run directory identity changed before archive: $run_entry"
        return 1
    fi
    if ! current_signature="$(path_signature . 2>/dev/null)" || [[ "$current_signature" != "$initial_signature" ]]; then
        archive_error "run directory object changed before archive: $run_entry"
        return 1
    fi
    if ! cd -P -- "$physical_example_root"; then
        archive_error "cannot return to physical example root: $physical_example_root"
        return 1
    fi
    if [[ "$(pwd -P)" != "$physical_example_root" ]]; then
        archive_error "example root identity changed before archive: $physical_example_root"
        return 1
    fi

    timestamp="$(date -u '+%Y%m%dT%H%M%SZ')"
    archive_base="run.saved.${timestamp}.$$"
    archive_name="$archive_base"
    archive_entry="$physical_example_root/$archive_name"
    while [[ -e "$archive_entry" || -L "$archive_entry" ]]; do
        counter=$((counter + 1))
        archive_name="${archive_base}.${counter}"
        archive_entry="$physical_example_root/$archive_name"
    done

    if [[ -L "$run_entry" || ! -d "$run_entry" ]]; then
        archive_error "run entry changed immediately before archive: $run_entry"
        return 1
    fi
    if ! current_signature="$(path_signature "$run_entry" 2>/dev/null)" || [[ "$current_signature" != "$initial_signature" ]]; then
        archive_error "run directory object changed immediately before archive: $run_entry"
        return 1
    fi

    # Use the rename(2) syscall directly. Unlike a general-purpose move command,
    # this never falls back to copy-and-delete across filesystems.
    if ! perl -e '
        use strict;
        use warnings;
        my ($source, $destination) = @ARGV;
        die "archive destination already exists: $destination\n" if lstat($destination);
        rename($source, $destination)
            or die "atomic archive rename failed for $source: $!\n";
    ' "$run_entry" "$archive_entry"; then
        archive_error "atomic archive rename failed; the original workspace was not intentionally deleted: $run_entry"
        return 1
    fi

    if [[ -e "$run_entry" || -L "$run_entry" ]]; then
        archive_error "source run entry still exists after archive: $run_entry"
        return 1
    fi
    if [[ -L "$archive_entry" || ! -d "$archive_entry" ]]; then
        archive_error "archive destination is not the expected directory: $archive_entry"
        return 1
    fi
    if ! physical_archive="$(cd -P -- "$archive_entry" 2>/dev/null && pwd -P)" || [[ "$physical_archive" != "$archive_entry" ]]; then
        archive_error "cannot verify physical archive path: $archive_entry"
        return 1
    fi
    if ! archived_signature="$(path_signature "$archive_entry" 2>/dev/null)" || [[ "$archived_signature" != "$initial_signature" ]]; then
        archive_error "archived directory identity does not match the original workspace: $archive_entry"
        return 1
    fi

    printf '[INFO][Example archive] archived %s workspace: %s\n' "$MODE" "$archive_entry"
    printf '[INFO][Example archive] inspect this archive before any manual deletion; no workspace contents were deleted.\n'
}

case "$MODE" in
    toy)
        archive_example_workspace "$REPO_ROOT/examples/302mammal"
        ;;
    preprint)
        archive_example_workspace "$REPO_ROOT/examples/preprint_302mammal"
        ;;
    *)
        echo "Usage: bash examples/archive_run.sh [toy|preprint]" >&2
        exit 1
        ;;
esac
