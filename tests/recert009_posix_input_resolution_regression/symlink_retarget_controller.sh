#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
SA="$ROOT/SplitAligner.pl"
TMP="$(mktemp -d "${TMPDIR:-/tmp}/splitaligner-recert009-retarget-XXXXXX")"

cleanup() {
    chmod -R u+w "$TMP" 2>/dev/null || true
    rm -rf "$TMP"
}
trap cleanup EXIT

fail() { printf '[FAIL] %s\n' "$*" >&2; exit 1; }

snapshot_outputs() {
    local dir="$1" output="$2"
    (
        cd "$dir"
        find . -type f ! -name '*.stdout.txt' ! -name '*.stderr.txt' \
            -print0 | sort -z | xargs -0 shasum -a 256
    ) >"$output"
}

wait_for_workdir_and_stop() {
    local pid="$1" dir="$2"
    for _ in $(seq 1 8000); do
        if find "$dir" -maxdepth 1 -type d -name '.splitaligner-*' -print -quit | grep -q .; then
            kill -STOP "$pid" 2>/dev/null || return 1
            return 0
        fi
        kill -0 "$pid" 2>/dev/null || break
        sleep 0.001
    done
    return 1
}

CASE="$TMP/case"
mkdir -p "$CASE/base/first/sub" "$CASE/base/second/sub"
printf '%s\n' '((A:1,B:1):2,(C:1,D:1):3);' >"$CASE/base/first/species.nwk"
printf '%s\n' '((A:1,C:1):2,(B:1,D:1):3);' >"$CASE/base/second/species.nwk"
printf '%s\n' 'g((A:1,B:2):3,(C:4,D:5):6);' >"$CASE/genes.nwk"
ln -s first/sub "$CASE/base/link"
PROVIDED="$CASE/base/link/../species.nwk"

(
    cd "$CASE"
    perl "$SA" --mode matrix --species "$PROVIDED" --gene genes.nwk --label x
) >"$CASE/initial.stdout.txt" 2>"$CASE/initial.stderr.txt"
snapshot_outputs "$CASE" "$TMP/before.sha256"

(
    cd "$CASE"
    exec perl "$SA" --mode matrix --species "$PROVIDED" --gene genes.nwk --label x --force
) >"$CASE/retarget.stdout.txt" 2>"$CASE/retarget.stderr.txt" &
pid=$!
wait_for_workdir_and_stop "$pid" "$CASE" || fail 'could not pause symlink-retarget fixture'
rm "$CASE/base/link"
ln -s second/sub "$CASE/base/link"
kill -CONT "$pid"

set +e
wait "$pid"
rc=$?
set -e
[[ $rc -ne 0 ]] || fail 'retargeted input unexpectedly published'
grep -Fq "Input role 'matrix species tree'" "$CASE/retarget.stderr.txt" ||
    fail 'retarget diagnostic omitted input role'
grep -Fq 'reason: resolved-path-changed' "$CASE/retarget.stderr.txt" ||
    fail 'retarget diagnostic omitted resolved-path reason'

snapshot_outputs "$CASE" "$TMP/after.sha256"
grep -vE './base/link|./retarget\.(stdout|stderr)\.txt' "$TMP/before.sha256" >"$TMP/before.outputs.sha256"
grep -vE './base/link|./retarget\.(stdout|stderr)\.txt' "$TMP/after.sha256" >"$TMP/after.outputs.sha256"
cmp "$TMP/before.outputs.sha256" "$TMP/after.outputs.sha256" ||
    fail 'retargeted input changed previously published outputs'
if find "$CASE" -maxdepth 1 -type d -name '.splitaligner-*' -print -quit | grep -q .; then
    fail 'retarget failure left transaction debris'
fi

echo '[PASS] PATH-009-15 controller symlink retarget rejected before publication'
