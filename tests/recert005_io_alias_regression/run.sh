#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
SA="$ROOT/SplitAligner.pl"
PROPERTY="$ROOT/tests/recert005_io_alias_regression/io_safety_property_test.pl"
REWRITE="$ROOT/tests/recert005_io_alias_regression/rewrite_state_manifest.py"
TMP="$(mktemp -d "${TMPDIR:-/tmp}/splitaligner-recert005-XXXXXX")"

cleanup() {
    chmod -R u+w "$TMP" 2>/dev/null || true
    rm -rf "$TMP"
}
trap cleanup EXIT

sha256() {
    shasum -a 256 "$1" | awk '{print $1}'
}

fail() {
    echo "[FAIL] $*" >&2
    exit 1
}

pass() {
    echo "[PASS] $*"
}

assert_no_workdir() {
    local dir="$1"
    if find "$dir" -maxdepth 1 -name '.splitaligner-*' -print -quit | grep -q .; then
        fail "temporary workspace remained in $dir"
    fi
}

expect_rejection() {
    local name="$1"
    local dir="$2"
    local reason="$3"
    shift 3
    set +e
    (
        cd "$dir"
        "$@"
    ) >"$dir/$name.stdout.txt" 2>"$dir/$name.stderr.txt"
    local rc=$?
    set -e
    [[ $rc -ne 0 ]] || fail "$name unexpectedly succeeded"
    grep -Fq -- "$reason" "$dir/$name.stderr.txt" || {
        cat "$dir/$name.stderr.txt" >&2
        fail "$name did not report '$reason'"
    }
    grep -Fq -- '--force cannot override' "$dir/$name.stderr.txt" ||
        fail "$name did not state that --force cannot override protection"
    assert_no_workdir "$dir"
    pass "$name"
}

copy_seed() {
    local destination="$1"
    mkdir -p "$destination"
    cp -R "$SEED"/. "$destination"/
}

snapshot_outputs() {
    local dir="$1"
    local output="$2"
    (
        cd "$dir"
        find . -type f \
            ! -name '*.stdout.txt' \
            ! -name '*.stderr.txt' \
            -print0 | sort -z | xargs -0 shasum -a 256
    ) >"$output"
}

wait_for_workdir_and_stop() {
    local pid="$1"
    local dir="$2"
    local found=0
    for _ in $(seq 1 5000); do
        if find "$dir" -maxdepth 1 -type d -name '.splitaligner-*' -print -quit | grep -q .; then
            kill -STOP "$pid" 2>/dev/null || return 1
            found=1
            break
        fi
        kill -0 "$pid" 2>/dev/null || break
        sleep 0.001
    done
    [[ $found -eq 1 ]]
}

perl "$PROPERTY"

SEED="$TMP/seed"
mkdir -p "$SEED"
cat >"$SEED/species.nwk" <<'EOF'
((A:1,B:1):2,(C:1,D:1):3);
EOF
cat >"$SEED/genes.nwk" <<'EOF'
g((A:1,B:1):2,(C:1,D:1):3);
EOF
(
    cd "$SEED"
    perl "$SA" --mode matrix --species species.nwk --gene genes.nwk --label free
    perl "$SA" --mode matrix --species species.nwk --gene genes.nwk --label fix
) >"$TMP/seed.stdout.txt" 2>"$TMP/seed.stderr.txt"

# RECERT005-01: exact species alias, with and without --force.
for force_mode in default force; do
    CASE="$TMP/01_species_$force_mode"
    mkdir -p "$CASE"
    cp "$SEED/species.nwk" "$CASE/species_tree.forSplit.nwk"
    cp "$SEED/genes.nwk" "$CASE/genes.nwk"
    before="$(sha256 "$CASE/species_tree.forSplit.nwk")"
    args=(perl "$SA" --mode matrix --species species_tree.forSplit.nwk --gene genes.nwk --label x)
    [[ $force_mode == force ]] && args+=(--force)
    expect_rejection "species_$force_mode" "$CASE" "reason: exact-path" "${args[@]}"
    [[ "$(sha256 "$CASE/species_tree.forSplit.nwk")" == "$before" ]] || fail "species input changed"
    [[ ! -e "$CASE/x.matrix_with_fuse.txt" ]] || fail "species alias published a matrix"
done

# RECERT005-02: exact gene alias.
CASE="$TMP/02_gene"
mkdir -p "$CASE"
cp "$SEED/species.nwk" "$CASE/species.nwk"
cp "$SEED/genes.nwk" "$CASE/x.matrix_with_fuse.txt"
before="$(sha256 "$CASE/x.matrix_with_fuse.txt")"
expect_rejection gene_alias "$CASE" "reason: exact-path" \
    perl "$SA" --mode matrix --species species.nwk --gene x.matrix_with_fuse.txt --label x --force
[[ "$(sha256 "$CASE/x.matrix_with_fuse.txt")" == "$before" ]] || fail "gene input changed"

# RECERT005-03: both matrix-owned output directories protect nested inputs.
for suffix in splits split_branch_label; do
    CASE="$TMP/03_directory_$suffix"
    mkdir -p "$CASE/x_$suffix"
    cp "$SEED/species.nwk" "$CASE/species.nwk"
    cp "$SEED/genes.nwk" "$CASE/x_$suffix/genes.nwk"
    before="$(sha256 "$CASE/x_$suffix/genes.nwk")"
    expect_rejection "directory_$suffix" "$CASE" "reason: inside-output-directory" \
        perl "$SA" --mode matrix --species species.nwk --gene "x_$suffix/genes.nwk" --label x --force
    [[ "$(sha256 "$CASE/x_$suffix/genes.nwk")" == "$before" ]] || fail "nested gene input changed"
done

# RECERT005-04: destination symlink resolves to an input.
CASE="$TMP/04_symlink"
mkdir -p "$CASE"
cp "$SEED/species.nwk" "$CASE/species.nwk"
cp "$SEED/genes.nwk" "$CASE/genes.nwk"
ln -s species.nwk "$CASE/species_tree.forSplit.nwk"
before="$(sha256 "$CASE/species.nwk")"
expect_rejection symlink_alias "$CASE" "reason: resolved-path" \
    perl "$SA" --mode matrix --species species.nwk --gene genes.nwk --label x --force
[[ -L "$CASE/species_tree.forSplit.nwk" ]] || fail "destination symlink was replaced"
[[ "$(sha256 "$CASE/species.nwk")" == "$before" ]] || fail "symlink target input changed"

# RECERT005-05: destination hardlink shares an input inode.
CASE="$TMP/05_hardlink"
mkdir -p "$CASE"
cp "$SEED/species.nwk" "$CASE/species.nwk"
cp "$SEED/genes.nwk" "$CASE/genes.nwk"
ln "$CASE/species.nwk" "$CASE/species_tree.forSplit.nwk"
before="$(sha256 "$CASE/species.nwk")"
expect_rejection hardlink_alias "$CASE" "reason: same-inode" \
    perl "$SA" --mode matrix --species species.nwk --gene genes.nwk --label x --force
[[ "$(sha256 "$CASE/species_tree.forSplit.nwk")" == "$before" ]] || fail "hardlink input changed"

# RECERT005-06: FREE matrix equals the final FREE classified destination.
CASE="$TMP/06_finalize_free"
copy_seed "$CASE"
cp "$CASE/free.matrix_with_fuse.txt" "$CASE/out.free.na_classified.txt"
before="$(sha256 "$CASE/out.free.na_classified.txt")"
expect_rejection finalize_free_alias "$CASE" "reason: exact-path" \
    perl "$SA" --mode finalize \
    --free out.free.na_classified.txt --free_manifest free.run_manifest.json \
    --fix fix.matrix_with_fuse.txt --fix_manifest fix.run_manifest.json \
    --final_label out --force
[[ "$(sha256 "$CASE/out.free.na_classified.txt")" == "$before" ]] || fail "FREE input changed"

# RECERT005-07: a derived FIX .na_fuse destination equals the FREE input.
CASE="$TMP/07_cross_role"
copy_seed "$CASE"
cp "$CASE/free.matrix_with_fuse.txt" "$CASE/fix.matrix_with_fuse.na_fuse.txt"
before="$(sha256 "$CASE/fix.matrix_with_fuse.na_fuse.txt")"
expect_rejection finalize_cross_role "$CASE" "reason: exact-path" \
    perl "$SA" --mode finalize \
    --free fix.matrix_with_fuse.na_fuse.txt --free_manifest free.run_manifest.json \
    --fix fix.matrix_with_fuse.txt --fix_manifest fix.run_manifest.json \
    --final_label out --force
[[ "$(sha256 "$CASE/fix.matrix_with_fuse.na_fuse.txt")" == "$before" ]] || fail "cross-role input changed"

# RECERT005-08: explicit and inferred manifests are protected.
CASE="$TMP/08_manifest_explicit"
copy_seed "$CASE"
cp "$CASE/free.run_manifest.json" "$CASE/out.finalize_manifest.json"
before="$(sha256 "$CASE/out.finalize_manifest.json")"
expect_rejection explicit_manifest_alias "$CASE" "reason: exact-path" \
    perl "$SA" --mode finalize \
    --free free.matrix_with_fuse.txt --free_manifest out.finalize_manifest.json \
    --fix fix.matrix_with_fuse.txt --fix_manifest fix.run_manifest.json \
    --final_label out --force
[[ "$(sha256 "$CASE/out.finalize_manifest.json")" == "$before" ]] || fail "explicit manifest changed"

CASE="$TMP/08_manifest_inferred"
copy_seed "$CASE"
ln "$CASE/free.run_manifest.json" "$CASE/out.finalize_manifest.json"
before="$(sha256 "$CASE/free.run_manifest.json")"
expect_rejection inferred_manifest_alias "$CASE" "reason: same-inode" \
    perl "$SA" --mode finalize \
    --free free.matrix_with_fuse.txt --fix fix.matrix_with_fuse.txt \
    --final_label out --force
[[ "$(sha256 "$CASE/free.run_manifest.json")" == "$before" ]] || fail "inferred manifest changed"

# RECERT005-09: manifest-discovered coordinate-state sidecar is protected.
CASE="$TMP/09_state"
copy_seed "$CASE"
cp "$CASE/free.primitive_state.tsv" "$CASE/out.free.na_classified.txt"
python3 "$REWRITE" \
    "$CASE/free.run_manifest.json" "$CASE/out.free.na_classified.txt" \
    out.free.na_classified.txt "$CASE/free.alias.run_manifest.json"
before="$(sha256 "$CASE/out.free.na_classified.txt")"
expect_rejection state_alias "$CASE" "reason: exact-path" \
    perl "$SA" --mode finalize \
    --free free.matrix_with_fuse.txt --free_manifest free.alias.run_manifest.json \
    --fix fix.matrix_with_fuse.txt --fix_manifest fix.run_manifest.json \
    --final_label out --force
[[ "$(sha256 "$CASE/out.free.na_classified.txt")" == "$before" ]] || fail "state sidecar changed"

# RECERT005-10: Support tree and derived branch-map inputs are protected.
CASE="$TMP/10_support_tree"
copy_seed "$CASE"
cp "$CASE/species_tree.forSplit.nwk" "$CASE/out.support_b.txt"
before="$(sha256 "$CASE/out.support_b.txt")"
expect_rejection support_tree_alias "$CASE" "reason: exact-path" \
    perl "$SA" --mode finalize \
    --free free.matrix_with_fuse.txt --fix fix.matrix_with_fuse.txt \
    --final_label out --species_tree out.support_b.txt --force
[[ "$(sha256 "$CASE/out.support_b.txt")" == "$before" ]] || fail "Support tree input changed"

CASE="$TMP/10_support_map"
copy_seed "$CASE"
ln "$CASE/species_tree.branch_map.txt" "$CASE/out.support_b.txt"
before="$(sha256 "$CASE/species_tree.branch_map.txt")"
expect_rejection support_map_alias "$CASE" "reason: same-inode" \
    perl "$SA" --mode finalize \
    --free free.matrix_with_fuse.txt --fix fix.matrix_with_fuse.txt \
    --final_label out --species_tree species_tree.forSplit.nwk --force
[[ "$(sha256 "$CASE/species_tree.branch_map.txt")" == "$before" ]] || fail "Support branch-map input changed"

# RECERT005-11: finalize_fix input aliases are protected.
CASE="$TMP/11_finalize_fix"
copy_seed "$CASE"
cp "$CASE/fix.matrix_with_fuse.txt" "$CASE/out.fix.na_classified.txt"
before="$(sha256 "$CASE/out.fix.na_classified.txt")"
expect_rejection finalize_fix_alias "$CASE" "reason: exact-path" \
    perl "$SA" --mode finalize_fix \
    --fix out.fix.na_classified.txt --fix_manifest fix.run_manifest.json \
    --final_label out --force
[[ "$(sha256 "$CASE/out.fix.na_classified.txt")" == "$before" ]] || fail "finalize_fix input changed"

# RECERT005-12: publication destinations may not share one inode.
CASE="$TMP/12_duplicate_destinations"
copy_seed "$CASE"
printf 'prior\n' >"$CASE/free.matrix_with_fuse.na_fuse.txt"
ln "$CASE/free.matrix_with_fuse.na_fuse.txt" "$CASE/fix.matrix_with_fuse.na_fuse.txt"
expect_rejection duplicate_destinations "$CASE" "reason: duplicate-destination" \
    perl "$SA" --mode finalize \
    --free free.matrix_with_fuse.txt --fix fix.matrix_with_fuse.txt \
    --final_label out --force

# RECERT005-13: Unicode CWD, label, absolute path, and decoded diagnostic.
CASE="$TMP/13_ユニコード_路径"
mkdir -p "$CASE/标签_splits"
cp "$SEED/species.nwk" "$CASE/物种.nwk"
cp "$SEED/genes.nwk" "$CASE/标签_splits/基因.nwk"
before="$(sha256 "$CASE/标签_splits/基因.nwk")"
expect_rejection unicode_directory_alias "$CASE" "reason: inside-output-directory" \
    perl "$SA" --mode matrix \
    --species "$CASE/物种.nwk" --gene "$CASE/标签_splits/基因.nwk" \
    --label 标签 --force
grep -Fq '标签_splits/基因.nwk' "$CASE/unicode_directory_alias.stderr.txt" ||
    fail "Unicode diagnostic was not decoded exactly"
[[ "$(sha256 "$CASE/标签_splits/基因.nwk")" == "$before" ]] || fail "Unicode input changed"

# RECERT005-14: benign --force replaces outputs while preserving inputs.
CASE="$TMP/14_benign_force"
mkdir -p "$CASE"
cp "$SEED/species.nwk" "$CASE/species.nwk"
cp "$SEED/genes.nwk" "$CASE/genes.nwk"
(
    cd "$CASE"
    perl "$SA" --mode matrix --species species.nwk --gene genes.nwk --label x
) >"$CASE/first.stdout.txt" 2>"$CASE/first.stderr.txt"
first_matrix="$(sha256 "$CASE/x.matrix_with_fuse.txt")"
cat >"$CASE/genes.nwk" <<'EOF'
g((A:0.1,B:0.2):0.4,(C:0.3,D:0.5):0.7);
EOF
species_before="$(sha256 "$CASE/species.nwk")"
gene_before="$(sha256 "$CASE/genes.nwk")"
(
    cd "$CASE"
    perl "$SA" --mode matrix --species species.nwk --gene genes.nwk --label x --force
) >"$CASE/force.stdout.txt" 2>"$CASE/force.stderr.txt"
[[ "$(sha256 "$CASE/species.nwk")" == "$species_before" ]] || fail "benign force changed species input"
[[ "$(sha256 "$CASE/genes.nwk")" == "$gene_before" ]] || fail "benign force changed gene input"
[[ "$(sha256 "$CASE/x.matrix_with_fuse.txt")" != "$first_matrix" ]] || fail "benign force did not replace outputs"
pass "benign --force replacement"

# RECERT005-15: a publication-time permission failure rolls back all prior outputs.
CASE="$TMP/15_publish_rollback"
mkdir -p "$CASE/free_input" "$CASE/fix_input"
cp "$SEED/free.matrix_with_fuse.txt" "$CASE/free_input/free.matrix_with_fuse.txt"
cp "$SEED/free.run_manifest.json" "$CASE/free_input/free.run_manifest.json"
cp "$SEED/free.primitive_state.tsv" "$CASE/free_input/free.primitive_state.tsv"
cp "$SEED/fix.matrix_with_fuse.txt" "$CASE/fix_input/fix.matrix_with_fuse.txt"
cp "$SEED/fix.run_manifest.json" "$CASE/fix_input/fix.run_manifest.json"
cp "$SEED/fix.primitive_state.tsv" "$CASE/fix_input/fix.primitive_state.tsv"
(
    cd "$CASE"
    perl "$SA" --mode finalize \
        --free free_input/free.matrix_with_fuse.txt \
        --fix fix_input/fix.matrix_with_fuse.txt --final_label out
) >"$CASE/initial.stdout.txt" 2>"$CASE/initial.stderr.txt"
# Preserve scientific content while changing one input-manifest byte so that
# the new finalization manifest is not byte-identical to the prior output.
printf ' ' >>"$CASE/free_input/free.run_manifest.json"
snapshot_outputs "$CASE" "$TMP/15_before.sha256"
set +e
(
    cd "$CASE"
    SPLITALIGNER_INTERNAL_TESTING=1 \
    SPLITALIGNER_TEST_FAIL_AFTER_PUBLISH=1 \
    perl "$SA" --mode finalize \
        --free free_input/free.matrix_with_fuse.txt \
        --fix fix_input/fix.matrix_with_fuse.txt --final_label out --force
) >"$CASE/failed.stdout.txt" 2>"$CASE/failed.stderr.txt"
rc=$?
set -e
[[ $rc -ne 0 ]] || fail "publication rollback fixture unexpectedly succeeded"
snapshot_outputs "$CASE" "$TMP/15_after.sha256"
cmp "$TMP/15_before.sha256" "$TMP/15_after.sha256" || fail "failed publication changed prior outputs or inputs"
assert_no_workdir "$CASE"
pass "publication failure rollback"

# RECERT005-16: mutation after preflight aborts before publication and preserves old outputs.
CASE="$TMP/16_input_mutation"
mkdir -p "$CASE"
cp "$SEED/species.nwk" "$CASE/species.nwk"
cp "$SEED/genes.nwk" "$CASE/genes.nwk"
(
    cd "$CASE"
    perl "$SA" --mode matrix --species species.nwk --gene genes.nwk --label x
) >"$CASE/initial.stdout.txt" 2>"$CASE/initial.stderr.txt"
snapshot_outputs "$CASE" "$TMP/16_before.sha256"
(
    cd "$CASE"
    exec perl "$SA" --mode matrix --species species.nwk --gene genes.nwk --label x --force
) >"$CASE/mutated.stdout.txt" 2>"$CASE/mutated.stderr.txt" &
pid=$!
wait_for_workdir_and_stop "$pid" "$CASE" || fail "could not pause input-mutation fixture"
printf '\n' >>"$CASE/genes.nwk"
kill -CONT "$pid"
set +e
wait "$pid"
rc=$?
set -e
[[ $rc -ne 0 ]] || fail "input mutation fixture unexpectedly succeeded"
grep -Fq "Input role 'matrix gene trees'" "$CASE/mutated.stderr.txt" || fail "mutation diagnostic omitted input role"
grep -Fq 'reason: content-sha256' "$CASE/mutated.stderr.txt" || fail "mutation diagnostic omitted SHA reason"
snapshot_outputs "$CASE" "$TMP/16_after.sha256"
grep -v './genes.nwk' "$TMP/16_before.sha256" >"$TMP/16_before_outputs.sha256"
grep -v './genes.nwk' "$TMP/16_after.sha256" >"$TMP/16_after_outputs.sha256"
cmp "$TMP/16_before_outputs.sha256" "$TMP/16_after_outputs.sha256" || fail "input mutation published new outputs"
assert_no_workdir "$CASE"
pass "pre-publication input mutation rejection"

echo "[PASS] RECERT-005 I/O alias regression suite"
