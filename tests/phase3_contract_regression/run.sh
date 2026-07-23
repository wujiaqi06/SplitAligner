#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
SA="$ROOT/SplitAligner.pl"
if [[ -n "${SPLITALIGNER_TEST_TMP:-}" ]]; then
    TMP="$SPLITALIGNER_TEST_TMP"
    mkdir -p "$TMP"
    trap 'printf "[INFO] Preserved Phase 3 evidence: %s\n" "$TMP"' EXIT
else
    TMP="$(mktemp -d "${TMPDIR:-/tmp}/splitaligner_phase3_contract_XXXXXX")"
    trap 'rm -rf "$TMP"' EXIT
fi

fail() { printf '[FAIL] %s\n' "$1" >&2; exit 1; }
pass() { printf '[PASS] %s\n' "$1"; }
expect_failure() {
    local label="$1"
    shift
    if "$@" >"$TMP/${label}.stdout" 2>"$TMP/${label}.stderr"; then
        fail "$label unexpectedly exited 0"
    fi
}

mkdir -p "$TMP/fixtures"
cat > "$TMP/fixtures/species.nwk" <<'EOF'
((A,B),(C,D));
EOF

# SA-001: the standalone discordant quartet must retain the complete B1..B5 axis.
cat > "$TMP/fixtures/free_single.nwk" <<'EOF'
g_discordant((A:1,C:1):2,(B:1,D:1):3);
EOF
cat > "$TMP/fixtures/free_batch.nwk" <<'EOF'
g_discordant((A:1,C:1):2,(B:1,D:1):3);
g_support((A:1,B:1):2,(C:1,D:1):3);
EOF
cat > "$TMP/fixtures/fix.nwk" <<'EOF'
g_discordant((A:1,B:1):2,(C:1,D:1):3);
EOF
for run in single batch; do
    mkdir -p "$TMP/sa001/$run"
    gene="$TMP/fixtures/free_single.nwk"
    [[ "$run" == batch ]] && gene="$TMP/fixtures/free_batch.nwk"
    (
        cd "$TMP/sa001/$run"
        perl "$SA" --mode matrix --species "$TMP/fixtures/species.nwk" --gene "$gene" --label free
        perl "$SA" --mode matrix --species "$TMP/fixtures/species.nwk" --gene "$TMP/fixtures/fix.nwk" --label fix
        perl "$SA" --mode finalize --free free.matrix_with_fuse.txt --fix fix.matrix_with_fuse.txt --final_label final
    ) >/dev/null 2>&1
    grep -qx $'gene\tB1\tB2\tB3\tB4\tB5' "$TMP/sa001/$run/free.matrix_no_fuse.txt"
done
grep '^g_discordant' "$TMP/sa001/single/free.matrix_no_fuse.txt" > "$TMP/sa001/single.row"
grep '^g_discordant' "$TMP/sa001/batch/free.matrix_no_fuse.txt" > "$TMP/sa001/batch.row"
cmp -s "$TMP/sa001/single.row" "$TMP/sa001/batch.row"
grep -q $'^g_discordant\t1\t1\t1\t1\tNA_topo$' "$TMP/sa001/single/final.free.na_classified.txt"
pass 'SA-001 exact primitive axis and batch invariance'

# SA-002: conflict by default; --force replaces only after a successful clean run.
cat > "$TMP/fixtures/first.nwk" <<'EOF'
g1((A:1,B:1):2,(C:1,D:1):3);
g2((A:4,C:4):5,(B:4,D:4):6);
EOF
cat > "$TMP/fixtures/second.nwk" <<'EOF'
g1((A:7,B:7):8,(C:7,D:7):9);
EOF
mkdir -p "$TMP/sa002"
(
    cd "$TMP/sa002"
    perl "$SA" --mode matrix --species "$TMP/fixtures/species.nwk" --gene "$TMP/fixtures/first.nwk" --label x
) >/dev/null 2>&1
cp "$TMP/sa002/x.matrix_with_fuse.txt" "$TMP/sa002/before.txt"
expect_failure sa002_conflict bash -c "cd '$TMP/sa002' && perl '$SA' --mode matrix --species '$TMP/fixtures/species.nwk' --gene '$TMP/fixtures/second.nwk' --label x"
cmp -s "$TMP/sa002/before.txt" "$TMP/sa002/x.matrix_with_fuse.txt"
cat > "$TMP/fixtures/invalid_force.nwk" <<'EOF'
g1((A:1,B:1):2,(C:1,X:1):3);
EOF
expect_failure sa002_failed_force bash -c "cd '$TMP/sa002' && perl '$SA' --mode matrix --species '$TMP/fixtures/species.nwk' --gene '$TMP/fixtures/invalid_force.nwk' --label x --force"
cmp -s "$TMP/sa002/before.txt" "$TMP/sa002/x.matrix_with_fuse.txt"
(
    cd "$TMP/sa002"
    perl "$SA" --mode matrix --species "$TMP/fixtures/species.nwk" --gene "$TMP/fixtures/second.nwk" --label x --force
) >/dev/null 2>&1
grep -q '^g1' "$TMP/sa002/x.matrix_with_fuse.txt"
! grep -q '^g2' "$TMP/sa002/x.matrix_with_fuse.txt"
[[ "$(find "$TMP/sa002/x_splits" -name '*.split.txt' | wc -l | tr -d ' ')" == 1 ]]
pass 'SA-002 transactional same-label replacement'

# SA-003: duplicate exact IDs fail; formerly colliding filenames remain distinct.
cat > "$TMP/fixtures/duplicate_ids.nwk" <<'EOF'
g((A:1,B:1):2,(C:1,D:1):3);
g((A:4,C:4):5,(B:4,D:4):6);
EOF
mkdir -p "$TMP/sa003_duplicate"
expect_failure sa003_duplicate bash -c "cd '$TMP/sa003_duplicate' && perl '$SA' --mode matrix --species '$TMP/fixtures/species.nwk' --gene '$TMP/fixtures/duplicate_ids.nwk' --label x"
[[ ! -e "$TMP/sa003_duplicate/x.matrix_with_fuse.txt" ]]
cat > "$TMP/fixtures/collision_ids.nwk" <<'EOF'
g/a((A:1,B:1):2,(C:1,D:1):3);
g_a((A:4,C:4):5,(B:4,D:4):6);
EOF
mkdir -p "$TMP/sa003_collision"
(
    cd "$TMP/sa003_collision"
    perl "$SA" --mode matrix --species "$TMP/fixtures/species.nwk" --gene "$TMP/fixtures/collision_ids.nwk" --label x
) >/dev/null 2>&1
grep -q '^g/a' "$TMP/sa003_collision/x.matrix_with_fuse.txt"
grep -q '^g_a' "$TMP/sa003_collision/x.matrix_with_fuse.txt"
[[ "$(find "$TMP/sa003_collision/x_splits" -name '*.split.txt' | wc -l | tr -d ' ')" == 2 ]]
pass 'SA-003 gene identity uniqueness and injective storage'

# SA-004/005: one numeric grammar and numeric-only terminal labels.
cat > "$TMP/fixtures/numeric_forms.nwk" <<'EOF'
g((A:1.,B:+1e0):-2.5E-1,(C:.5,D:0.5):3.);
EOF
mkdir -p "$TMP/sa004"
(
    cd "$TMP/sa004"
    perl "$SA" --mode matrix --species "$TMP/fixtures/species.nwk" --gene "$TMP/fixtures/numeric_forms.nwk" --label x
) >/dev/null 2>&1
grep -qx $'gene\tB1\tB2\tB3\tB4\tB5' "$TMP/sa004/x.matrix_no_fuse.txt"
grep -q $'^g\t1\.\t+1e0\t\.5\t0\.5\t2\.75$' "$TMP/sa004/x.matrix_no_fuse.txt"
cat > "$TMP/fixtures/numeric_species.nwk" <<'EOF'
((1,2),(3,4));
EOF
cat > "$TMP/fixtures/numeric_taxa.nwk" <<'EOF'
g((1:1,2:2):5,(3:3,4:4):6);
EOF
mkdir -p "$TMP/sa005"
(
    cd "$TMP/sa005"
    perl "$SA" --mode matrix --species "$TMP/fixtures/numeric_species.nwk" --gene "$TMP/fixtures/numeric_taxa.nwk" --label x
) >/dev/null 2>&1
grep -q $'^g\t1\t2\t3\t4\t11$' "$TMP/sa005/x.matrix_no_fuse.txt"
pass 'SA-004/005 shared numeric grammar and numeric taxa'

# SA-006/007: bad taxa or grammar fail before label-owned outputs exist.
cat > "$TMP/fixtures/unknown.nwk" <<'EOF'
g((A:1,B:1):2,(C:1,X:1):3);
EOF
cat > "$TMP/fixtures/duplicate_taxon.nwk" <<'EOF'
g((A:1,A:2):3,(C:1,D:1):4);
EOF
cat > "$TMP/fixtures/unbalanced.nwk" <<'EOF'
g((A:1,B:1),(C:1,D:1);
EOF
: > "$TMP/fixtures/empty.nwk"
cat > "$TMP/fixtures/two_trees.nwk" <<'EOF'
g((A:1,B:1),(C:1,D:1));((A:1,C:1),(B:1,D:1));
EOF
for case in unknown duplicate_taxon unbalanced empty two_trees; do
    mkdir -p "$TMP/invalid/$case"
    expect_failure "invalid_$case" bash -c "cd '$TMP/invalid/$case' && perl '$SA' --mode matrix --species '$TMP/fixtures/species.nwk' --gene '$TMP/fixtures/$case.nwk' --label x"
    [[ ! -e "$TMP/invalid/$case/x.matrix_with_fuse.txt" ]]
    [[ ! -e "$TMP/invalid/$case/x_splits" ]]
done
grep -q "Unknown taxon 'X'" "$TMP/invalid_unknown.stderr"
grep -q "Duplicate terminal taxon 'A'" "$TMP/invalid_duplicate_taxon.stderr"
pass 'SA-006/007 taxon-set and grammar validation'

# SA-008: missing root half-edge evidence cannot collapse to a numeric value.
cat > "$TMP/fixtures/root_evidence.nwk" <<'EOF'
g_partial((A:1,B:1):2,(C:1,D:1));
g_degree3(A:1,B:1,(C:1,D:1));
g_sentinel((A:1,B:1):2,(C:1,D:1):3);
EOF
mkdir -p "$TMP/sa008"
(
    cd "$TMP/sa008"
    perl "$SA" --mode matrix --species "$TMP/fixtures/species.nwk" --gene "$TMP/fixtures/root_evidence.nwk" --label x
) >/dev/null 2>&1
grep -q $'^g_partial\t1\t1\t1\t1\tNA$' "$TMP/sa008/x.matrix_no_fuse.txt"
grep -q $'^g_degree3\t1\t1\t1\t1\tNA$' "$TMP/sa008/x.matrix_no_fuse.txt"
grep -q $'^g_sentinel\t1\t1\t1\t1\t5$' "$TMP/sa008/x.matrix_no_fuse.txt"
pass 'SA-008 complete-evidence root contraction'

# SA-009: aliases may change with serialization; canonical dictionaries remain explicit.
cat > "$TMP/fixtures/species_order_a.nwk" <<'EOF'
((A,B),((C,D),(E,F)));
EOF
cat > "$TMP/fixtures/species_order_b.nwk" <<'EOF'
(((F,E),(D,C)),(B,A));
EOF
cat > "$TMP/fixtures/six_taxa_gene.nwk" <<'EOF'
g((A:1,B:2):10,((C:3,D:4):20,(E:5,F:6):30):40);
EOF
for variant in a b; do
    mkdir -p "$TMP/sa009/$variant"
    (
        cd "$TMP/sa009/$variant"
        perl "$SA" --mode matrix --species "$TMP/fixtures/species_order_${variant}.nwk" --gene "$TMP/fixtures/six_taxa_gene.nwk" --label x
    ) >/dev/null 2>&1
    cut -f2 "$TMP/sa009/$variant/x.primitive_axis.tsv" | sort > "$TMP/sa009/$variant/splits.sorted"
done
cmp -s "$TMP/sa009/a/splits.sorted" "$TMP/sa009/b/splits.sorted"
if cmp -s "$TMP/sa009/a/x.primitive_axis.tsv" "$TMP/sa009/b/x.primitive_axis.tsv"; then
    fail 'SA-009 fixture did not exercise serialization-local aliases'
fi
grep -q 'serialization-local aliases' "$ROOT/README.md"
pass 'SA-009 honest serialization-local B-alias contract'

# SOL-010: textual B headers are insufficient; canonical ledgers must match.
cat > "$TMP/fixtures/species_reordered.nwk" <<'EOF'
((C,D),(A,B));
EOF
cat > "$TMP/fixtures/shared_genes.nwk" <<'EOF'
g(B:2,(C:4,D:5):6);
s((A:1,B:2):3,(C:4,D:5):6);
EOF
mkdir -p "$TMP/sol010/free" "$TMP/sol010/fix" "$TMP/sol010/final"
(
    cd "$TMP/sol010/free"
    perl "$SA" --mode matrix --species "$TMP/fixtures/species.nwk" --gene "$TMP/fixtures/shared_genes.nwk" --label free
    cd "$TMP/sol010/fix"
    perl "$SA" --mode matrix --species "$TMP/fixtures/species_reordered.nwk" --gene "$TMP/fixtures/shared_genes.nwk" --label fix
) >/dev/null 2>&1
expect_failure sol010 bash -c "cd '$TMP/sol010/final' && perl '$SA' --mode finalize --free '$TMP/sol010/free/free.matrix_with_fuse.txt' --fix '$TMP/sol010/fix/fix.matrix_with_fuse.txt' --free_manifest '$TMP/sol010/free/free.run_manifest.json' --fix_manifest '$TMP/sol010/fix/fix.run_manifest.json' --final_label final"
grep -q 'Primitive axis mismatch' "$TMP/sol010.stderr"
[[ ! -e "$TMP/sol010/final/final.free.na_classified.txt" ]]
[[ ! -e "$TMP/sol010/free/free.matrix_with_fuse.na_fuse.txt" ]]

# A Support tree from the wrong ledger is rejected before output publication.
mkdir -p "$TMP/support_mismatch"
expect_failure support_mismatch bash -c "cd '$TMP/support_mismatch' && perl '$SA' --mode finalize --free '$TMP/sa001/single/free.matrix_with_fuse.txt' --fix '$TMP/sa001/single/fix.matrix_with_fuse.txt' --free_manifest '$TMP/sa001/single/free.run_manifest.json' --fix_manifest '$TMP/sa001/single/fix.run_manifest.json' --species_tree '$TMP/sol010/fix/species_tree.forSplit.nwk' --final_label final"
[[ ! -e "$TMP/support_mismatch/final.support_b.txt" ]]
pass 'SOL-010 ledger binding and Support-tree verification'

printf '[PASS] Phase 3 contract regression suite\n'
