#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
SA="$ROOT/SplitAligner.pl"
TMP="$(mktemp -d "${TMPDIR:-/tmp}/splitaligner_phase3_meta_XXXXXX")"
trap 'rm -rf "$TMP"' EXIT

cat > "$TMP/species.nwk" <<'EOF'
((A,B),((C,D),(E,F)));
EOF
cat > "$TMP/order1.nwk" <<'EOF'
gene_z((A:1,B:2):3,((C:4,D:5):6,(E:7,F:8):9):10);
gene_a((A:2,C:3):4,((B:5,D:6):7,(E:8,F:9):10):11);
EOF
cat > "$TMP/order2.nwk" <<'EOF'
gene_a((A:2,C:3):4,((B:5,D:6):7,(E:8,F:9):10):11);
gene_z((A:1,B:2):3,((C:4,D:5):6,(E:7,F:8):9):10);
EOF

for seed in 1 2 3 4 5 6 7 8; do
    for order in 1 2; do
        out="$TMP/seed_${seed}_order_${order}"
        mkdir -p "$out"
        (
            cd "$out"
            PERL_HASH_SEED="$seed" PERL_PERTURB_KEYS=2 \
                perl "$SA" --mode matrix --species "$TMP/species.nwk" --gene "$TMP/order${order}.nwk" --label x
        ) >/dev/null 2>&1
    done
done

REF="$TMP/seed_1_order_1"
for seed in 1 2 3 4 5 6 7 8; do
    for order in 1 2; do
        out="$TMP/seed_${seed}_order_${order}"
        cmp -s "$REF/x.matrix_no_fuse.txt" "$out/x.matrix_no_fuse.txt"
        cmp -s "$REF/x.matrix_with_fuse.txt" "$out/x.matrix_with_fuse.txt"
        cmp -s "$REF/x.primitive_state.tsv" "$out/x.primitive_state.tsv"
        cmp -s "$REF/x.primitive_axis.tsv" "$out/x.primitive_axis.tsv"
        cmp -s "$REF/species_tree.splits.txt" "$out/species_tree.splits.txt"
    done
done

cat > "$TMP/numeric_a.nwk" <<'EOF'
g((A:1.0,B:2.0):3.0,((C:4.0,D:5.0):6.0,(E:7.0,F:8.0):9.0):10.0);
EOF
cat > "$TMP/numeric_b.nwk" <<'EOF'
g((A:1.,B:+2e0):3e0,((C:.4e1,D:5):+6.0,(E:7e0,F:8.):9):1e1);
EOF
for variant in a b; do
    mkdir -p "$TMP/numeric_$variant"
    (
        cd "$TMP/numeric_$variant"
        perl "$SA" --mode matrix --species "$TMP/species.nwk" --gene "$TMP/numeric_${variant}.nwk" --label x
    ) >/dev/null 2>&1
done

perl -e '
    use strict; use warnings;
    sub row { my ($p)=@_; open my $f,"<",$p or die $!; <$f>; my $r=<$f>; chomp $r; my @x=split /\t/,$r,-1; shift @x; return @x; }
    my @a=row($ARGV[0]); my @b=row($ARGV[1]);
    die "column count mismatch\n" unless @a==@b;
    for my $i (0..$#a) { die "numeric mismatch at $i: $a[$i] vs $b[$i]\n" if abs((0+$a[$i])-(0+$b[$i])) > 1e-12; }
' "$TMP/numeric_a/x.matrix_no_fuse.txt" "$TMP/numeric_b/x.matrix_no_fuse.txt"
cmp -s "$TMP/numeric_a/x.primitive_state.tsv" "$TMP/numeric_b/x.primitive_state.tsv"

printf '[PASS] Phase 3 gene-order, hash-seed, numeric-format, and state-sidecar metamorphic regression\n'
