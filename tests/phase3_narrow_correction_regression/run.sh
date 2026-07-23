#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
SA="$ROOT/SplitAligner.pl"
TMP="$(mktemp -d "${TMPDIR:-/tmp}/splitaligner_phase3_narrow_XXXXXX")"
trap 'rm -rf "$TMP"' EXIT

fail() { printf '[FAIL] %s\n' "$1" >&2; exit 1; }
pass() { printf '[PASS] %s\n' "$1"; }

run_delimiter_case() {
    local name="$1"
    local special_taxon="$2"
    local dir="$TMP/$name"
    mkdir -p "$dir"
    printf '((A,B),(%s,C));\n' "$special_taxon" > "$dir/species.nwk"
    printf 'g((A:1,B:2):10,(%s:3,C:4):20);\n' "$special_taxon" > "$dir/genes.nwk"
    (
        cd "$dir"
        env LC_ALL=C LANG=C perl "$SA" --mode matrix \
            --species species.nwk --gene genes.nwk --label x
    ) >"$dir/stdout.txt" 2>"$dir/stderr.txt"
    [[ "$(wc -l < "$dir/x.primitive_axis.tsv" | tr -d ' ')" == 5 ]] \
        || fail "$name did not retain five quartet coordinates"
    grep -qx $'g\t1\t2\t3\t4\t30' "$dir/x.matrix_no_fuse.txt" \
        || fail "$name merged or misplaced numeric edge values"
}

# P3-011 minimal collision and additional legal delimiter-bearing labels.
run_delimiter_case dotdot 'A..B'
run_delimiter_case double_pipe 'A||B'
run_delimiter_case single_pipe 'A|B'
run_delimiter_case unicode_comma_like 'α，β'
pass 'P3-011 delimiter-bearing and UTF-8 taxon labels remain injective'

# Ordinary single periods remain in the legacy-readable canonical format.
run_delimiter_case ordinary_period 'NC_045512.2'
if grep -q $'\tHX1:' "$TMP/ordinary_period/x.primitive_axis.tsv"; then
    fail 'ordinary single-period label unnecessarily changed canonical format'
fi
pass 'ordinary single-period taxon labels remain backward-compatible'

# Canonical set keys are injective over the tested legal-label universe, and
# canonical split orientation is invariant to side and child order.
env LC_ALL=C LANG=C perl -I"$ROOT/lib" -CS -Mutf8 \
    -MSplitAligner::Newick=canonical_split,canonical_taxon_set_key,split_taxon_sets \
    -e '
        my @u = ("A", "A..B", "A||B", "A|B", "NC_045512.2", "α，β", "雪");
        my %seen;
        for my $mask (0 .. (1 << @u) - 1) {
            my @set = map { $u[$_] } grep { $mask & (1 << $_) } 0 .. $#u;
            my $key = canonical_taxon_set_key(@set);
            die "set-key collision: $key\n" if exists $seen{$key};
            $seen{$key} = $mask;
        }
        my $a = canonical_split(["A..B", "A"], ["C", "B"]);
        my $b = canonical_split(["B", "C"], ["A", "A..B"]);
        die "unrooted orientation mismatch\n" unless $a eq $b;
        my ($left, $right) = split_taxon_sets($a);
        my @roundtrip = sort (@{$left}, @{$right});
        die "round-trip mismatch\n" unless join("\0", @roundtrip) eq join("\0", sort("A", "A..B", "B", "C"));
    '
pass 'property-style canonical set-key and split round-trip checks'

# Colon and ASCII whitespace remain grammar delimiters, not taxon content.
mkdir -p "$TMP/rejected_colon" "$TMP/rejected_space"
cat > "$TMP/rejected_colon/species.nwk" <<'NWK'
((A,B),(A:bad,C));
NWK
cat > "$TMP/rejected_space/species.nwk" <<'NWK'
((A,B),(A B,C));
NWK
cat > "$TMP/rejected_colon/genes.nwk" <<'NWK'
g((A:1,B:2):3,(C:4,D:5):6);
NWK
cp "$TMP/rejected_colon/genes.nwk" "$TMP/rejected_space/genes.nwk"
for case in rejected_colon rejected_space; do
    if (
        cd "$TMP/$case"
        env LC_ALL=C LANG=C perl "$SA" --mode matrix \
            --species species.nwk --gene genes.nwk --label x
    ) >"$TMP/$case/stdout.txt" 2>"$TMP/$case/stderr.txt"; then
        fail "$case was unexpectedly accepted"
    fi
    [[ ! -e "$TMP/$case/x.matrix_no_fuse.txt" ]] || fail "$case created final output"
    [[ ! -d "$TMP/$case/x_splits" ]] || fail "$case created split output"
done
pass 'colon and ASCII whitespace taxon content fail during preflight'

# P3-012 reserved label fails before a temporary workspace or final output.
mkdir -p "$TMP/reserved_label"
cat > "$TMP/reserved_label/species.nwk" <<'NWK'
((A,B),(C,D));
NWK
cat > "$TMP/reserved_label/genes.nwk" <<'NWK'
g((A:1,B:2):3,(C:4,D:5):6);
NWK
if (
    cd "$TMP/reserved_label"
    env LC_ALL=C LANG=C perl "$SA" --mode matrix \
        --species species.nwk --gene genes.nwk --label species_tree
) >"$TMP/reserved_label/stdout.txt" 2>"$TMP/reserved_label/stderr.txt"; then
    fail 'reserved matrix label species_tree was unexpectedly accepted'
fi
grep -q -- "--label 'species_tree' is reserved" "$TMP/reserved_label/stderr.txt" \
    || fail 'reserved-label diagnostic was not explicit'
[[ "$(find "$TMP/reserved_label" -maxdepth 1 -name '.splitaligner-*' | wc -l | tr -d ' ')" == 0 ]] \
    || fail 'reserved label created a temporary computation directory'
[[ "$(find "$TMP/reserved_label" -maxdepth 1 -type f | wc -l | tr -d ' ')" == 4 ]] \
    || fail 'reserved label created unexpected final files'
pass 'P3-012 reserved matrix label fails before output creation'

printf '[PASS] Phase 3 narrow-correction regression suite\n'
