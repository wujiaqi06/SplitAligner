#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
SA="$ROOT/SplitAligner.pl"
TMP="$(mktemp -d "${TMPDIR:-/tmp}/splitaligner_phase3_codec_XXXXXX")"
trap 'rm -rf "$TMP"' EXIT

fail() { printf '[FAIL] %s\n' "$1" >&2; exit 1; }
pass() { printf '[PASS] %s\n' "$1"; }

run_period_case() {
    local name="$1"
    local special_taxon="$2"
    local key_form="${3:-structured}"
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
    if [[ "$key_form" == structured ]]; then
        grep -q $'\tHX1:' "$dir/x.primitive_axis.tsv" \
            || fail "$name did not use the structured canonical split form"
    elif grep -q $'\tHX1:' "$dir/x.primitive_axis.tsv"; then
        fail "$name unnecessarily changed canonical format"
    fi
}

# Direct P3-013 counterexample: these are different unordered bipartitions.
env LC_ALL=C LANG=C perl -I"$ROOT/lib" \
    -MSplitAligner::Newick=canonical_split,canonicalize_split,split_taxon_sets \
    -e '
        my @universe = (".A", ".A.", ".B", "B", "BA");
        my @left1 = (".A.", "B");
        my %left1 = map { $_ => 1 } @left1;
        my @right1 = grep { !$left1{$_} } @universe;
        my @left2 = (".A", ".B");
        my %left2 = map { $_ => 1 } @left2;
        my @right2 = grep { !$left2{$_} } @universe;
        my $k1 = canonical_split(\@left1, \@right1);
        my $k2 = canonical_split(\@left2, \@right2);
        die "P3-013 collision remains\n" if $k1 eq $k2;
        for my $key ($k1, $k2) {
            die "counterexample did not use HX1\n" unless index($key, "HX1:") == 0;
            my ($left, $right) = split_taxon_sets($key);
            die "counterexample round-trip failed\n"
                unless canonical_split($left, $right) eq $key;
            die "counterexample side reversal failed\n"
                unless canonical_split($right, $left) eq $key;
            die "canonicalize_split is unstable\n"
                unless canonicalize_split($key) eq $key;
        }
    '
pass 'P3-013 direct canonical split collision is closed'

run_period_case trailing_period 'A.'
run_period_case leading_period '.A'
run_period_case period_only '.'
run_period_case internal_and_trailing_period 'NC.045.'
run_period_case leading_period_accession '.NC045'
pass 'leading and trailing period taxa retain five independent quartet coordinates'

# A single internal period is provably legacy-safe and remains readable.
run_period_case ordinary_internal_period 'NC_045512.2' legacy
pass 'ordinary internal period remains backward-compatible'

# Exhaust every nontrivial subset over a finite legal-label universe. The
# signature uses integer membership, independently of the production codec.
env LC_ALL=C LANG=C perl -I"$ROOT/lib" \
    -MSplitAligner::Newick=canonical_split,canonicalize_split,split_taxon_sets \
    -e '
        my @u = ("A", "B", ".C", "D.", "N.E", "NC.045.", ".NC045", ".");
        my %seen;
        my $all = (1 << @u) - 1;
        for my $mask (1 .. $all - 1) {
            my $other = $all ^ $mask;
            my @left  = map { $u[$_] } grep { $mask  & (1 << $_) } 0 .. $#u;
            my @right = map { $u[$_] } grep { $other & (1 << $_) } 0 .. $#u;
            my $low = $mask < $other ? $mask : $other;
            my $high = $mask < $other ? $other : $mask;
            my $signature = "$low/$high";
            my $key = canonical_split(\@left, \@right);
            die "distinct bipartitions collided at $key\n"
                if exists $seen{$key} && $seen{$key} ne $signature;
            $seen{$key} = $signature;
            my ($decoded_left, $decoded_right) = split_taxon_sets($key);
            die "decode/re-encode instability at $key\n"
                unless canonical_split($decoded_left, $decoded_right) eq $key;
            die "side reversal instability at $key\n"
                unless canonical_split(\@right, \@left) eq $key;
            die "canonicalize instability at $key\n"
                unless canonicalize_split($key) eq $key;
        }
        die "unexpected bipartition count\n" unless scalar(keys %seen) == 127;
    '
pass 'finite exhaustive hybrid-split injectivity and round-trip checks'

printf '[PASS] Phase 3 codec-closure regression suite\n'
