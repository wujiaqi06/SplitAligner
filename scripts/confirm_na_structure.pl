#!/usr/bin/env perl

# ==============================================================================
# Script:      confirm_na_structure.pl
# Author:      Jiaqi Wu (wujiaqi@hiroshima-u.ac.jp)
# Affiliation: Graduate School of Integrated Sciences for Life, Hiroshima University, Japan
# Copyright:   (c) 2026 Jiaqi Wu. All rights reserved.
#
# Description:
#   Reconciles the "free-topology" and "fixed-topology" gene branch matrices to
#   distinguish the biological nature of missing data (NA).
#
#   Differential diagnosis uses explicit S/D/F/U coordinate-state sidecars.
#   Numeric missingness alone is never interpreted as structural absence.
#
# Inputs:
#   --fix  <fix_matrix>   : matrix generated from fixed-topology gene trees
#   --free <free_matrix>  : matrix generated from free-topology gene trees
#   --species_tree <species_tree.forSplit.nwk> : optional species tree for
#                                                branch-wise Support annotation
#
# Outputs (prefix = -o):
#   <prefix>.fix.na_classified.txt   : fixed matrix with NA_struct applied
#   <prefix>.free.na_classified.txt  : free matrix with NA_struct / NA_topo applied
#   <prefix>.support_b.txt           : optional branch-wise Support summary
#   <species_prefix>.support_b.nwk   : optional species tree in standard Newick
#                                      format with internal-node Support values
#
# Usage:
#   perl confirm_na_structure.pl --fix fix.matrix.txt --free free.matrix.txt \
#     --species_tree species_tree.forSplit.nwk -o out_prefix
# ==============================================================================

use strict;
use warnings;
use Getopt::Long;
use File::Basename qw(basename);
use File::Spec;
use FindBin qw($RealBin);
use lib "$RealBin/../lib";
use SplitAligner::Newick qw(is_numeric_value);
use SplitAligner::CoordinateState qw(
    read_matrix_table
    read_state_matrix
    validate_state_for_matrix
);
use SplitAligner::TextIO qw(
    configure_utf8_stdio
    decode_argv_utf8
    read_utf8_file
    read_utf8_lines
    write_utf8_file
);

decode_argv_utf8();
configure_utf8_stdio();

my ($fix_matrix_path, $free_matrix_path, $fix_state_path, $free_state_path,
    $species_tree_path, $out_prefix);

GetOptions(
    'fix=s'          => \$fix_matrix_path,
    'free=s'         => \$free_matrix_path,
    'fix_state=s'    => \$fix_state_path,
    'free_state=s'   => \$free_state_path,
    'species_tree=s' => \$species_tree_path,
    'o=s'            => \$out_prefix,
) or die "[ERROR] Invalid command line arguments.\n";

die "Usage: $0 --fix <fix_matrix> --free <free_matrix> --fix_state <fix_state.tsv> --free_state <free_state.tsv> [--species_tree <species_tree.forSplit.nwk>] -o <out_prefix>\n"
    unless defined $fix_matrix_path && defined $free_matrix_path
        && defined $fix_state_path && defined $free_state_path && defined $out_prefix;

my $out_fix  = "$out_prefix.fix.na_classified.txt";
my $out_free = "$out_prefix.free.na_classified.txt";
my $out_support = "$out_prefix.support_b.txt";

my $fix_matrix = read_matrix_table($fix_matrix_path);
my $free_matrix = read_matrix_table($free_matrix_path);
my $fix_state = read_state_matrix($fix_state_path);
my $free_state = read_state_matrix($free_state_path);
validate_state_for_matrix(
    state => $fix_state, matrix => $fix_matrix, axis => $fix_matrix->{primitive},
    role => 'FIX finalization input', allow_na_fuse => 1, allow_classified => 1,
);
validate_state_for_matrix(
    state => $free_state, matrix => $free_matrix, axis => $free_matrix->{primitive},
    role => 'FREE finalization input', allow_na_fuse => 1, allow_classified => 1,
);

my @fix_branches = @{$fix_matrix->{branches}};
my @free_branches = @{$free_matrix->{branches}};
my $fix_header_line = join("\t", 'gene', @fix_branches);
my $free_header_line = join("\t", 'gene', @free_branches);
my %fix_row = map {
    $_ => join("\t", @{$fix_matrix->{rows}{$_}})
} @{$fix_matrix->{genes}};
my %free_row = map {
    $_ => join("\t", @{$free_matrix->{rows}{$_}})
} @{$free_matrix->{genes}};

# -------------------------
# Sanity check: branch axis
# -------------------------
my @fix_primitive_branches  = primitive_branch_axis(\@fix_branches);
my @free_primitive_branches = primitive_branch_axis(\@free_branches);

die "[ERROR] No primitive species-tree branches were found in FIX matrix header.\n"
    unless @fix_primitive_branches;
die "[ERROR] No primitive species-tree branches were found in FREE matrix header.\n"
    unless @free_primitive_branches;

my $fix_axis  = join("\t", @fix_primitive_branches);
my $free_axis = join("\t", @free_primitive_branches);

if ($fix_axis eq $free_axis) {
    print STDERR "[INFO] Primitive branch axes matched between fix and free matrices.\n";
} else {
    my $fix_count  = scalar(@fix_primitive_branches);
    my $free_count = scalar(@free_primitive_branches);
    die "[ERROR] Primitive branch axes DO NOT match between fix and free matrices. ".
        "FIX has $fix_count primitive branches, FREE has $free_count. ".
        "Support/discordance classification requires the same primitive species-tree branch axis and order.\n";
}

# -------------------------
# Shared gene check
# -------------------------
my @shared_genes = grep { exists $free_row{$_} } @{$fix_matrix->{genes}};

my $n_fix       = scalar(keys %fix_row);
my $n_free      = scalar(keys %free_row);
my $n_shared    = scalar(@shared_genes);
my $n_fix_only  = $n_fix  - $n_shared;
my $n_free_only = $n_free - $n_shared;

print STDERR "[INFO] Genes in FIX matrix : $n_fix\n";
print STDERR "[INFO] Genes in FREE matrix: $n_free\n";
print STDERR "[INFO] Shared genes        : $n_shared\n";
print STDERR "[INFO] FIX-only genes      : $n_fix_only\n";
print STDERR "[INFO] FREE-only genes     : $n_free_only\n";

die "[ERROR] No shared genes were found between FIX and FREE matrices. Please check whether the two inputs were generated from comparable gene sets and whether gene IDs are consistent.\n"
    if $n_shared == 0;

my @fix_only_genes  = grep { !exists $free_row{$_} } @{$fix_matrix->{genes}};
my @free_only_genes = grep { !exists $fix_row{$_} } @{$free_matrix->{genes}};
my @gene_id_alias_matches = detect_hyphen_underscore_aliases(\@fix_only_genes, \@free_only_genes);

if (@gene_id_alias_matches) {
    my $preview_count = @gene_id_alias_matches < 10 ? scalar(@gene_id_alias_matches) : 10;
    my @preview = map {
        $_->{fix_gene} . " <-> " . $_->{free_gene}
    } @gene_id_alias_matches[0 .. $preview_count - 1];

    die "[ERROR] FIX and FREE matrices contain gene IDs that differ only by underscore/hyphen normalization. ".
        "This would silently drop genes during finalize because only exact shared IDs are retained. ".
        "Detected ".scalar(@gene_id_alias_matches)." likely alias pair(s), e.g. ".
        join(", ", @preview).
        ". Please harmonize gene IDs in the input gene-tree files before rerunning SplitAligner.\n";
}

my %support_stats = initialize_support_stats(
    fix_rows     => \%fix_row,
    free_rows    => \%free_row,
    fix_branches => \@fix_branches,
    free_branches => \@free_branches,
    shared_genes => \@shared_genes,
);

my %fix_state_index = map { $fix_state->{axis}[$_] => $_ } 0 .. $#{$fix_state->{axis}};
my %free_state_index = map { $free_state->{axis}[$_] => $_ } 0 .. $#{$free_state->{axis}};
my @out_fix_lines = ($fix_header_line);
my @out_free_lines = ($free_header_line);

for my $gene (@shared_genes) {
    my @free_vals = split(/\t/, $free_row{$gene}, -1);
    my %free = map { $free_branches[$_] => $free_vals[$_] } 0 .. $#free_branches;

    my @fix_vals = split(/\t/, $fix_row{$gene}, -1);
    my %fix = map { $fix_branches[$_] => $fix_vals[$_] } 0 .. $#fix_branches;
    my @fix_states = @{$fix_state->{rows}{$gene}};
    my @free_states = @{$free_state->{rows}{$gene}};

    for my $b (@free_branches) {
        next unless exists $fix{$b};

        my $fix_code = $fix_states[ $fix_state_index{$b} ];
        my $free_code = $free_states[ $free_state_index{$b} ];
        $fix{$b} = 'NA_struct' if $fix{$b} eq 'NA' && $fix_code eq 'S';

        if ($free{$b} eq 'NA' && $free_code eq 'S') {
            $free{$b} = 'NA_struct';
        }
        elsif ($free{$b} eq 'NA' && $free_code eq 'U'
            && $fix_code eq 'D' && is_numeric_branch_evidence($fix{$b})) {
            $free{$b} = 'NA_topo';
        }
    }

    push @out_fix_lines, join("\t", $gene, map { $fix{$_} // 'NA' } @fix_branches);
    push @out_free_lines, join("\t", $gene, map { $free{$_} // 'NA' } @free_branches);
}

write_utf8_file($out_fix, join("\n", @out_fix_lines) . "\n");
write_utf8_file($out_free, join("\n", @out_free_lines) . "\n");

print STDERR "[INFO] Wrote: $out_fix\n";
print STDERR "[INFO] Wrote: $out_free\n";

if (defined $species_tree_path && $species_tree_path ne '') {
    die "[ERROR] Species tree file not found: $species_tree_path\n" unless -e $species_tree_path;
    write_support_outputs(
        stats             => \%support_stats,
        branches          => \@fix_primitive_branches,
        out_support       => $out_support,
        species_tree_path => $species_tree_path,
    );
}

sub write_support_outputs {
    my %arg = @_;

    my $stats_ref         = $arg{stats};
    my $branches_ref      = $arg{branches};
    my $out_support_path  = $arg{out_support};
    my $tree_path         = $arg{species_tree_path};
    my ($branch_type_ref, $duplicate_loser_ref) = read_branch_support_metadata($tree_path);
    my %branch_type_for = %{$branch_type_ref};
    my $tree_text = read_tree_text($tree_path);
    my @support_lines = (join(
        "\t",
        qw(branch_id branch_type n_shared_genes n_fix_non_na n_free_non_na support_percent discordance_percent)
    ));

    my %support_value_for;
    for my $b (@{$branches_ref}) {
        my $s = $stats_ref->{$b} || {};
        my $n_fix_non_na  = $s->{n_fix_non_na}  || 0;
        my $n_free_non_na = $s->{n_free_non_na} || 0;
        my $support       = $n_fix_non_na > 0 ? (100 * $n_free_non_na / $n_fix_non_na) : 0;
        my $discordance   = $n_fix_non_na > 0 ? (100 - $support) : 0;

        $support_value_for{$b} = sprintf('%.10f', $support);

        push @support_lines, join(
            "\t",
            $b,
            (exists $branch_type_for{$b} ? $branch_type_for{$b} : 'NA'),
            $s->{n_shared_genes} || 0,
            $n_fix_non_na,
            $n_free_non_na,
            sprintf('%.10f', $support),
            sprintf('%.10f', $discordance),
        );
    }

    for my $loser (sort keys %{$duplicate_loser_ref}) {
        my $winner = $duplicate_loser_ref->{$loser};
        next unless defined $winner && exists $support_value_for{$winner};
        $support_value_for{$loser} = $support_value_for{$winner};
    }

    $tree_text = standardize_support_tree($tree_text, \%support_value_for);

    my $tree_base = basename($tree_path);
    if ($tree_base =~ /\.forSplit\.nwk$/) {
        $tree_base =~ s/\.forSplit\.nwk$/.support_b.nwk/;
    }
    elsif ($tree_base =~ /\.nwk$/) {
        $tree_base =~ s/\.nwk$/.support_b.nwk/;
    }
    else {
        $tree_base .= '.support_b.nwk';
    }

    write_utf8_file($out_support_path, join("\n", @support_lines) . "\n");
    write_utf8_file($tree_base, $tree_text);

    print STDERR "[INFO] Wrote: $out_support_path\n";
    print STDERR "[INFO] Wrote: $tree_base\n";
}

sub primitive_branch_axis {
    my ($branches_ref) = @_;
    return grep { defined $_ && $_ =~ /^B\d+$/ } @{$branches_ref};
}

sub read_branch_support_metadata {
    my ($tree_path) = @_;

    my $map_path = derive_branch_map_path($tree_path);
    return ({}, {}) unless defined $map_path && -e $map_path;

    my %branch_type_for;
    my %duplicate_loser_of;
    my $line_no = 0;
    for my $line (@{ read_utf8_lines($map_path) }) {
        next if $line =~ /^\s*$/;
        $line_no++;
        next if $line_no == 1; # header

        my @f = split(/\t/, $line, -1);
        next unless @f >= 3;
        my ($branch_id, undef, $branch_type, $note) = @f[0, 1, 2, 3];
        next unless defined $branch_id && $branch_id ne '';
        $branch_type_for{$branch_id} = $branch_type ne '' ? $branch_type : 'NA';
        if (defined $note && $note =~ /duplicate_unrooted_split_loser_of=(B\d+)/) {
            $duplicate_loser_of{$branch_id} = $1;
        }
    }

    print STDERR "[INFO] Loaded branch types from $map_path\n";
    return (\%branch_type_for, \%duplicate_loser_of);
}

sub derive_branch_map_path {
    my ($tree_path) = @_;

    my ($vol, $dir, $file) = File::Spec->splitpath($tree_path);
    my @candidates;
    if ($file =~ /(.*)\.forSplit\.nwk$/) {
        push @candidates, File::Spec->catpath($vol, $dir, "$1.branch_map.txt");
    }
    if ($file =~ /(.*)\.nwk$/) {
        push @candidates, File::Spec->catpath($vol, $dir, "$1.branch_map.txt");
    }
    push @candidates, File::Spec->catpath($vol, $dir, 'species_tree.branch_map.txt');

    my %seen;
    for my $candidate (@candidates) {
        next if $seen{$candidate}++;
        return $candidate if -e $candidate;
    }

    print STDERR "[WARN] Branch map not found alongside species tree; branch_type will be written as NA in support table.\n";
    return undef;
}

sub read_tree_text {
    my ($path) = @_;
    my $text = read_utf8_file($path);
    die "[ERROR] Species tree file is empty: $path\n" unless defined $text && $text ne '';
    return $text;
}

sub standardize_support_tree {
    my ($tree_text, $support_ref) = @_;

    # Remove the optional leading record id before the first '('.
    $tree_text =~ s/^\s*[^\(]+(?=\()//;

    # Tips keep only the taxon name; the forSplit branch ids are removed.
    $tree_text =~ s/([^\s(),:;]+):B\d+/$1/g;

    # Internal branches are converted from forSplit labels like "):B305"
    # into standard Newick internal-node labels like ")95.000000".
    $tree_text =~ s/\):([A-Z]\d+)(?![\w|])/')' . support_label($1, $support_ref)/ge;

    return $tree_text;
}

sub support_label {
    my ($branch_id, $support_ref) = @_;
    my $support = exists $support_ref->{$branch_id} ? $support_ref->{$branch_id} : '0.0000000000';
    return $support;
}

sub initialize_support_stats {
    my %arg = @_;

    my $fix_rows_ref      = $arg{fix_rows};
    my $free_rows_ref     = $arg{free_rows};
    my $fix_branches_ref  = $arg{fix_branches};
    my $free_branches_ref = $arg{free_branches};
    my $shared_genes_ref  = $arg{shared_genes};

    my %stats;
    for my $b (@{$fix_branches_ref}) {
        $stats{$b} = {
            n_shared_genes => scalar(@{$shared_genes_ref}),
            n_fix_non_na => 0,
            n_free_non_na => 0,
        };
    }

    my %free_index = map { $free_branches_ref->[$_] => $_ } 0 .. $#{$free_branches_ref};

    for my $gene (@{$shared_genes_ref}) {
        next unless exists $fix_rows_ref->{$gene} && exists $free_rows_ref->{$gene};

        my @fix_vals  = split(/\t/, $fix_rows_ref->{$gene}, -1);
        my @free_vals = split(/\t/, $free_rows_ref->{$gene}, -1);

        for my $i (0 .. $#{$fix_branches_ref}) {
            my $branch    = $fix_branches_ref->[$i];
            my $fix_value = defined $fix_vals[$i] ? $fix_vals[$i] : 'NA';

            $stats{$branch}{n_fix_non_na}++ if is_numeric_branch_evidence($fix_value);

            next unless exists $free_index{$branch};
            my $free_idx   = $free_index{$branch};
            my $free_value = defined $free_vals[$free_idx] ? $free_vals[$free_idx] : 'NA';
            $stats{$branch}{n_free_non_na}++ if is_numeric_branch_evidence($free_value);
        }
    }

    return %stats;
}

sub detect_hyphen_underscore_aliases {
    my ($fix_only_ref, $free_only_ref) = @_;

    my %free_by_normalized;
    for my $gene (@{$free_only_ref}) {
        my $normalized = normalize_gene_id_for_alias_check($gene);
        next unless defined $normalized;
        push @{ $free_by_normalized{$normalized} }, $gene;
    }

    my @matches;
    for my $gene (@{$fix_only_ref}) {
        my $normalized = normalize_gene_id_for_alias_check($gene);
        next unless defined $normalized;
        next unless exists $free_by_normalized{$normalized};
        for my $free_gene (@{ $free_by_normalized{$normalized} }) {
            next if $gene eq $free_gene;
            push @matches, {
                fix_gene  => $gene,
                free_gene => $free_gene,
            };
        }
    }

    return @matches;
}

sub normalize_gene_id_for_alias_check {
    my ($gene) = @_;
    return undef unless defined $gene;
    my $normalized = $gene;
    $normalized =~ tr/_-/__/;
    return $normalized;
}

sub is_numeric_branch_evidence {
    my ($value) = @_;
    return is_numeric_value($value);
}
