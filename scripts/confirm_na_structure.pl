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
#   Differential diagnosis (per gene × branch cell):
#     - NA in BOTH free and fixed matrices  -> NA_struct (structural / missing signal)
#     - value in fixed but NA in free       -> NA_topo  (topological discordance; e.g., ILS)
#     - otherwise                           -> keep original value/NA
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

my ($fix_matrix_path, $free_matrix_path, $species_tree_path, $out_prefix);

GetOptions(
    'fix=s'          => \$fix_matrix_path,
    'free=s'         => \$free_matrix_path,
    'species_tree=s' => \$species_tree_path,
    'o=s'            => \$out_prefix,
) or die "[ERROR] Invalid command line arguments.\n";

die "Usage: $0 --fix <fix_matrix> --free <free_matrix> [--species_tree <species_tree.forSplit.nwk>] -o <out_prefix>\n"
    unless defined $fix_matrix_path && defined $free_matrix_path && defined $out_prefix;

my $out_fix  = "$out_prefix.fix.na_classified.txt";
my $out_free = "$out_prefix.free.na_classified.txt";
my $out_support = "$out_prefix.support_b.txt";

# -------------------------
# Read FIX matrix
# -------------------------
my (%fix_row, @fix_branches, $fix_header_line);
open(my $FIX, '<', $fix_matrix_path) or die "[ERROR] Cannot open $fix_matrix_path: $!\n";

my $line_no = 0;
while (my $line = <$FIX>) {
    chomp $line;
    next if $line =~ /^\s*$/;
    $line_no++;

    my @f = split(/\t/, $line, -1);

    if ($line_no == 1) {
        $fix_header_line = $line;
        @fix_branches = @f[1 .. $#f];
        next;
    }

    my $gene = $f[0];
    next unless defined $gene && $gene ne '';
    $fix_row{$gene} = join("\t", @f[1 .. $#f]);
}
close $FIX;

die "[ERROR] FIX matrix appears empty or missing header: $fix_matrix_path\n"
    unless defined $fix_header_line;

# -------------------------
# Read FREE matrix
# -------------------------
my (%free_row, @free_branches, $free_header_line);
open(my $FREE, '<', $free_matrix_path) or die "[ERROR] Cannot open $free_matrix_path: $!\n";

$line_no = 0;
while (my $line = <$FREE>) {
    chomp $line;
    next if $line =~ /^\s*$/;
    $line_no++;

    my @f = split(/\t/, $line, -1);

    if ($line_no == 1) {
        $free_header_line = $line;
        @free_branches = @f[1 .. $#f];
        next;
    }

    my $gene = $f[0];
    next unless defined $gene && $gene ne '';
    $free_row{$gene} = join("\t", @f[1 .. $#f]);
}
close $FREE;

die "[ERROR] FREE matrix appears empty or missing header: $free_matrix_path\n"
    unless defined $free_header_line;

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
my @shared_genes = grep { exists $free_row{$_} } sort keys %fix_row;

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

my @fix_only_genes  = grep { !exists $free_row{$_} } sort keys %fix_row;
my @free_only_genes = grep { !exists $fix_row{$_} } sort keys %free_row;
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

# -------------------------
# Open outputs only after validation
# -------------------------
open(my $OUT_FIX,  '>', $out_fix)  or die "[ERROR] Cannot write $out_fix: $!\n";
open(my $OUT_FREE, '>', $out_free) or die "[ERROR] Cannot write $out_free: $!\n";

print {$OUT_FIX}  "$fix_header_line\n";
print {$OUT_FREE} "$free_header_line\n";

my %support_stats = initialize_support_stats(
    fix_rows     => \%fix_row,
    free_rows    => \%free_row,
    fix_branches => \@fix_branches,
    free_branches => \@free_branches,
    shared_genes => \@shared_genes,
);

# -------------------------
# Classify NA types
# -------------------------
for my $gene (@shared_genes) {
    my @free_vals = split(/\t/, $free_row{$gene}, -1);
    my %free = map { $free_branches[$_] => $free_vals[$_] } 0 .. $#free_branches;

    my @fix_vals = split(/\t/, $fix_row{$gene}, -1);
    my %fix = map { $fix_branches[$_] => $fix_vals[$_] } 0 .. $#fix_branches;

    for my $b (@free_branches) {
        next unless exists $fix{$b};

        if ($fix{$b} eq 'NA' && $free{$b} eq 'NA') {
            $fix{$b}  = 'NA_struct';
            $free{$b} = 'NA_struct';
        }
        elsif (is_numeric_branch_evidence($fix{$b}) && $free{$b} eq 'NA') {
            $free{$b} = 'NA_topo';
        }
    }

    print {$OUT_FIX} $gene;
    for my $b (@fix_branches) {
        my $v = exists $fix{$b} ? $fix{$b} : 'NA';
        print {$OUT_FIX} "\t$v";
    }
    print {$OUT_FIX} "\n";

    print {$OUT_FREE} $gene;
    for my $b (@free_branches) {
        my $v = exists $free{$b} ? $free{$b} : 'NA';
        print {$OUT_FREE} "\t$v";
    }
    print {$OUT_FREE} "\n";
}

close $OUT_FIX;
close $OUT_FREE;

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

    open(my $SUP, '>', $out_support_path) or die "[ERROR] Cannot write $out_support_path: $!\n";
    print {$SUP} join(
        "\t",
        qw(branch_id branch_type n_shared_genes n_fix_non_na n_free_non_na support_percent discordance_percent)
    ), "\n";

    my %support_value_for;
    for my $b (@{$branches_ref}) {
        my $s = $stats_ref->{$b} || {};
        my $n_fix_non_na  = $s->{n_fix_non_na}  || 0;
        my $n_free_non_na = $s->{n_free_non_na} || 0;
        my $support       = $n_fix_non_na > 0 ? (100 * $n_free_non_na / $n_fix_non_na) : 0;
        my $discordance   = $n_fix_non_na > 0 ? (100 - $support) : 0;

        $support_value_for{$b} = sprintf('%.10f', $support);

        print {$SUP} join(
            "\t",
            $b,
            (exists $branch_type_for{$b} ? $branch_type_for{$b} : 'NA'),
            $s->{n_shared_genes} || 0,
            $n_fix_non_na,
            $n_free_non_na,
            sprintf('%.10f', $support),
            sprintf('%.10f', $discordance),
        ), "\n";
    }
    close $SUP;

    for my $loser (sort keys %{$duplicate_loser_ref}) {
        my $winner = $duplicate_loser_ref->{$loser};
        next unless defined $winner && exists $support_value_for{$winner};
        $support_value_for{$loser} = $support_value_for{$winner};
    }

    my $tree_text = read_tree_text($tree_path);
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

    open(my $TREE_OUT, '>', $tree_base) or die "[ERROR] Cannot write $tree_base: $!\n";
    print {$TREE_OUT} $tree_text;
    close $TREE_OUT;

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

    open(my $MAP, '<', $map_path) or die "[ERROR] Cannot open branch map $map_path: $!\n";
    my %branch_type_for;
    my %duplicate_loser_of;
    my $line_no = 0;
    while (my $line = <$MAP>) {
        chomp $line;
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
    close $MAP;

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
    open(my $IN, '<', $path) or die "[ERROR] Cannot open $path: $!\n";
    local $/ = undef;
    my $text = <$IN>;
    close $IN;
    die "[ERROR] Species tree file is empty: $path\n" unless defined $text && $text ne '';
    return $text;
}

sub standardize_support_tree {
    my ($tree_text, $support_ref) = @_;

    # Remove the optional leading record id before the first '('.
    $tree_text =~ s/^\s*[^\(]+(?=\()//;

    # Tips keep only the taxon name; the forSplit branch ids are removed.
    $tree_text =~ s/([A-Za-z0-9_.-]+):B\d+/$1/g;

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
    return 0 unless defined $value;
    return 0 if $value =~ /^\s*$/;
    return 0 if $value =~ /^NA/;
    return $value =~ /^-?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;
}
