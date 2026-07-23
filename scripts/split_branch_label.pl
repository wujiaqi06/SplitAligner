#!/usr/bin/env perl
# ==============================================================================
# Script:      split_branch_label.pl
# Author:      Jiaqi Wu (wujiaqi@hiroshima-u.ac.jp, wujiaqi06@gmail.com)
# Affiliation: Graduate School of Integrated Sciences for Life, Hiroshima University, Japan
# Copyright:   (c) 2026 Jiaqi Wu. All rights reserved.
#
# Description:
#   Map bipartite splits extracted from individual gene trees onto the
#   standardized branch coordinate system defined by the species tree.
#   Missing taxa in a gene tree are handled by dynamically pruning the
#   species-tree split set before the 1-to-1 topological matching.
#
# Input / Output conventions (publish-grade):
#   -i : directory produced by tree_to_splits.pl in gene mode: <label>_splits/
#        containing many '*.split.txt'
#   -j : species split coordinate axis: species_tree.splits.txt
#   -o : label (NOT a folder). Output directory is <label>_split_branch_label/
#
# Outputs:
#   <label>_split_branch_label/*.split.txt
#   <label>_split_branch_label/<label>.split_branch_label.errors.log
#   <label>_split_branch_label/<label>.branch_patterns.txt
#   <label>.primitive_state.tsv
#
# Usage:
#   perl split_branch_label.pl -i <label>_splits -j species_tree.splits.txt -o <label>
# ==============================================================================

use 5.010;
use strict;
use warnings;
use Getopt::Std;
use FindBin qw($RealBin);
use File::Spec;
use File::Path qw(make_path);
use lib "$RealBin/../lib";
use SplitAligner::Newick qw(
    canonical_split
    canonicalize_split
    split_taxon_sets
);
use SplitAligner::CoordinateState qw(write_state_matrix);
use SplitAligner::Provenance qw(read_axis_ledger);
use SplitAligner::TextIO qw(
    close_utf8_writer
    configure_utf8_stdio
    decode_argv_utf8
    open_utf8_writer
    print_utf8
    read_utf8_lines
    write_utf8_file
);

decode_argv_utf8();
configure_utf8_stdio();

my %opt;
getopts('i:j:o:a:g:', \%opt);

my $gene_splits_dir = $opt{i} // '';
my $species_axis_file = $opt{j} // '';
my $label = $opt{o} // '';
my $primitive_axis_file = $opt{a} // '';
my $gene_map_file = $opt{g} // '';

if (!$gene_splits_dir || !$species_axis_file || !$label || !$primitive_axis_file || !$gene_map_file) {
    die "Usage: $0 -i <label>_splits -j species_tree.splits.txt -a <primitive_axis.tsv> -g <gene_id_map.tsv> -o <label>\n";
}

-d $gene_splits_dir or die "[ERROR] Cannot find input directory: $gene_splits_dir\n";
-e $species_axis_file or die "[ERROR] Cannot read species split axis file: $species_axis_file\n";
-e $primitive_axis_file or die "[ERROR] Cannot read primitive axis file: $primitive_axis_file\n";
-e $gene_map_file or die "[ERROR] Cannot read gene identity map: $gene_map_file\n";

my $axis = read_axis_ledger($primitive_axis_file);
my @primitive_axis = map { $_->{branch_id} } @{$axis};
my ($storage_to_gene, $gene_order) = read_gene_id_map($gene_map_file);

my $out_dir = "${label}_split_branch_label";
make_path($out_dir);

my $error_log = File::Spec->catfile($out_dir, "${label}.split_branch_label.errors.log");
my $pattern_out = File::Spec->catfile($out_dir, "${label}.branch_patterns.txt");
my $state_out = "${label}.primitive_state.tsv";

my $ERR = open_utf8_writer($error_log);

# -------------------------------
# Read species-tree split axis
# -------------------------------
my (%branch_to_split, %split_to_branch, %terminal_to_branch, %is_terminal_branch);

for my $line (@{ read_utf8_lines($species_axis_file) }) {
    next if $line eq '';

    my @fields = split(/\t/, $line);
    next unless @fields >= 2;

    my $raw_split = $fields[0];
    my $branch_id = $fields[1];

    my $canonical_split = reorder_split($raw_split);
    if ($canonical_split eq 'Error') {
        print_utf8($ERR, "[ERROR] Bad split line in species axis: $line\n");
        next;
    }

    if (exists $branch_to_split{$canonical_split}) {
        my $current_branch = $branch_to_split{$canonical_split};
        my $winner = preferred_branch_id($current_branch, $branch_id);
        my $loser  = ($winner eq $current_branch) ? $branch_id : $current_branch;

        print_utf8($ERR, "[WARN] Duplicate canonical species split detected: $canonical_split\twinner=$winner\tloser=$loser\n");

        if ($winner eq $current_branch) {
            next;
        }

        delete $split_to_branch{$current_branch};
        delete $is_terminal_branch{$current_branch};
        for my $tip (keys %terminal_to_branch) {
            delete $terminal_to_branch{$tip} if $terminal_to_branch{$tip} eq $current_branch;
        }
    }

    $branch_to_split{$canonical_split} = $branch_id;
    $split_to_branch{$branch_id} = $canonical_split;

    # terminal split bookkeeping for missing-taxa pruning
    my ($left_taxa, $right_taxa) = split_taxon_sets($canonical_split);
    if (@{$left_taxa} == 1) {
            $terminal_to_branch{$left_taxa->[0]} = $branch_id;
            $is_terminal_branch{$branch_id} = 1;
    } elsif (@{$right_taxa} == 1) {
            $terminal_to_branch{$right_taxa->[0]} = $branch_id;
            $is_terminal_branch{$branch_id} = 1;
    }
}

die "[ERROR] species split axis seems empty: $species_axis_file\n" unless %branch_to_split;
my %axis_branch = map { $_ => 1 } @primitive_axis;
die "[ERROR] Species split axis and primitive ledger contain different retained branch IDs.\n"
    unless keys(%axis_branch) == keys(%split_to_branch)
        && !grep { !$split_to_branch{$_} } @primitive_axis;

# Derive the full species set from one complete split.
my @axis_splits = sort keys %branch_to_split;
my $first_split = $axis_splits[0];
my ($first_left, $first_right) = split_taxon_sets($first_split);
my @all_species = sort(@{$first_left}, @{$first_right});

# -------------------------------
# Iterate over gene split files
# -------------------------------
opendir(my $DH, $gene_splits_dir) or die "[ERROR] Cannot open directory $gene_splits_dir: $!\n";
my @split_files = sort grep { /\.split\.txt$/ && -f File::Spec->catfile($gene_splits_dir, $_) } readdir($DH);
closedir $DH;

my %observed_branch_patterns;  # for summary output
my %state_for_gene;
my %seen_storage;

for my $fname (@split_files) {
    my $in_path  = File::Spec->catfile($gene_splits_dir, $fname);
    my $out_path = File::Spec->catfile($out_dir, $fname);

    (my $storage_key = $fname) =~ s/\.split\.txt\z//;
    die "[ERROR] Split file has no gene identity mapping: $fname\n"
        unless exists $storage_to_gene->{$storage_key};
    die "[ERROR] Duplicate split file for storage key '$storage_key'.\n"
        if $seen_storage{$storage_key}++;
    my $gene_id = $storage_to_gene->{$storage_key};
    my $OUT = open_utf8_writer($out_path);

    my %gene_split_to_values;

    for my $line (@{ read_utf8_lines($in_path) }) {
        next if $line eq '';
        my @fields = split(/\t/, $line);
        next unless @fields >= 2;

        my $raw_split = $fields[0];
        my $val = $fields[1];

        my $canonical_split = reorder_split($raw_split);
        if ($canonical_split eq 'Error') {
            print_utf8($ERR, "[ERROR] Bad split in $fname: $line\n");
            next;
        }

        push @{ $gene_split_to_values{$canonical_split} }, $val;
    }

    my @gene_splits = sort keys %gene_split_to_values;
    if (!@gene_splits) {
        print_utf8($ERR, "[WARN] Empty gene split file: $fname\n");
        close_utf8_writer($OUT, $out_path);
        next;
    }

    # Derive species present in this gene from one complete split.
    my $g0 = $gene_splits[0];
    my ($gene_left, $gene_right) = split_taxon_sets($g0);
    my @gene_species = sort(@{$gene_left}, @{$gene_right});

    my @missing_species = array_split(\@all_species, \@gene_species);
    my %state = map { $_ => 'U' } @primitive_axis;

    # Build projected species axis under missing taxa
    my %projected_branch_meta = map {
        $_ => {
            split      => $split_to_branch{$_},
            observable => 0,
            fuse_ok    => 0,
        }
    } keys %split_to_branch;

    if (@missing_species) {
        # Remove terminal branches that correspond to missing taxa
        for my $sp (@missing_species) {
            my $terminal_branch = $terminal_to_branch{$sp};
            if (defined $terminal_branch) {
                $state{$terminal_branch} = 'S';
                delete $projected_branch_meta{$terminal_branch};
            }
        }

        # Prune missing taxa from every split and classify each projected branch
        # into:
        #   observable => can retain an independent numeric value
        #   fuse_ok    => can still participate in fused-path bookkeeping
        for my $branch_id (keys %projected_branch_meta) {
            my $s = $projected_branch_meta{$branch_id}{split};
            $s = prune_taxa_from_split($s, \@missing_species);
            my ($observable, $fuse_ok) = classify_projected_split($s, $is_terminal_branch{$branch_id});
            if (!$fuse_ok) {
                $state{$branch_id} = 'S';
                delete $projected_branch_meta{$branch_id};
                next;
            }

            $projected_branch_meta{$branch_id}{split}      = reorder_split($s);
            $projected_branch_meta{$branch_id}{observable} = $observable;
            $projected_branch_meta{$branch_id}{fuse_ok}    = $fuse_ok;
        }

        # Collect mapping: gene split -> one or multiple projected branch ids.
        # Endpoint-collapsed internal branches (e.g. 1|k) remain eligible for
        # fuse bookkeeping, but they should not be emitted as singleton numeric
        # branches when they are not independently observable.
        my %split_to_projected_branches;
        my ($hit, $miss) = (0, 0);

        for my $branch_id (sort keys %projected_branch_meta) {
            my $proj_split = $projected_branch_meta{$branch_id}{split};
            if (exists $gene_split_to_values{$proj_split}) {
                $hit++;
                push @{ $split_to_projected_branches{$proj_split} }, $branch_id;
            } else {
                $miss++;
                print_utf8($ERR, "[MISS] $fname\t$branch_id\t$proj_split\n");
            }
        }

        print_utf8($ERR, "[SUMMARY] $fname\thit=$hit\tmiss=$miss\tmissing_taxa=" . join(',', @missing_species) . "\n");

        for my $gene_split (sort keys %split_to_projected_branches) {
            my @members = @{ $split_to_projected_branches{$gene_split} };
            my $has_observable = 0;
            for my $branch_id (@members) {
                if ($projected_branch_meta{$branch_id}{observable}) {
                    $has_observable = 1;
                    last;
                }
            }

            # Suppress singleton endpoint-collapsed internal branches. They are
            # valid bookkeeping objects for fused-path construction, but they
            # must not survive as independently observed numeric branches.
            if (@members >= 2) {
                $state{$_} = 'F' for @members;
            } elsif ($has_observable) {
                $state{$members[0]} = 'D';
            }

            next if @members == 1 && !$has_observable;

            my $branch_pattern = join('|', sort_branch_ids(@members));
            for my $val (@{ $gene_split_to_values{$gene_split} }) {
                print_utf8($OUT, "$gene_split\t$branch_pattern\t$val\n");
                $observed_branch_patterns{$branch_pattern} = 1;
            }
        }

    } else {
        # No missing taxa: direct match against full axis (canonical splits)
        for my $axis_split (sort keys %branch_to_split) {
            if (exists $gene_split_to_values{$axis_split}) {
                my $branch_id = $branch_to_split{$axis_split};
                $state{$branch_id} = 'D';
                for my $val (@{ $gene_split_to_values{$axis_split} }) {
                    print_utf8($OUT, "$axis_split\t$branch_id\t$val\n");
                    $observed_branch_patterns{$branch_id} = 1;
                }
            } else {
                # keep silent (original script printed to STDOUT); logging would be noisy for large data
            }
        }
    }

    $state_for_gene{$gene_id} = [ map { $state{$_} } @primitive_axis ];
    close_utf8_writer($OUT, $out_path);
}

# Write branch pattern summary
write_utf8_file($pattern_out, join('', map { "$_\n" } sort keys %observed_branch_patterns));
close_utf8_writer($ERR, $error_log);

for my $storage_key (keys %{$storage_to_gene}) {
    die "[ERROR] Missing split file for gene '$storage_to_gene->{$storage_key}' (storage key $storage_key).\n"
        unless $seen_storage{$storage_key};
}
my @state_rows = map {
    my $gene_id = $_;
    die "[ERROR] Missing coordinate state row for gene '$gene_id'.\n"
        unless exists $state_for_gene{$gene_id};
    { gene_id => $gene_id, states => $state_for_gene{$gene_id} }
} @{$gene_order};
write_state_matrix($state_out, \@primitive_axis, \@state_rows);

# -------------------------------
# Utilities
# -------------------------------
sub array_split {
    my ($universe, $subset) = @_;
    my %present = map { $_ => 1 } @{$subset};
    return grep { !$present{$_} } @{$universe};
}

sub reorder_split {
    my ($split) = @_;
    my $canonical = eval { canonicalize_split($split) };
    return $@ ? 'Error' : $canonical;
}

sub prune_taxa_from_split {
    my ($split, $missing_ref) = @_;
    my %missing = map { $_ => 1 } @{$missing_ref};

    my ($left_taxa, $right_taxa) = split_taxon_sets($split);
    my @left = grep { !$missing{$_} } @{$left_taxa};
    my @right = grep { !$missing{$_} } @{$right_taxa};
    return canonical_split(\@left, \@right);
}

sub classify_projected_split {
    my ($split, $is_terminal) = @_;

    my ($left_taxa, $right_taxa) = split_taxon_sets($split);
    my @left  = @{$left_taxa};
    my @right = @{$right_taxa};

    # Structural absence: one projected side is empty.
    return (0, 0) if !@left || !@right;

    # Terminal branches remain independently observable when the terminal taxon
    # is present and the opposite side is non-empty.
    return (1, 1) if $is_terminal;

    # Internal branches are independently observable only when the projected
    # unrooted split remains non-trivial.
    return (1, 1) if @left >= 2 && @right >= 2;

    # Endpoint-collapsed internal branches (e.g. 1|k or 1|1) are not
    # independently observable, but they still participate in fused-path
    # bookkeeping and must not be discarded as structural absence.
    return (0, 1);
}

sub sort_branch_ids {
    return sort {
        my ($an) = $a =~ /^B(\d+)$/;
        my ($bn) = $b =~ /^B(\d+)$/;
        return $an <=> $bn;
    } @_;
}

sub preferred_branch_id {
    my ($left, $right) = @_;
    my ($ln) = defined $left  ? ($left  =~ /^B(\d+)$/) : ();
    my ($rn) = defined $right ? ($right =~ /^B(\d+)$/) : ();

    return $left  if defined $ln && defined $rn && $ln <= $rn;
    return $right if defined $ln && defined $rn;
    return defined $left ? $left : $right;
}

sub read_gene_id_map {
    my ($path) = @_;
    my $lines = read_utf8_lines($path);
    die "[ERROR] Gene identity map is empty: $path\n" unless @{$lines};
    die "[ERROR] Unexpected gene identity map header in $path.\n"
        unless $lines->[0] eq "storage_key\tgene_id\tinput_line";

    my (%storage_to_gene, %seen_gene, @genes);
    for my $index (1 .. $#{$lines}) {
        my $line = $lines->[$index];
        next if $line eq '';
        my ($storage_key, $gene_id, $input_line) = split(/\t/, $line, -1);
        die "[ERROR][$path line " . ($index + 1) . "] Malformed gene identity row.\n"
            unless defined $input_line && defined $storage_key && $storage_key ne ''
                && defined $gene_id && $gene_id ne '';
        die "[ERROR][$path line " . ($index + 1) . "] Duplicate storage key '$storage_key'.\n"
            if exists $storage_to_gene{$storage_key};
        die "[ERROR][$path line " . ($index + 1) . "] Duplicate gene ID '$gene_id'.\n"
            if $seen_gene{$gene_id}++;
        $storage_to_gene{$storage_key} = $gene_id;
        push @genes, $gene_id;
    }
    die "[ERROR] Gene identity map has no records: $path\n" unless @genes;
    return (\%storage_to_gene, \@genes);
}
