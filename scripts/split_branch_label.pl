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
#
# Usage:
#   perl split_branch_label.pl -i <label>_splits -j species_tree.splits.txt -o <label>
# ==============================================================================

use 5.010;
use strict;
use warnings;
use Getopt::Std;
use File::Spec;
use File::Path qw(make_path);

my %opt;
getopts('i:j:o:', \%opt);

my $gene_splits_dir = $opt{i} // '';
my $species_axis_file = $opt{j} // '';
my $label = $opt{o} // '';

if (!$gene_splits_dir || !$species_axis_file || !$label) {
    die "Usage: $0 -i <label>_splits -j species_tree.splits.txt -o <label>\n";
}

-d $gene_splits_dir or die "[ERROR] Cannot find input directory: $gene_splits_dir\n";
-e $species_axis_file or die "[ERROR] Cannot read species split axis file: $species_axis_file\n";

my $out_dir = "${label}_split_branch_label";
make_path($out_dir);

my $error_log = File::Spec->catfile($out_dir, "${label}.split_branch_label.errors.log");
my $pattern_out = File::Spec->catfile($out_dir, "${label}.branch_patterns.txt");

open(my $ERR, '>', $error_log) or die "[ERROR] Cannot write $error_log: $!\n";

# -------------------------------
# Read species-tree split axis
# -------------------------------
open(my $AXIS, '<', $species_axis_file) or die "[ERROR] Cannot read $species_axis_file: $!\n";

my (%branch_to_split, %split_to_branch, %terminal_to_branch, %is_terminal_branch);

while (my $line = <$AXIS>) {
    chomp $line;
    next if $line eq '';

    my @fields = split(/\t/, $line);
    next unless @fields >= 2;

    my $raw_split = $fields[0];
    my $branch_id = $fields[1];

    my $canonical_split = reorder_split($raw_split);
    if ($canonical_split eq 'Error') {
        print {$ERR} "[ERROR] Bad split line in species axis: $line\n";
        next;
    }

    if (exists $branch_to_split{$canonical_split}) {
        my $current_branch = $branch_to_split{$canonical_split};
        my $winner = preferred_branch_id($current_branch, $branch_id);
        my $loser  = ($winner eq $current_branch) ? $branch_id : $current_branch;

        print {$ERR} "[WARN] Duplicate canonical species split detected: $canonical_split\twinner=$winner\tloser=$loser\n";

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
    my @parts = split(/\|\|/, $canonical_split);
    if (@parts == 2) {
        if ($parts[0] !~ /\.\./) {
            $terminal_to_branch{$parts[0]} = $branch_id;
            $is_terminal_branch{$branch_id} = 1;
        } elsif ($parts[1] !~ /\.\./) {
            $terminal_to_branch{$parts[1]} = $branch_id;
            $is_terminal_branch{$branch_id} = 1;
        }
    }
}
close $AXIS;

die "[ERROR] species split axis seems empty: $species_axis_file\n" unless %branch_to_split;

# Derive full species set string (..sp1..sp2..)
my @axis_splits = sort keys %branch_to_split;
my $first_split = $axis_splits[0];
(my $tmp = $first_split) =~ s/\|\|/../;
my @all_species = split(/\.\./, $tmp);
my $species_universe = '..' . join('..', @all_species) . '..';

# -------------------------------
# Iterate over gene split files
# -------------------------------
opendir(my $DH, $gene_splits_dir) or die "[ERROR] Cannot open directory $gene_splits_dir: $!\n";
my @split_files = sort grep { /\.split\.txt$/ && -f File::Spec->catfile($gene_splits_dir, $_) } readdir($DH);
closedir $DH;

my %observed_branch_patterns;  # for summary output

for my $fname (@split_files) {
    my $in_path  = File::Spec->catfile($gene_splits_dir, $fname);
    my $out_path = File::Spec->catfile($out_dir, $fname);

    open(my $IN,  '<', $in_path)  or die "[ERROR] Cannot read $in_path: $!\n";
    open(my $OUT, '>', $out_path) or die "[ERROR] Cannot write $out_path: $!\n";

    my %gene_split_to_values;

    while (my $line = <$IN>) {
        chomp $line;
        next if $line eq '';
        my @fields = split(/\t/, $line);
        next unless @fields >= 2;

        my $raw_split = $fields[0];
        my $val = $fields[1];

        my $canonical_split = reorder_split($raw_split);
        if ($canonical_split eq 'Error') {
            print {$ERR} "[ERROR] Bad split in $fname: $line\n";
            next;
        }

        push @{ $gene_split_to_values{$canonical_split} }, $val;
    }

    my @gene_splits = sort keys %gene_split_to_values;
    if (!@gene_splits) {
        print {$ERR} "[WARN] Empty gene split file: $fname\n";
        close $IN;
        close $OUT;
        next;
    }

    # derive species present in this gene (from first split)
    my $g0 = $gene_splits[0];
    (my $g0tmp = $g0) =~ s/\|\|/../;
    my @gene_species = split(/\.\./, $g0tmp);
    my $gene_species_str = '..' . join('..', @gene_species) . '..';

    my @missing_species = array_split($species_universe, $gene_species_str);

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
            delete $projected_branch_meta{$terminal_branch} if defined $terminal_branch;
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
                print {$ERR} "[MISS] $fname\t$branch_id\t$proj_split\n";
            }
        }

        print {$ERR} "[SUMMARY] $fname\thit=$hit\tmiss=$miss\tmissing_taxa=" . join(',', @missing_species) . "\n";

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
            next if @members == 1 && !$has_observable;

            my $branch_pattern = join('|', sort_branch_ids(@members));
            for my $val (@{ $gene_split_to_values{$gene_split} }) {
                print {$OUT} "$gene_split\t$branch_pattern\t$val\n";
                $observed_branch_patterns{$branch_pattern} = 1;
            }
        }

    } else {
        # No missing taxa: direct match against full axis (canonical splits)
        for my $axis_split (sort keys %branch_to_split) {
            if (exists $gene_split_to_values{$axis_split}) {
                my $branch_id = $branch_to_split{$axis_split};
                for my $val (@{ $gene_split_to_values{$axis_split} }) {
                    print {$OUT} "$axis_split\t$branch_id\t$val\n";
                    $observed_branch_patterns{$branch_id} = 1;
                }
            } else {
                # keep silent (original script printed to STDOUT); logging would be noisy for large data
            }
        }
    }

    close $IN;
    close $OUT;
}

# Write branch pattern summary
open(my $PAT, '>', $pattern_out) or die "[ERROR] Cannot write $pattern_out: $!\n";
print {$PAT} "$_\n" for sort keys %observed_branch_patterns;
close $PAT;

close $ERR;

# -------------------------------
# Utilities
# -------------------------------
sub array_split {
    my ($universe, $subset) = @_;
    my @universe_species = split(/\.\./, $universe);
    my @missing;
    for my $sp (@universe_species) {
        next if $sp eq '';
        if ($subset !~ /\.\Q$sp\E\./) {
            push @missing, $sp;
        }
    }
    return @missing;
}

sub reorder_split {
    my ($split) = @_;
    $split =~ s/^\.\.//;
    $split =~ s/\.\.$//;

    my @parts = split(/\|\|/, $split);
    if (@parts != 2) {
        return 'Error';
    }

    my @left  = sort split(/\.\./, $parts[0]);
    my @right = sort split(/\.\./, $parts[1]);

    my $left  = join('..', grep { $_ ne '' } @left);
    my $right = join('..', grep { $_ ne '' } @right);

    return ($left le $right) ? "$left||$right" : "$right||$left";
}

sub prune_taxa_from_split {
    my ($split, $missing_ref) = @_;
    my %missing = map { $_ => 1 } @{$missing_ref};

    my @parts = split(/\|\|/, $split, -1);
    return $split unless @parts == 2;

    my @left = grep { $_ ne '' && !$missing{$_} } split(/\.\./, $parts[0], -1);
    my @right = grep { $_ ne '' && !$missing{$_} } split(/\.\./, $parts[1], -1);

    return join('||', join('..', @left), join('..', @right));
}

sub classify_projected_split {
    my ($split, $is_terminal) = @_;

    my @parts = split(/\|\|/, $split, -1);
    return (0, 0) unless @parts == 2;

    my @left  = grep { $_ ne '' } split(/\.\./, $parts[0], -1);
    my @right = grep { $_ ne '' } split(/\.\./, $parts[1], -1);

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
