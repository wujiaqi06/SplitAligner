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

my (%branch_to_split, %split_to_branch, %terminal_to_branch);

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

    $branch_to_split{$canonical_split} = $branch_id;
    $split_to_branch{$branch_id} = $canonical_split;

    # terminal split bookkeeping for missing-taxa pruning
    my @parts = split(/\|\|/, $canonical_split);
    if (@parts == 2) {
        if ($parts[0] !~ /\.\./) {
            $terminal_to_branch{$parts[0]} = $branch_id;
        } elsif ($parts[1] !~ /\.\./) {
            $terminal_to_branch{$parts[1]} = $branch_id;
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
    my %projected_branch_to_split = %split_to_branch; # branch_id -> split

    if (@missing_species) {
        # Remove terminal branches that correspond to missing taxa
        for my $sp (@missing_species) {
            my $terminal_branch = $terminal_to_branch{$sp};
            delete $projected_branch_to_split{$terminal_branch} if defined $terminal_branch;
        }

        # Prune missing taxa from every split; drop invalid splits
        for my $branch_id (keys %projected_branch_to_split) {
            my $s = $projected_branch_to_split{$branch_id};
            $s = prune_taxa_from_split($s, \@missing_species);
            my @parts = split(/\|\|/, $s);
            if (@parts != 2 || $parts[0] !~ /\w+/ || $parts[1] !~ /\w+/) {
                delete $projected_branch_to_split{$branch_id};
            } else {
                $projected_branch_to_split{$branch_id} = reorder_split($s);
            }
        }

        # Collect mapping: gene split -> one or multiple projected branch ids (fused)
        my %split_to_projected_branches;
        my ($hit, $miss) = (0, 0);

        for my $branch_id (sort keys %projected_branch_to_split) {
            my $proj_split = $projected_branch_to_split{$branch_id};
            if (exists $gene_split_to_values{$proj_split}) {
                $hit++;
                if (!exists $split_to_projected_branches{$proj_split}) {
                    $split_to_projected_branches{$proj_split} = $branch_id;
                } else {
                    $split_to_projected_branches{$proj_split} .= "|$branch_id";
                }
            } else {
                $miss++;
                print {$ERR} "[MISS] $fname\t$branch_id\t$proj_split\n";
            }
        }

        print {$ERR} "[SUMMARY] $fname\thit=$hit\tmiss=$miss\tmissing_taxa=" . join(',', @missing_species) . "\n";

        for my $gene_split (sort keys %split_to_projected_branches) {
            my $branch_pattern = $split_to_projected_branches{$gene_split};
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
