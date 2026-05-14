#!/usr/bin/env perl
# ==============================================================================
# Script:      extract_na_fuse.pl
# Author:      Jiaqi Wu (wujiaqi@hiroshima-u.ac.jp)
# Affiliation: Graduate School of Integrated Sciences for Life, Hiroshima University, Japan
# Copyright:   (c) 2026 Jiaqi Wu. All rights reserved.
#
# Description:
#   Marks primitive branches (e.g., B12) as NA_fuse when the corresponding fused
#   branch (e.g., B12|B47) carries numeric fused signal for a gene, indicating that
#   the gene-tree signal cannot resolve the species-tree nodes separately.
#
#   This script keeps ONLY primitive branches in the output table, but updates
#   NA values to NA_fuse when supported by a numeric fused signal.
#
# Inputs:
#   -i / --input  <matrix_with_fuse.txt>
#
# Outputs (default):
#   <input_basename>.na_fuse.txt
#   Example: free.matrix_with_fuse.txt -> free.matrix_with_fuse.na_fuse.txt
#
# Usage:
#   perl extract_na_fuse.pl -i free.matrix_with_fuse.txt
#   perl extract_na_fuse.pl -i free.matrix_with_fuse.txt -o free.matrix.na_fuse.txt
# ==============================================================================

use 5.010;
use strict;
use warnings;
use Getopt::Long;
use File::Basename;

my ($input_matrix, $output_matrix);

GetOptions(
    'i|input=s'  => \$input_matrix,
    'o|output=s' => \$output_matrix,
) or die "[ERROR] Invalid command line arguments.\n";

die "Usage: $0 -i <matrix_with_fuse.txt> [-o <output.txt>]\n"
    unless defined $input_matrix;

# Default output name: insert ".na_fuse" before final ".txt", otherwise append.
if (!defined $output_matrix) {
    if ($input_matrix =~ /\.txt$/) {
        ($output_matrix = $input_matrix) =~ s/\.txt$/.na_fuse.txt/;
    } else {
        $output_matrix = $input_matrix . ".na_fuse.txt";
    }
}

open(my $IN,  '<', $input_matrix)  or die "[ERROR] Cannot open $input_matrix: $!\n";
open(my $OUT, '>', $output_matrix) or die "[ERROR] Cannot write $output_matrix: $!\n";

my (@header_cols, @primitive_branches, @fused_branches);

my $row_no = 0;
while (my $line = <$IN>) {
    chomp $line;
    $row_no++;
    my @f = split(/\t/, $line);

    # Header
    if ($row_no == 1) {
        @header_cols = @f;

        # Identify primitive vs fused branches.
        for my $col (@header_cols[1 .. $#header_cols]) {
            if ($col =~ /^B\d+(?:\|B\d+)+$/) {
                push @fused_branches, $col;
            } elsif ($col =~ /^B\d+$/) {
                push @primitive_branches, $col;
            }
        }

        # Write header: keep row ID column + primitive branches only
        print {$OUT} $header_cols[0];
        for my $b (@primitive_branches) {
            print {$OUT} "\t$b";
        }
        print {$OUT} "\n";
        next;
    }

    # Data row
    my $gene_id = $f[0];

    # Map column -> value for this row
    my %val_for;
    for my $i (1 .. $#header_cols) {
        $val_for{$header_cols[$i]} = $f[$i];
    }

    # If a fused branch has a numeric signal, mark any NA primitive component as NA_fuse
    for my $fb (@fused_branches) {
        my $fb_val = $val_for{$fb};

        next unless is_numeric_branch_value($fb_val);

        my @parts = split(/\|/, $fb);
        for my $p (@parts) {
            next unless exists $val_for{$p};
            if (defined $val_for{$p} && $val_for{$p} =~ /^NA/) {
                $val_for{$p} = 'NA_fuse';
            }
        }
    }

    # Write primitive-only output row
    print {$OUT} $gene_id;
    for my $b (@primitive_branches) {
        my $v = exists $val_for{$b} ? $val_for{$b} : 'NA';
        print {$OUT} "\t$v";
    }
    print {$OUT} "\n";
}

close $IN;
close $OUT;

print STDERR "[INFO] Wrote: $output_matrix\n";

sub is_numeric_branch_value {
    my ($value) = @_;
    return 0 unless defined $value;
    return 0 if $value =~ /^\s*$/;
    return 0 if $value =~ /^NA/;
    return $value =~ /^-?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;
}
