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
    write_utf8_file
);

decode_argv_utf8();
configure_utf8_stdio();

my ($input_matrix, $output_matrix, $state_matrix);

GetOptions(
    'i|input=s'  => \$input_matrix,
    'o|output=s' => \$output_matrix,
    's|state=s'  => \$state_matrix,
) or die "[ERROR] Invalid command line arguments.\n";

die "Usage: $0 -i <matrix_with_fuse.txt> --state <primitive_state.tsv> [-o <output.txt>]\n"
    unless defined $input_matrix && defined $state_matrix;

# Default output name: insert ".na_fuse" before final ".txt", otherwise append.
if (!defined $output_matrix) {
    if ($input_matrix =~ /\.txt$/) {
        ($output_matrix = $input_matrix) =~ s/\.txt$/.na_fuse.txt/;
    } else {
        $output_matrix = $input_matrix . ".na_fuse.txt";
    }
}

my $matrix = read_matrix_table($input_matrix);
my $state = read_state_matrix($state_matrix);
validate_state_for_matrix(
    state  => $state,
    matrix => $matrix,
    axis   => $matrix->{primitive},
    role   => 'NA_fuse input',
);

my @primitive_branches = @{$matrix->{primitive}};
my @fused_branches = grep { /\|/ } @{$matrix->{branches}};
my %matrix_index = map { $matrix->{branches}[$_] => $_ } 0 .. $#{$matrix->{branches}};
my %state_index = map { $state->{axis}[$_] => $_ } 0 .. $#{$state->{axis}};
my @output = (join("\t", 'gene', @primitive_branches));

for my $gene_id (@{$matrix->{genes}}) {
    my @values = @{$matrix->{rows}{$gene_id}};
    my @states = @{$state->{rows}{$gene_id}};

    # If a fused branch has a numeric signal, mark any NA primitive component as NA_fuse
    for my $fb (@fused_branches) {
        my $fb_val = $values[ $matrix_index{$fb} ];

        next unless is_numeric_branch_value($fb_val);

        my @parts = split(/\|/, $fb);
        for my $p (@parts) {
            die "[ERROR][gene $gene_id] Fused coordinate '$fb' contains unknown primitive '$p'.\n"
                unless exists $matrix_index{$p} && exists $state_index{$p};
            my $value_index = $matrix_index{$p};
            next unless $values[$value_index] eq 'NA';
            my $code = $states[ $state_index{$p} ];
            die "[ERROR][gene $gene_id branch $p] Numeric fused evidence is incompatible with state '$code'.\n"
                unless $code eq 'F';
            $values[$value_index] = 'NA_fuse';
            }
    }

    my @primitive_values = map { $values[ $matrix_index{$_} ] } @primitive_branches;
    push @output, join("\t", $gene_id, @primitive_values);
}

write_utf8_file($output_matrix, join("\n", @output) . "\n");

print STDERR "[INFO] Wrote: $output_matrix\n";

sub is_numeric_branch_value {
    my ($value) = @_;
    return is_numeric_value($value);
}
