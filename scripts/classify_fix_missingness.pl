#!/usr/bin/env perl
# ==============================================================================
# Script:      classify_fix_missingness.pl
# Author:      Jiaqi Wu (wujiaqi@hiroshima-u.ac.jp)
# Description:
#   Finalize a fixed-topology matrix after NA_fuse marking using explicit
#   per-cell coordinate state. Only STRUCT_ABSENT cells become NA_struct.
#
# Inputs:
#   --fix / -i     <fix.matrix_with_fuse.na_fuse.txt>
#   --output / -o  <final_fix.fix.na_classified.txt>
#
# Output:
#   A primitive-branch matrix where structural cells are NA_struct and
#   nonstructural unavailable values remain residual NA.
# ==============================================================================

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use FindBin qw($RealBin);
use lib "$RealBin/../lib";
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

my ($fix_in, $state_in, $out);

GetOptions(
    'fix|i=s'    => \$fix_in,
    'state|s=s'  => \$state_in,
    'output|o=s' => \$out,
) or die usage();

die usage() unless defined $fix_in && $fix_in ne '';
die "[ERROR] Cannot find input file: $fix_in\n" unless -e $fix_in;
die usage() unless defined $state_in && $state_in ne '';
die "[ERROR] Cannot find state file: $state_in\n" unless -e $state_in;
die usage() unless defined $out && $out ne '';

my $matrix = read_matrix_table($fix_in);
my $state = read_state_matrix($state_in);
validate_state_for_matrix(
    state          => $state,
    matrix         => $matrix,
    axis           => $matrix->{primitive},
    role           => 'fix-only classification input',
    allow_na_fuse  => 1,
);

my %matrix_index = map { $matrix->{branches}[$_] => $_ } 0 .. $#{$matrix->{branches}};
my %state_index = map { $state->{axis}[$_] => $_ } 0 .. $#{$state->{axis}};
my @output = (join("\t", 'gene', @{$matrix->{primitive}}));

for my $gene (@{$matrix->{genes}}) {
    my @values = @{$matrix->{rows}{$gene}};
    my @states = @{$state->{rows}{$gene}};
    my @classified;
    for my $branch (@{$matrix->{primitive}}) {
        my $value = $values[ $matrix_index{$branch} ];
        my $code = $states[ $state_index{$branch} ];
        $value = 'NA_struct' if $value eq 'NA' && $code eq 'S';
        push @classified, $value;
    }
    push @output, join("\t", $gene, @classified);
}

write_utf8_file($out, join("\n", @output) . "\n");

print STDERR "[INFO] Wrote: $out\n";
exit 0;

sub usage {
    return <<"USAGE";
Usage:
  perl $0 --fix <fix.matrix_with_fuse.na_fuse.txt> --state <fix.primitive_state.tsv> --output <final_fix.fix.na_classified.txt>
USAGE
}
