#!/usr/bin/env perl
# ==============================================================================
# Script:      classify_fix_missingness.pl
# Author:      Jiaqi Wu (wujiaqi@hiroshima-u.ac.jp)
# Description:
#   Finalize a fixed-topology matrix after NA_fuse marking. In the fix-only
#   setting, any remaining generic NA is interpreted as NA_struct.
#
# Inputs:
#   --fix / -i     <fix.matrix_with_fuse.na_fuse.txt>
#   --output / -o  <final_fix.fix.na_classified.txt>
#
# Output:
#   A primitive-branch matrix where remaining NA cells are rewritten as
#   NA_struct and existing numeric / NA_fuse values are preserved.
# ==============================================================================

use strict;
use warnings;
use Getopt::Long qw(GetOptions);

my ($fix_in, $out);

GetOptions(
    'fix|i=s'    => \$fix_in,
    'output|o=s' => \$out,
) or die usage();

die usage() unless defined $fix_in && $fix_in ne '';
die "[ERROR] Cannot find input file: $fix_in\n" unless -e $fix_in;
die usage() unless defined $out && $out ne '';

open(my $IN, '<', $fix_in) or die "[ERROR] Cannot open $fix_in: $!\n";
open(my $OUT, '>', $out) or die "[ERROR] Cannot write $out: $!\n";

my $row_no = 0;
while (my $line = <$IN>) {
    chomp $line;
    $row_no++;

    my @fields = split(/\t/, $line, -1);
    if ($row_no == 1) {
        print {$OUT} join("\t", @fields) . "\n";
        next;
    }

    for my $i (1 .. $#fields) {
        next unless defined $fields[$i];
        if ($fields[$i] eq 'NA') {
            $fields[$i] = 'NA_struct';
        }
    }

    print {$OUT} join("\t", @fields) . "\n";
}

close $IN;
close $OUT;

print STDERR "[INFO] Wrote: $out\n";
exit 0;

sub usage {
    return <<"USAGE";
Usage:
  perl $0 --fix <fix.matrix_with_fuse.na_fuse.txt> --output <final_fix.fix.na_classified.txt>
USAGE
}
