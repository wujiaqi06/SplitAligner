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
#
# Outputs (prefix = -o):
#   <prefix>.fix.na_classified.txt   : fixed matrix with NA_struct applied
#   <prefix>.free.na_classified.txt  : free matrix with NA_struct / NA_topo applied
#
# Usage:
#   perl confirm_na_structure.pl --fix fix.matrix.txt --free free.matrix.txt -o out_prefix
# ==============================================================================

use strict;
use warnings;
use Getopt::Long;

my ($fix_matrix_path, $free_matrix_path, $out_prefix);

GetOptions(
    'fix=s'  => \$fix_matrix_path,
    'free=s' => \$free_matrix_path,
    'o=s'    => \$out_prefix,
) or die "[ERROR] Invalid command line arguments.\n";

die "Usage: $0 --fix <fix_matrix> --free <free_matrix> -o <out_prefix>\n"
    unless defined $fix_matrix_path && defined $free_matrix_path && defined $out_prefix;

my $out_fix  = "$out_prefix.fix.na_classified.txt";
my $out_free = "$out_prefix.free.na_classified.txt";

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
my $fix_axis  = join("_", sort @fix_branches);
my $free_axis = join("_", sort @free_branches);

if ($fix_axis eq $free_axis) {
    print STDERR "[INFO] Branch axes matched between fix and free matrices.\n";
} else {
    print STDERR "[WARN] Branch axes DO NOT match between fix and free matrices.\n";
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

# -------------------------
# Open outputs only after validation
# -------------------------
open(my $OUT_FIX,  '>', $out_fix)  or die "[ERROR] Cannot write $out_fix: $!\n";
open(my $OUT_FREE, '>', $out_free) or die "[ERROR] Cannot write $out_free: $!\n";

print {$OUT_FIX}  "$fix_header_line\n";
print {$OUT_FREE} "$free_header_line\n";

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
        elsif ($fix{$b} ne 'NA' && $free{$b} eq 'NA') {
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
