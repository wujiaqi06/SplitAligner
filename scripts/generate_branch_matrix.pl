#!/usr/bin/env perl
# ==============================================================================
# Script:      generate_branch_matrix.pl
# Author:      Jiaqi Wu (wujiaqi@hiroshima-u.ac.jp, wujiaqi06@gmail.com)
# Affiliation: Graduate School of Integrated Sciences for Life, Hiroshima University, Japan
#
# Description:
#   Generate gene x branch matrices from mapped splits (output of split_branch_label.pl).
#   This script writes TWO matrices:
#     1) <label>.matrix_no_fuse.txt   : primitive branches only (canonical axis)
#     2) <label>.matrix_with_fuse.txt : primitive axis + fused branches (e.g., B12|B47)
#
# Usage:
#   perl generate_branch_matrix.pl -i <label>_split_branch_label -o <label>
#
# Semantics:
#   -i : input directory created by split_branch_label.pl (<label>_split_branch_label/)
#   -o : label (NOT a folder). Output files are written to current directory.
# ==============================================================================

use 5.010;
use strict;
use warnings;
use Getopt::Std;
use File::Spec;

my %opts;
getopts('i:o:', \%opts);

my $in_dir = $opts{'i'} or die "Usage: $0 -i <label>_split_branch_label -o <label>\n";
my $label  = $opts{'o'} or die "Usage: $0 -i <label>_split_branch_label -o <label>\n";

die "[ERROR] Input directory not found: $in_dir\n" unless -d $in_dir;

# Collect input split files
opendir(my $DH, $in_dir) or die "[ERROR] Cannot open directory $in_dir: $!\n";
my @files = sort grep { /\.split\.txt$/ && -f File::Spec->catfile($in_dir, $_) } readdir($DH);
closedir($DH);

# gene_id -> { pattern => value }
my %gene_branch;
my %genes;

# Track primitive branch axis (B1..Bmax)
my $max_b = 0;

# Collect fused patterns (contain '|')
my %fused_patterns;

foreach my $fn (@files) {
    my $path = File::Spec->catfile($in_dir, $fn);
    (my $gene_id = $fn) =~ s/\.split\.txt$//;
    $genes{$gene_id} = 1;

    open(my $IN, '<', $path) or die "[ERROR] Cannot read $path: $!\n";
    while (my $line = <$IN>) {
        chomp $line;
        next if $line =~ /^\s*$/;

        my @f = split(/\t/, $line);
        next unless @f >= 3; # expected: split \t branch_pattern \t value

        my $pattern = $f[1];
        my $value   = $f[2];

        # Canonicalize fused order: B3|B1 -> B1|B3
        if ($pattern =~ /\|/) {
            my @parts = sort_branch_ids(split(/\|/, $pattern));
            $pattern = join('|', @parts);
            $fused_patterns{$pattern} = 1;
        }

        if (exists $gene_branch{$gene_id}{$pattern}) {
            $gene_branch{$gene_id}{$pattern} = _merge_pattern_values(
                $gene_branch{$gene_id}{$pattern},
                $value,
            );
        } else {
            $gene_branch{$gene_id}{$pattern} = $value;
        }

        # Update primitive max branch id
        for my $b (split(/\|/, $pattern)) {
            if ($b =~ /^B(\d+)$/) {
                $max_b = $1 if $1 > $max_b;
            }
        }
    }
    close $IN;
}

# Canonical primitive axis
my @primitive = map { "B$_" } (1 .. $max_b);

# Sort fused patterns by their numeric components (stable, human-friendly)
my @fused_sorted = sort {
    _cmp_branch_pattern($a, $b)
} keys %fused_patterns;

# Output files
my $out_no   = "${label}.matrix_no_fuse.txt";
my $out_with = "${label}.matrix_with_fuse.txt";

open(my $NO,   '>', $out_no)   or die "[ERROR] Cannot write $out_no: $!\n";
open(my $WITH, '>', $out_with) or die "[ERROR] Cannot write $out_with: $!\n";

# Headers
print {$NO}   join("\t", "gene", @primitive), "\n";
print {$WITH} join("\t", "gene", @primitive, @fused_sorted), "\n";

# Rows
for my $gene_id (sort keys %genes) {

    # no_fuse
    my @row_no = ($gene_id);
    for my $b (@primitive) {
        if (exists $gene_branch{$gene_id}{$b}) {
            push @row_no, $gene_branch{$gene_id}{$b};
        } else {
            push @row_no, "NA";
        }
    }
    print {$NO} join("\t", @row_no), "\n";

    # with_fuse
    my @row_with = @row_no; # starts with gene + primitive axis
    for my $pat (@fused_sorted) {
        if (exists $gene_branch{$gene_id}{$pat}) {
            push @row_with, $gene_branch{$gene_id}{$pat};
            next;
        }

        # Otherwise, attempt to synthesize from primitive components.
        my $synth = _synthesize_fused_value($gene_branch{$gene_id}, $pat);
        push @row_with, $synth;
    }
    print {$WITH} join("\t", @row_with), "\n";
}

close $NO;
close $WITH;

exit 0;

# ------------------------------------------------------------------------------
# Internal helpers
# ------------------------------------------------------------------------------

sub _synthesize_fused_value {
    my ($gene_href, $pattern) = @_;
    my @parts = split(/\|/, $pattern);

    my $sum = 0;
    for my $b (@parts) {
        return "NA" unless exists $gene_href->{$b};
        my $v = $gene_href->{$b};

        # Only accept numeric values; otherwise NA (avoid implicit 0)
        return "NA" unless is_numeric_branch_value($v);
        $sum += $v;
    }
    return $sum;
}

sub _merge_pattern_values {
    my ($old, $new) = @_;

    my $is_old_num = is_numeric_branch_value($old);
    my $is_new_num = is_numeric_branch_value($new);

    return $old + $new if $is_old_num && $is_new_num;
    return $old if defined $old && (!defined $new || $new eq '');
    return $new;
}

sub is_numeric_branch_value {
    my ($value) = @_;
    return 0 unless defined $value;
    return 0 if $value =~ /^\s*$/;
    return 0 if $value =~ /^NA/;
    return $value =~ /^-?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;
}

sub _cmp_branch_pattern {
    my ($a, $b) = @_;

    my @aa = map { /^B(\d+)$/ ? $1 : $_ } split(/\|/, $a);
    my @bb = map { /^B(\d+)$/ ? $1 : $_ } split(/\|/, $b);

    my $n = @aa < @bb ? @aa : @bb;
    for (my $i = 0; $i < $n; $i++) {
        my $x = $aa[$i];
        my $y = $bb[$i];
        if ($x =~ /^\d+$/ && $y =~ /^\d+$/) {
            return $x <=> $y if $x != $y;
        } else {
            my $c = "$x" cmp "$y";
            return $c if $c != 0;
        }
    }
    return @aa <=> @bb;
}

sub sort_branch_ids {
    return sort {
        my ($an) = $a =~ /^B(\d+)$/;
        my ($bn) = $b =~ /^B(\d+)$/;
        return $a cmp $b unless defined $an && defined $bn;
        return $an <=> $bn;
    } @_;
}
