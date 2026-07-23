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
#     -a <primitive_axis.tsv> -g <gene_id_map.tsv>
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
use FindBin qw($RealBin);
use lib "$RealBin/../lib";
use SplitAligner::Newick qw(is_numeric_value);
use SplitAligner::Provenance qw(read_axis_ledger);
use SplitAligner::TextIO qw(
    close_utf8_writer
    configure_utf8_stdio
    decode_argv_utf8
    open_utf8_writer
    print_utf8
    read_utf8_lines
);

decode_argv_utf8();
configure_utf8_stdio();

my %opts;
getopts('i:o:a:g:', \%opts);

my $usage = "Usage: $0 -i <label>_split_branch_label -o <label> -a <primitive_axis.tsv> -g <gene_id_map.tsv>\n";
my $in_dir = $opts{'i'} or die $usage;
my $label  = $opts{'o'} or die $usage;
my $axis_file = $opts{'a'} or die $usage;
my $gene_map_file = $opts{'g'} or die $usage;

die "[ERROR] Input directory not found: $in_dir\n" unless -d $in_dir;
die "[ERROR] Primitive axis ledger not found: $axis_file\n" unless -e $axis_file;
die "[ERROR] Gene identity map not found: $gene_map_file\n" unless -e $gene_map_file;

my $axis = read_axis_ledger($axis_file);
my @primitive = map { $_->{branch_id} } @{$axis};
my %primitive_id = map { $_ => 1 } @primitive;
my ($storage_to_gene, $gene_to_storage) = read_gene_id_map($gene_map_file);

# Collect input split files
opendir(my $DH, $in_dir) or die "[ERROR] Cannot open directory $in_dir: $!\n";
my @files = sort grep { /\.split\.txt$/ && -f File::Spec->catfile($in_dir, $_) } readdir($DH);
closedir($DH);

# gene_id -> { pattern => value }
my %gene_branch;
my %genes;

# Collect fused patterns (contain '|')
my %fused_patterns;

my %seen_storage;

foreach my $fn (@files) {
    my $path = File::Spec->catfile($in_dir, $fn);
    (my $storage_key = $fn) =~ s/\.split\.txt$//;
    die "[ERROR] Split file has no entry in gene identity map: $fn\n"
        unless exists $storage_to_gene->{$storage_key};
    die "[ERROR] Duplicate split file for storage key '$storage_key'.\n"
        if $seen_storage{$storage_key}++;
    my $gene_id = $storage_to_gene->{$storage_key};
    $genes{$gene_id} = 1;

    for my $line (@{ read_utf8_lines($path) }) {
        next if $line =~ /^\s*$/;

        my @f = split(/\t/, $line);
        next unless @f >= 3; # expected: split \t branch_pattern \t value

        my $pattern = $f[1];
        my $value   = $f[2];

        # Canonicalize fused order: B3|B1 -> B1|B3
        if ($pattern =~ /\|/) {
            my @parts = sort_branch_ids(split(/\|/, $pattern));
            validate_pattern_members(\@parts, \%primitive_id, $path, $pattern);
            $pattern = join('|', @parts);
            $fused_patterns{$pattern} = 1;
        } else {
            validate_pattern_members([$pattern], \%primitive_id, $path, $pattern);
        }

        if (exists $gene_branch{$gene_id}{$pattern}) {
            $gene_branch{$gene_id}{$pattern} = _merge_pattern_values(
                $gene_branch{$gene_id}{$pattern},
                $value,
            );
        } else {
            $gene_branch{$gene_id}{$pattern} = $value;
        }

    }
}

for my $storage_key (sort keys %{$storage_to_gene}) {
    die "[ERROR] Missing mapped split file for gene '$storage_to_gene->{$storage_key}' (storage key $storage_key).\n"
        unless $seen_storage{$storage_key};
}

# Sort fused patterns by their numeric components (stable, human-friendly)
my @fused_sorted = sort {
    _cmp_branch_pattern($a, $b)
} keys %fused_patterns;

# Output files
my $out_no   = "${label}.matrix_no_fuse.txt";
my $out_with = "${label}.matrix_with_fuse.txt";

my $NO = open_utf8_writer($out_no);
my $WITH = open_utf8_writer($out_with);

# Headers
print_utf8($NO, join("\t", "gene", @primitive), "\n");
print_utf8($WITH, join("\t", "gene", @primitive, @fused_sorted), "\n");

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
    print_utf8($NO, join("\t", @row_no), "\n");

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
    print_utf8($WITH, join("\t", @row_with), "\n");
}

close_utf8_writer($NO, $out_no);
close_utf8_writer($WITH, $out_with);

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
    return is_numeric_value($value);
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

sub validate_pattern_members {
    my ($parts, $axis_ids, $path, $pattern) = @_;
    for my $branch_id (@{$parts}) {
        die "[ERROR][$path] Invalid branch pattern '$pattern'.\n"
            unless defined $branch_id && $branch_id =~ /^B\d+$/;
        die "[ERROR][$path] Branch '$branch_id' in pattern '$pattern' is absent from the retained primitive axis ledger.\n"
            unless $axis_ids->{$branch_id};
    }
}

sub read_gene_id_map {
    my ($path) = @_;
    my $lines = read_utf8_lines($path);
    die "[ERROR] Gene identity map is empty: $path\n" unless @{$lines};
    my $header = $lines->[0];
    die "[ERROR] Unexpected gene identity map header in $path.\n"
        unless $header eq "storage_key\tgene_id\tinput_line";

    my (%storage_to_gene, %gene_to_storage);
    my $line_no = 1;
    for my $index (1 .. $#{$lines}) {
        my $line = $lines->[$index];
        $line_no++;
        next if $line eq '';
        my ($storage_key, $gene_id, $input_line) = split(/\t/, $line, -1);
        die "[ERROR][$path line $line_no] Malformed gene identity row.\n"
            unless defined $input_line && $storage_key ne '' && $gene_id ne '';
        die "[ERROR][$path line $line_no] Duplicate storage key '$storage_key'.\n"
            if exists $storage_to_gene{$storage_key};
        die "[ERROR][$path line $line_no] Duplicate gene ID '$gene_id'.\n"
            if exists $gene_to_storage{$gene_id};
        $storage_to_gene{$storage_key} = $gene_id;
        $gene_to_storage{$gene_id} = $storage_key;
    }
    die "[ERROR] Gene identity map has no records: $path\n" unless %storage_to_gene;
    return (\%storage_to_gene, \%gene_to_storage);
}
