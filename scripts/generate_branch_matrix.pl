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
use POSIX qw(isfinite);
use File::Spec;
use FindBin qw($RealBin);
use lib "$RealBin/../lib";
use SplitAligner::Newick qw(format_pooled_length is_numeric_value sum_component_lengths);
use SplitAligner::Provenance qw(read_axis_ledger);
use SplitAligner::TextIO qw(
    close_utf8_writer
    configure_utf8_stdio
    decode_argv_utf8
    open_utf8_writer
    print_utf8
    for_each_utf8_line
    read_utf8_lines
);

decode_argv_utf8();
configure_utf8_stdio();

my %opts;
getopts('i:o:a:g:v:', \%opts) or die "Invalid matrix options.\n";

my $usage = "Usage: $0 -i <label>_split_branch_label -o <label> -a <primitive_axis.tsv> -g <gene_id_map.tsv>\n";
my $in_dir = $opts{'i'} or die $usage;
my $label  = $opts{'o'} or die $usage;
my $axis_file = $opts{'a'} or die $usage;
my $gene_map_file = $opts{'g'} or die $usage;

my $view = $opts{v} // 'whole';
die "[ERROR] Fusion view must be whole or direct.\n" unless $view =~ /\A(?:whole|direct)\z/;

die "[ERROR] Input directory not found: $in_dir\n" unless -d $in_dir;
die "[ERROR] Primitive axis ledger not found: $axis_file\n" unless -e $axis_file;
die "[ERROR] Gene identity map not found: $gene_map_file\n" unless -e $gene_map_file;

my $axis = read_axis_ledger($axis_file);
my @primitive = map { $_->{branch_id} } @{$axis};
my %primitive_id = map { $_ => 1 } @primitive;
my (undef, $gene_to_storage) = read_gene_id_map($gene_map_file);

# The shared axis comes from reference restriction, even for a group whose
# projected split is absent from every gene tree where it forms.
my %targets;
my $groups_path = File::Spec->catfile($in_dir, "$label.projected_groups.tsv");
for_each_utf8_line($groups_path, sub {
    my ($line, $number) = @_;
    return if $number == 1;
    my (undef, undef, $pattern) = split(/\t/, $line);
    return if exists $targets{$pattern};
    my @parts = split(/\|/, $pattern);
    validate_pattern_members(\@parts, \%primitive_id, $groups_path, $pattern);
    $targets{$pattern} = \@parts;
});
my @fused_sorted = sort { _cmp_branch_pattern($a, $b) } keys %targets;

my $out_no = "${label}.matrix_no_fuse.txt";
my $out_with = "${label}.matrix_with_fuse.txt";
my $out_cells = "${label}.fusion_cells.tsv";
my $NO = open_utf8_writer($out_no);
my $WITH = open_utf8_writer($out_with);
my $CELLS = open_utf8_writer($out_cells);
print_utf8($NO, join("\t", "gene", @primitive), "\n");
print_utf8($WITH, join("\t", "gene", @primitive, @fused_sorted), "\n");
print_utf8($CELLS, join("\t", qw(gene column reference_state recovery_status numeric_status
    construction_mode components deleted_members cut_components length)), "\n");

for my $gene_id (sort keys %{$gene_to_storage}) {
    my $storage = $gene_to_storage->{$gene_id};
    my $mapped = File::Spec->catfile($in_dir, "$storage.split.txt");
    my %primitive_values;
    for_each_utf8_line($mapped, sub {
        my ($line) = @_;
        my (undef, $pattern, $length) = split(/\t/, $line);
        $primitive_values{$pattern} = $length if $primitive_id{$pattern};
    });
    my @row = ($gene_id, map { $primitive_values{$_} // 'NA' } @primitive);
    print_utf8($NO, join("\t", @row), "\n");
    if (@fused_sorted) {
        my $component_path = File::Spec->catfile($in_dir, "$storage.components.tsv");
        my $owners = read_components($component_path);
        for my $target (@fused_sorted) {
            my @cell = construct_cell($targets{$target}, $owners, $view);
            push @row, $cell[-1];
            print_utf8($CELLS, join("\t", $gene_id, $target, @cell), "\n");
        }
    }
    print_utf8($WITH, join("\t", @row), "\n");
}
close_utf8_writer($NO, $out_no);
close_utf8_writer($WITH, $out_with);
close_utf8_writer($CELLS, $out_cells);

sub read_components {
    my ($path) = @_;
    my %owner;
    for_each_utf8_line($path, sub {
        my ($line, $number) = @_;
        return if $number == 1;
        my ($pattern, $source_row, $unmatched_split, $length) = split(/\t/, $line, -1);
        my @parts = split(/\|/, $pattern);
        my $component = { id => $pattern, size => scalar(@parts),
            recovered => length($source_row), length => $length };
        for my $part (@parts) {
            die "[ERROR][$path line $number] Repeated component membership for $part.\n"
                if exists $owner{$part};
            $owner{$part} = $component;
        }
    });
    return \%owner;
}

sub construct_cell {
    my ($target, $owner, $view) = @_;
    my (@parts, @deleted, %intersection_size);
    for my $branch (@{$target}) {
        my $c = $owner->{$branch};
        if (!defined $c) {
            push @deleted, $branch;
            next;
        }
        push @parts, $c unless $intersection_size{$c->{id}}++;
    }
    my @cut = map { $_->{id} } grep { $intersection_size{$_->{id}} != $_->{size} } @parts;
    my $ids = join(',', map { $_->{id} } @parts);
    if (@deleted || @cut) {
        return ('unavailable', 'not_applicable', 'not_applicable', 'none',
            $ids, join(',', @deleted), join(',', @cut), 'NA');
    }

    my $active = @parts == 1;
    my $reference = $active ? 'active' : 'decomposable';
    my $mode = $active ? 'single_edge_lookup' : 'selected_component_edge_sum';
    return ($reference, 'not_recovered', 'not_recovered', $mode, $ids, '', '', 'NA')
        if grep { !$_->{recovered} } @parts;
    # Direct view suppresses cross-component construction before inspecting
    # lengths. The current-component evidence is preserved independently.
    return ($reference, 'recovered', 'direct_only', $mode, $ids, '', '', 'NA')
        if !$active && $view eq 'direct';
    return ($reference, 'recovered', 'missing_length', $mode, $ids, '', '', 'NA')
        if grep { !is_numeric_value($_->{length}) || !isfinite(0 + $_->{length}) } @parts;
    my $value = $parts[0]{length};
    if (!$active) {
        $value = sum_component_lengths(map { $_->{length} } @parts);
        return ($reference, 'recovered', 'nonfinite_sum', $mode, $ids, '', '', 'NA')
            unless isfinite($value);
        $value = format_pooled_length($value);
    }
    return ($reference, 'recovered', 'available', $mode, $ids, '', '', $value);
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
