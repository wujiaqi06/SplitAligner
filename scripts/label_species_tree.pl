#!/usr/bin/env perl -w
# ==============================================================================
# Script:      label_species_tree.pl
# Author:      Jiaqi Wu (wujiaqi@hiroshima-u.ac.jp, wujiaqi06@gmail.com)
# Affiliation: Graduate School of Integrated Sciences for Life, Hiroshima University, Japan
# Copyright:   (c) 2026 Jiaqi Wu. All rights reserved.
# 
# Description: Prepares the backbone species tree by assigning stable internal 
#              branch identifiers (e.g., B1, B2) to create the reference coordinate 
#              system for all downstream split-based mappings.
# 
# Usage:       perl label_species_tree.pl -i <species_tree.nwk> -o <output_prefix> [-l <record_id>]
# ==============================================================================

use 5.010;
use strict;
use Cwd;
use Getopt::Std;
use FindBin;
my $script_dir = $FindBin::RealBin;

# [Original clean code continues below...]
# Node labels are temporary and not topology-invariant.
# Stable branch identity will be defined in the split generation step.

my %opts;
getopts('i:o:l:', \%opts);
my $species_tree_file = $opts{'i'} or die usage();
my $output_prefix = $opts{'o'} or die usage();
my $record_id = $opts{'l'} // 'species_tree';  # record id used in *.forSplit.nwk

sub usage {
    return "Usage: perl $0 -i <species_tree.nwk> -o <output_prefix> [-l <record_id>]\n"
         . "  -i  Input species tree in Newick format (single tree)\n"
         . "  -o  Output prefix; writes <prefix>.forSplit.nwk and <prefix>.FigTree.tre\n"
         . "  -l  Optional record id used in <prefix>.forSplit.nwk (default: species_tree)\n";
}


my $raw = &read_single_newick_tree_or_die($species_tree_file);
my $species_nwk = &sanitize_newick_regex_full($raw);
&check_newick_basic_or_die($species_nwk, $species_tree_file);

my $tip_labels_ref = tip_node_labels($species_nwk);

die "[ERROR] Input tree already contains branch labels like :B123\n"
    if $species_nwk =~ /:B\d+/;

my $n_tip = scalar(keys %{$tip_labels_ref});
die "[ERROR] Species tree must contain at least 3 taxa.\n" if $n_tip < 3;

my $for_split = $species_nwk;
my $for_figtree = $species_nwk;

#add tip labels from 1 to n_tip in $for_split
my @tips = sort { length($b) <=> length($a) } keys %{$tip_labels_ref};

my %tmp2final;
my $k = 1;

for my $tax (@tips) {
    my $tmp = "__TIP__$k" . "__";
    $tmp2final{$tmp} = $tax . ":B" . $tip_labels_ref->{$tax};
    $for_split =~ s/\Q$tax\E(?=[,\)\;])/$tmp/g;
    $k++;
}

for my $tmp (sort keys %tmp2final) {
    $for_split =~ s/\Q$tmp\E/$tmp2final{$tmp}/g;
}

#add tip labels from 1 to n_tip in $for_figtree
my @tips2 = sort { length($b) <=> length($a) } keys %{$tip_labels_ref};

my %tmp2final_fig;
my $m = 1;

for my $tax (@tips2) {
    my $tmp = "__TIPFIG__$m" . "__";
    $tmp2final_fig{$tmp} = $tax . "_B" . $tip_labels_ref->{$tax};
    $for_figtree =~ s/\Q$tax\E(?=[,\)\;])/$tmp/g;
    $m++;
}

for my $tmp (sort keys %tmp2final_fig) {
    $for_figtree =~ s/\Q$tmp\E/$tmp2final_fig{$tmp}/g;
}


my $total_close = ($species_nwk =~ tr/)//);

my $seen = 0;
my $count = $n_tip + 1;
$for_split =~ s/\)/(++$seen < $total_close) ? "):B" . $count++ : ")"/eg;

$seen = 0;
$count = $n_tip + 1;
$for_figtree =~ s/\)/(++$seen < $total_close) ? ')"B' . $count++ . '"' : ')'/eg;

#print "$for_split\n";
#print "$for_figtree\n";
die "[ERROR] Placeholder remained in forSplit output.\n"
    if $for_split =~ /__TIP__\d+__/;
die "[ERROR] Placeholder remained in FigTree output.\n"
    if $for_figtree =~ /__TIPFIG__\d+__/;

open OUT1, ">$output_prefix.forSplit.nwk" or die $!;
open OUT2, ">$output_prefix.FigTree.tre" or die $!;
print OUT1 "$record_id$for_split\n";
print OUT2 "$for_figtree\n";

close OUT1;
close OUT2;



sub check_newick_file_or_die {
    my ($path) = @_;
    open my $fh, "<", $path or die "[ERROR] Cannot open $path: $!\n";
    local $/ = undef;
    my $s = <$fh>;
    close $fh;

    check_newick_basic_or_die($s, $path);
    return $s; 
}

sub check_newick_basic_or_die {
    my ($s, $name) = @_;
    $name //= 'Newick';

    die "[ERROR][$name] Empty string.\n" if !defined($s) || $s eq '';

    # no whitespace is expected (you already strip), but tolerate anyway
    $s =~ s/\s+//g;

    # must end with ';'
    die "[ERROR][$name] Newick must end with ';'.\n" if $s !~ /;$/;

    my $depth    = 0;
    my $in_quote = 0;

    my $len = length($s);
    for (my $i = 0; $i < $len; $i++) {
        my $ch = substr($s, $i, 1);

        if ($ch eq "'") { $in_quote = !$in_quote; next; }
        next if $in_quote;

        if ($ch eq '(') {
            $depth++;
        } elsif ($ch eq ')') {
            $depth--;
            die "[ERROR][$name] Unmatched ')' at index $i\n" if $depth < 0;
        }
    }

    die "[ERROR][$name] Unclosed single quote (').\n" if $in_quote;
    die "[ERROR][$name] Unmatched '(' (depth remaining=$depth).\n" if $depth != 0;

    die "[ERROR][$name] Empty subtree '()' detected.\n" if $s =~ /\(\)/;
    die "[ERROR][$name] Empty child '(,' detected.\n"  if $s =~ /\(,/;
    die "[ERROR][$name] Empty child ',)' detected.\n"  if $s =~ /,\)/;
    die "[ERROR][$name] Double comma ',,' detected.\n" if $s =~ /,,/;
    # crude extra guard: semicolon should appear exactly once at top level (after sanitization it should)
    my $top_semicolons = 0;
    $depth    = 0;
    $in_quote = 0;
    for (my $i = 0; $i < $len; $i++) {
        my $ch = substr($s, $i, 1);

        if ($ch eq "'") { $in_quote = !$in_quote; next; }
        next if $in_quote;

        $depth++ if $ch eq '(';
        $depth-- if $ch eq ')';

        if ($ch eq ';' && $depth == 0) { $top_semicolons++; }
    }
    die "[ERROR][$name] Found $top_semicolons top-level ';'. Expect exactly 1 tree.\n"
        if $top_semicolons != 1;

    return 1;
}

sub read_single_newick_tree_or_die {
    my ($path) = @_;
    open my $fh, "<", $path or die "[ERROR] Cannot open $path: $!\n";
    local $/ = undef;
    my $s = <$fh>;
    close $fh;

    die "[ERROR][$path] Empty file.\n" if !defined($s) || $s eq '';

    # remove whitespace
    $s =~ s/\s+//g;

    # Strip leading/trailing empties
    $s =~ s/^\s+//;
    $s =~ s/\s+$//;

    # Must contain at least one ';'
    die "[ERROR][$path] No ';' found. Not a valid Newick.\n" if $s !~ /;/;

    # Scan to find top-level tree boundaries: a ';' that occurs when depth==0 and not in quotes.
    my $depth = 0;
    my $in_quote = 0;
    my @top_level_semicolons;

    my $len = length($s);
    for (my $i = 0; $i < $len; $i++) {
        my $ch = substr($s, $i, 1);

        if ($ch eq "'") { $in_quote = !$in_quote; next; }
        next if $in_quote;

        if ($ch eq '(') { $depth++; }
        elsif ($ch eq ')') {
            $depth--;
            die "[ERROR][$path] Unmatched ')' at index $i\n" if $depth < 0;
        }
        elsif ($ch eq ';') {
            # semicolon ending a top-level tree must appear at depth 0
            push @top_level_semicolons, $i if $depth == 0;
        }
    }

    die "[ERROR][$path] Unclosed single quote.\n" if $in_quote;
    die "[ERROR][$path] Unmatched '(' (depth remaining=$depth).\n" if $depth != 0;

    # Must be exactly one top-level ';'
    if (@top_level_semicolons != 1) {
        my $n = scalar(@top_level_semicolons);
        die "[ERROR][$path] Species tree file must contain EXACTLY 1 Newick tree. Found $n top-level ';'.\n";
    }

    # Extract that one tree (everything up to the only top-level ';')
    my $end = $top_level_semicolons[0];
    my $tree = substr($s, 0, $end + 1);

    # Ensure there's no extra non-empty content after it
    my $tail = substr($s, $end + 1);
    die "[ERROR][$path] Extra content after the single Newick tree.\n" if defined($tail) && $tail ne '';

    return $tree;
}

sub sanitize_newick_regex_full {
    my ($s) = @_;
    $s = sanitize_newick_regex_basic($s);

    # 3) remove branch lengths: :<anything until delimiter>
    #    handles :0.01, :1e-3, :0.01B123, etc.
    $s =~ s/:[^,\)\;]+(?=[,\)\;])//g;

    return $s;
}

sub sanitize_newick_regex_basic {
    my ($s) = @_;
    $s =~ s/\s+//g;

    # 1) remove bracket comments: [...]
    $s =~ s/\[[^\[\]]*\]//g;

    # 2) remove internal node labels: )LABEL  where LABEL has no delimiters
    #    keep ')' itself, stop before : , ) ;
    $s =~ s/\)([^:,\)\;]+)(?=[:,\)\;])/\)/g;

    return $s;
}

sub tip_node_labels {
    my ($s) = @_;
    #print "$s\n";
    $s =~ s/;//;
    $s =~ s/[,\(\)]+/\t/g;
    #print "$s\n";
    my @tips = grep { $_ ne "" } split(/\t/, $s);
    my %node_numbers;
    my $count = 0;
    foreach my $i (@tips){
    die "[ERROR] Illegal character in taxon name: $i\n"
        if $i =~ /[(),:;]/;
        $count ++;
        $node_numbers{$i} = $count;
    }
    my @tip_uniq = &uniq(@tips);
    if (scalar(@tips) != scalar(@tip_uniq)){
        die "ERROR! Repeating terminal taxa detected!";
    }
    return \%node_numbers;
}

sub uniq{
    my %seen;
    return grep {!$seen{$_}++} @_;
}