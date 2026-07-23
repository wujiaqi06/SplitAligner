#!/usr/bin/env perl
# ==============================================================================
# Script:      label_species_tree.pl
# Author:      Jiaqi Wu (wujiaqi@hiroshima-u.ac.jp, wujiaqi06@gmail.com)
# Description: Parse one species tree structurally and assign stable B aliases
#              using the historical SplitAligner traversal order.
# ==============================================================================

use strict;
use warnings;
use Getopt::Std;
use FindBin qw($RealBin);
use lib "$RealBin/../lib";
use SplitAligner::Newick qw(
    assign_branch_ids
    parse_species_tree_file
    serialize_figtree
    serialize_for_split
    tip_labels
);
use SplitAligner::TextIO qw(
    configure_utf8_stdio
    decode_argv_utf8
    write_utf8_file
);

decode_argv_utf8();
configure_utf8_stdio();

my %opts;
getopts('i:o:l:', \%opts);
my $species_tree_file = $opts{i} or die usage();
my $output_prefix = $opts{o} or die usage();
my $record_id = defined $opts{l} ? $opts{l} : 'species_tree';

my $root = parse_species_tree_file($species_tree_file);
my @tips = tip_labels($root);
die "[ERROR] Species tree must contain at least 3 taxa.\n" if @tips < 3;
assign_branch_ids($root);

my $for_split = "$output_prefix.forSplit.nwk";
my $for_figtree = "$output_prefix.FigTree.tre";

write_utf8_file($for_split, serialize_for_split($root, $record_id) . "\n");
write_utf8_file($for_figtree, serialize_figtree($root) . "\n");

sub usage {
    return "Usage: perl $0 -i <species_tree.nwk> -o <output_prefix> [-l <record_id>]\n"
         . "  -i  Input species tree in Newick format (exactly one tree)\n"
         . "  -o  Output prefix; writes <prefix>.forSplit.nwk and <prefix>.FigTree.tre\n"
         . "  -l  Optional record id used in <prefix>.forSplit.nwk (default: species_tree)\n";
}
