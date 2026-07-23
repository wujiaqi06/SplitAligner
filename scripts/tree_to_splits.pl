#!/usr/bin/env perl
# ==============================================================================
# Script:      tree_to_splits.pl
# Author:      Jiaqi Wu (wujiaqi@hiroshima-u.ac.jp, wujiaqi06@gmail.com)
# Description: Convert structurally parsed species or gene trees to canonical
#              unrooted edge splits.
# ==============================================================================

use strict;
use warnings;
use Getopt::Std;
use FindBin qw($RealBin);
use File::Path qw(make_path);
use File::Spec;
use lib "$RealBin/../lib";
use SplitAligner::Newick qw(
    branch_id_num
    edge_records
    gene_split_records
    parse_gene_tree_file
    parse_labeled_species_tree_file
    serialize_branch_subtree
);
use SplitAligner::Provenance qw(write_axis_ledger);
use SplitAligner::TextIO qw(
    close_utf8_writer
    configure_utf8_stdio
    decode_argv_utf8
    open_utf8_writer
    print_utf8
);

decode_argv_utf8();
configure_utf8_stdio();

my %opt;
getopts('i:m:o:', \%opt);
my $input = $opt{i} // '';
my $mode  = $opt{m} // '';
my $label = defined $opt{o} && $opt{o} ne '' ? $opt{o} : ($mode eq 'species' ? 'species_tree' : 'gene');

die usage() unless $input ne '' && ($mode eq 'species' || $mode eq 'gene');
die "[ERROR] Input tree file not found: $input\n" unless -e $input;

if ($mode eq 'species') {
    write_species_axis($input, $label);
} else {
    write_gene_splits($input, $label);
}

exit 0;

sub write_species_axis {
    my ($path, $axis_label) = @_;
    my (undef, $root) = parse_labeled_species_tree_file($path);
    my $edges = edge_records($root);

    my %seen_branch;
    for my $edge (@{$edges}) {
        my $branch_id = $edge->{branch_id};
        die "[ERROR][$path] A non-root species-tree edge lacks a B alias.\n"
            unless defined $branch_id && $branch_id =~ /^B\d+$/;
        die "[ERROR][$path] Duplicate B alias '$branch_id'.\n" if $seen_branch{$branch_id}++;
    }

    my %by_split;
    push @{ $by_split{$_->{split}} }, $_ for @{$edges};

    my (%winner_for_split, %losers_for_split);
    for my $split (keys %by_split) {
        my @ordered = sort {
            branch_id_num($a->{branch_id}) <=> branch_id_num($b->{branch_id})
        } @{ $by_split{$split} };
        $winner_for_split{$split} = $ordered[0];
        $losers_for_split{$split} = [ @ordered[1 .. $#ordered] ] if @ordered > 1;
    }

    my $splits_path = "$axis_label.splits.txt";
    my $splits_fh = open_utf8_writer($splits_path);
    for my $split (sort keys %winner_for_split) {
        print_utf8($splits_fh, "$split\t$winner_for_split{$split}{branch_id}\n");
    }
    close_utf8_writer($splits_fh, $splits_path);

    my $map_path = "$axis_label.branch_map.txt";
    my $map_fh = open_utf8_writer($map_path);
    print_utf8($map_fh, "branch_id\tsub_tree\ttype\tnote\n");
    for my $edge (sort {
        branch_id_num($a->{branch_id}) <=> branch_id_num($b->{branch_id})
    } @{$edges}) {
        my $winner = $winner_for_split{$edge->{split}}{branch_id};
        my $note = '';
        if ($winner ne $edge->{branch_id}) {
            $note = "duplicate_unrooted_split_loser_of=$winner";
        } elsif ($losers_for_split{$edge->{split}}) {
            my @losers = map { $_->{branch_id} } @{ $losers_for_split{$edge->{split}} };
            $note = 'duplicate_unrooted_split_winner_over=' . join('|', @losers);
        }
        print_utf8($map_fh, join("\t",
            $edge->{branch_id},
            serialize_branch_subtree($edge->{node}),
            $edge->{type},
            $note,
        ), "\n");
    }
    close_utf8_writer($map_fh, $map_path);

    my @axis = map {
        {
            branch_id => $winner_for_split{$_}{branch_id},
            split     => $_,
        }
    } sort {
        branch_id_num($winner_for_split{$a}{branch_id})
            <=> branch_id_num($winner_for_split{$b}{branch_id})
    } keys %winner_for_split;
    write_axis_ledger("$axis_label.primitive_axis.tsv", \@axis);
}

sub write_gene_splits {
    my ($path, $output_label) = @_;

    # Whole-file preflight happens before any output path is created.
    my $records = parse_gene_tree_file($path);
    for my $record (@{$records}) {
        die "[ERROR][$path line $record->{line_no}] Gene ID contains a tab and cannot be represented safely: '$record->{gene_id}'.\n"
            if $record->{gene_id} =~ /\t/;
    }

    my @ordered = sort {
        $a->{gene_id} cmp $b->{gene_id} || $a->{line_no} <=> $b->{line_no}
    } @{$records};
    my $out_dir = "${output_label}_splits";
    make_path($out_dir);

    my $map_path = "$output_label.gene_id_map.tsv";
    my $map_fh = open_utf8_writer($map_path);
    print_utf8($map_fh, "storage_key\tgene_id\tinput_line\n");

    my $ordinal = 0;
    for my $record (@ordered) {
        my $storage_key = sprintf('g%08d', ++$ordinal);
        print_utf8($map_fh, join("\t", $storage_key, $record->{gene_id}, $record->{line_no}), "\n");

        my $out_path = File::Spec->catfile($out_dir, "$storage_key.split.txt");
        my $out_fh = open_utf8_writer($out_path);
        my $split_records = gene_split_records($record->{root});
        for my $split (@{$split_records}) {
            print_utf8($out_fh, "$split->{split}\t$split->{value}\n");
        }
        close_utf8_writer($out_fh, $out_path);
    }
    close_utf8_writer($map_fh, $map_path);
}

sub usage {
    return "Usage: perl $0 -i <tree_file> -m <species|gene> -o <label>\n";
}
