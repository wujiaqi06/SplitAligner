#!/usr/bin/env perl
# ==============================================================================
# Script:      tree_to_splits.pl
# Author:      Jiaqi Wu (wujiaqi@hiroshima-u.ac.jp, wujiaqi06@gmail.com)
# Affiliation: Graduate School of Integrated Sciences for Life, Hiroshima University, Japan
#
# Description:
#   Parse Newick trees and decompose each edge into a bipartite split.
#   - Species mode (-m species): writes global coordinate axis files to CURRENT directory:
#         species_tree.splits.txt
#         species_tree.branch_map.txt
#     If -o is provided in species mode, it is used as an axis label and outputs
#     are written to <axis_label>.splits.txt and <axis_label>.branch_map.txt.
#   - Gene mode (-m gene -o <label>): writes per-gene split files to directory:
#         <label>_splits/*.split.txt
#
# Usage:
#   perl tree_to_splits.pl -i <tree_file> -m <species|gene> [-o <label>]
# ==============================================================================

use strict;
use warnings;
use Getopt::Std;

my %opts;
getopts('i:m:o:', \%opts);

my $input_tree_file = $opts{'i'}
  or die "Usage: $0 -i <tree_file> -m <species|gene> [-o <label>]\n";

my $mode = $opts{'m'} || 'gene';
die "[ERROR] -m must be 'species' or 'gene'\n" unless $mode eq 'species' || $mode eq 'gene';

my $opt_label = $opts{'o'};  # label is optional; meaning depends on mode

open(my $IN, '<', $input_tree_file) or die "[ERROR] Cannot open $input_tree_file: $!\n";

if ($mode eq 'species') {

    # In species mode, -o (if provided) is treated as an *axis label* to avoid
    # filename collisions in parallel runs. Default keeps backward compatibility.
    my $axis_label = defined $opt_label && $opt_label ne '' ? $opt_label : 'species_tree';

    # Species mode produces the global coordinate axis in the CURRENT directory.
    # To avoid accidental overwrite due to multi-line inputs, we only process the first non-empty line.
    my $line;
    my $lineno = 0;
    while (my $l = <$IN>) {
        $lineno++;
        $l =~ s/^\s+//;
        $l =~ s/\s+$//;
        next if $l eq '';
        $line = $l;
        last;
    }
    die "[ERROR] No valid tree line found in $input_tree_file\n" unless defined $line;

    my ($tree_id, $tree) = parse_gene_tree_line($line, $lineno);

    my @all_tips  = collect_tip_taxa($tree);
    my %subtrees  = parse_tree_components($tree);

    my $splits_file = "$axis_label.splits.txt";
    my $map_file    = "$axis_label.branch_map.txt";

    open(my $SPEC,  '>', $splits_file) or die "[ERROR] Cannot write $splits_file: $!\n";
    open(my $MAP,   '>', $map_file)    or die "[ERROR] Cannot write $map_file: $!\n";
    print {$MAP} "branch_id\tsub_tree\ttype\tnote\n";

    my @axis_entries;
    for my $node (subtree_nodes_in_branch_order(\%subtrees)) {
        my $branch_id = extract_branch_id($subtrees{$node});
        next unless defined $branch_id;

        my (@left, @right);
        my ($type, $canonical_split);
        if ($node =~ /\(/) {
            @left  = collect_tip_taxa($node);
            @right = subtract_array_multiset(\@all_tips, \@left);
            my $left_s  = join("..", sort @left);
            my $right_s = join("..", sort @right);
            $canonical_split = join("||", sort ($left_s, $right_s));
            $type = ($left_s ne '' && $right_s ne '') ? 'internal' : 'root';
        } else {
            @left  = ($node);
            @right = subtract_array_multiset(\@all_tips, \@left);
            my $left_s  = join("..", sort @left);
            my $right_s = join("..", sort @right);
            $canonical_split = join("||", sort ($left_s, $right_s));
            $type = 'terminal';
        }

        push @axis_entries, {
            branch_id => $branch_id,
            node      => $node,
            type      => $type,
            split     => $canonical_split,
        };
    }

    my (%split_rep, %split_losers);
    for my $entry (@axis_entries) {
        next if $entry->{type} eq 'root';
        my $split = $entry->{split};
        if (!exists $split_rep{$split}) {
            $split_rep{$split} = $entry;
            next;
        }

        my $current = $split_rep{$split};
        if (branch_id_num($entry->{branch_id}) < branch_id_num($current->{branch_id})) {
            push @{ $split_losers{$split} }, $current->{branch_id};
            $split_rep{$split} = $entry;
        } else {
            push @{ $split_losers{$split} }, $entry->{branch_id};
        }
    }

    for my $split (sort split_sort keys %split_rep) {
        my $entry = $split_rep{$split};
        print {$SPEC} "$split\t$entry->{branch_id}\n";
    }

    for my $entry (@axis_entries) {
        my $note = '';
        if ($entry->{type} ne 'root' && exists $split_rep{$entry->{split}}) {
            my $winner = $split_rep{$entry->{split}}{branch_id};
            if ($winner ne $entry->{branch_id}) {
                $note = "duplicate_unrooted_split_loser_of=$winner";
            } elsif (exists $split_losers{$entry->{split}}) {
                my @losers = sort_branch_ids(@{ $split_losers{$entry->{split}} });
                $note = "duplicate_unrooted_split_winner_over=" . join('|', @losers);
            }
        }
        print {$MAP} join("\t", $entry->{branch_id}, $entry->{node}, $entry->{type}, $note), "\n";
    }

    close $SPEC;
    close $MAP;

    # If there are additional non-empty lines, warn (do not fail).
    while (my $rest = <$IN>) {
        $rest =~ s/^\s+//; $rest =~ s/\s+$//;
        next if $rest eq '';
        warn "[WARN] Additional tree lines detected in species mode; only the first tree is processed.\n";
        last;
    }

} else {
    # Gene mode: write each gene into <label>_splits/<gene>.split.txt
    my $label   = defined $opt_label && $opt_label ne '' ? $opt_label : 'gene';
    my $out_dir = "${label}_splits";

    if (!-d $out_dir) {
        mkdir $out_dir or die "[ERROR] Cannot create directory $out_dir: $!\n";
    }

    my $lineno = 0;
    while (my $line = <$IN>) {
        $lineno++;
        chomp $line;
        next if $line =~ /^\s*$/;

        my ($gene_id, $tree) = parse_gene_tree_line($line, $lineno);

        # Safe filename (keep behavior close; only remove path separators)
        (my $gene_fn = $gene_id) =~ s{[\/\\]}{_}g;

        open(my $OUT, '>', "$out_dir/$gene_fn.split.txt")
          or die "[ERROR] Cannot write $out_dir/$gene_fn.split.txt: $!\n";

        my @all_tips = collect_tip_taxa($tree);
        my %subtrees = parse_tree_components($tree);

        foreach my $node (keys %subtrees) {
            $subtrees{$node} =~ s/\S+://;

            my (@left, @right);
            if ($node =~ /\(/) {
                @left  = collect_tip_taxa($node);
                @right = subtract_array_multiset(\@all_tips, \@left);
                my $left_s  = join("..", sort @left);
                my $right_s = join("..", sort @right);
                my $split   = join("||", sort ($left_s, $right_s));

                if ($left_s ne '' && $right_s ne '') {
                    print {$OUT} "$split\t$subtrees{$node}\n";
                }
            } else {
                @left  = ($node);
                @right = subtract_array_multiset(\@all_tips, \@left);
                my $left_s  = join("..", sort @left);
                my $right_s = join("..", sort @right);
                my $split   = join("||", sort ($left_s, $right_s));
                print {$OUT} "$split\t$subtrees{$node}\n";
            }
        }

        close $OUT;
    }
}

close $IN;
exit 0;

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------

sub parse_gene_tree_line {
    my ($line, $lineno) = @_;
    $lineno //= '?';
    $line =~ s/^\s+//;
    $line =~ s/\s+$//;

    die "[ERROR] Line $lineno is empty.\n" if $line eq '';
    die "[ERROR] Line $lineno has no '(' to start a tree.\n" unless $line =~ /\(/;

    my $pos  = index($line, '(');
    my $gene = substr($line, 0, $pos);
    my $tree = substr($line, $pos);

    $gene =~ s/\s+$//;
    $tree =~ s/^\s+//;

    die "[ERROR] Line $lineno has empty gene ID.\n" unless defined $gene && $gene ne '';
    die "[ERROR] Line $lineno has empty tree string.\n" unless defined $tree && $tree ne '';
    die "[ERROR] Duplicated semicolon in line $lineno.\n" if $tree =~ /;;/;
    die "[ERROR] Tree in line $lineno must end with ';'.\n" unless $tree =~ /;\s*$/;

    return ($gene, $tree);
}

sub collect_tip_taxa {
    my ($body) = @_;
    $body =~ s/\s+//g;
    $body =~ s/;//g;
    $body = strip_internal_node_labels($body);
    $body =~ s/:B\d+//g;
    $body =~ s/:\d*\.?\d+(?:[eE][+-]?\d+)?//g;
    $body =~ s/\)\d*\.?\d+(?:[eE][+-]?\d+)?/)/g;
    $body =~ s/[,\(\)]+/\t/g;

    my @tips = grep { $_ ne "" } split(/\t/, $body);
    my @uniq = uniq(@tips);
    if (@tips != @uniq) {
        warn "[WARN] Repeating terminal taxa detected in subtree (possible duplicated species entries).\n";
    }
    return @tips;
}

sub parse_tree_components {
    my ($body) = @_;
    my %tree_hash;
    $body =~ s/\s+//g;
    $body =~ s/;$//;
    $body = strip_internal_node_labels($body);

    my $terminal = $body;
    $terminal =~ s/[(),]+/\t/g;
    my @all_cuts = grep { $_ ne "" } split(/\t/, $terminal);

    foreach my $i (@all_cuts) {
        if (($i !~ /^\d+:/) && ($i !~ /^:/)) {
            my ($term_taxa, $anno) = split(/:/, $i);
            $tree_hash{$term_taxa} = $anno;
        }
    }

    my $nwk = $body;
    my $len = length($nwk);

    for my $i (0 .. $len - 1) {
        next unless substr($nwk, $i, 1) eq '(';

        my $depth = 0;
        my $match_end = -1;

        for my $j ($i .. $len - 1) {
            my $char = substr($nwk, $j, 1);
            $depth++ if $char eq '(';
            $depth-- if $char eq ')';
            if ($depth == 0) {
                $match_end = $j;
                last;
            }
        }

        if ($match_end != -1) {
            my $subtree = substr($nwk, $i, $match_end - $i + 1);
            my $remainder = substr($nwk, $match_end + 1);
            $remainder =~ s/[(),]+.*$//;
            $remainder =~ s/^://;
            $tree_hash{$subtree} = $remainder;
        }
    }
    return %tree_hash;
}

sub strip_internal_node_labels {
    my ($body) = @_;
    # Remove string node labels that appear immediately after a closing
    # parenthesis, e.g. )N_11:0.5 -> ):0.5. This preserves branch lengths and
    # SplitAligner species-tree branch tags such as ):B11.
    $body =~ s/\)([^:;,\(\)]+)(?=[:),;]|$)/)/g;
    return $body;
}

sub extract_branch_id {
    my ($annotation) = @_;
    return undef unless defined $annotation;
    if ($annotation =~ /(?:^|:)(B\d+)(?:$|[:;])/ ) {
        return $1;
    }
    return undef;
}

sub subtract_array_multiset {
    my ($arr1_ref, $arr2_ref) = @_;
    my %count2;
    my @result;

    foreach my $taxa (@{$arr2_ref}) { $count2{$taxa}++; }
    foreach my $taxa (@{$arr1_ref}) {
        if (exists $count2{$taxa} && $count2{$taxa} > 0) {
            $count2{$taxa}--;
        } else {
            push @result, $taxa;
        }
    }
    return @result;
}

sub uniq {
    my %seen;
    return grep { !$seen{$_}++ } @_;
}

sub branch_id_num {
    my ($branch_id) = @_;
    return 9**9**9 unless defined $branch_id && $branch_id =~ /^B(\d+)$/;
    return $1;
}

sub sort_branch_ids {
    return sort {
        branch_id_num($a) <=> branch_id_num($b) || $a cmp $b
    } @_;
}

sub subtree_nodes_in_branch_order {
    my ($subtrees_ref) = @_;
    return sort {
        branch_id_num(extract_branch_id($subtrees_ref->{$a}))
            <=>
        branch_id_num(extract_branch_id($subtrees_ref->{$b}))
        ||
        $a cmp $b
    } keys %{$subtrees_ref};
}

sub split_sort {
    return $a cmp $b;
}
