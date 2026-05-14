#!/usr/bin/env perl
# ==============================================================================
# Script:      SplitAligner.pl
# Author:      Jiaqi Wu (wujiaqi@hiroshima-u.ac.jp)
# Description: Main controller for the SplitAligner workflow.
#
# Modes:
#   1) matrix
#      Generate gene × branch matrices from a species tree and a set of gene trees.
#
#      Required arguments:
#        --species <species_tree.nwk>
#        --gene    <gene_trees.nwk>
#        --label   <output_label>
#
#      Output files:
#        species_tree.forSplit.nwk
#        species_tree.FigTree.tre
#        species_tree.splits.txt
#        species_tree.branch_map.txt
#        <label>_splits/
#        <label>_split_branch_label/
#        <label>.matrix_no_fuse.txt
#        <label>.matrix_with_fuse.txt
#
#   2) finalize
#      Convert <label>.matrix_with_fuse.txt matrices into final NA-classified
#      outputs by applying NA_fuse followed by NA_struct / NA_topo diagnosis.
#
#      Required arguments:
#        --free        <free.matrix_with_fuse.txt>
#        --fix         <fix.matrix_with_fuse.txt>
#        --final_label <output_prefix>
#      Optional arguments:
#        --species_tree <species_tree.forSplit.nwk>
#
#      Output files:
#        <free>.na_fuse.txt
#        <fix>.na_fuse.txt
#        <final_label>.fix.na_classified.txt
#        <final_label>.free.na_classified.txt
#        <final_label>.support_b.txt
#        <species_prefix>.support_b.nwk
#
#   3) finalize_fix
#      Finalize a fixed-topology matrix when no free-topology comparator is
#      available. NA_fuse is identified first; remaining NA is written as
#      NA_struct.
#
#      Required arguments:
#        --fix         <fix.matrix_with_fuse.txt>
#        --final_label <output_prefix>
#
#      Output files:
#        <fix>.na_fuse.txt
#        <final_label>.fix.na_classified.txt
#
# Usage examples:
#   SplitAligner.pl --mode matrix --species speciesTree302.nwk \
#     --gene free_tree.examples.nwk --label free
#
#   SplitAligner.pl --mode finalize --free free.matrix_with_fuse.txt \
#     --fix fix.matrix_with_fuse.txt --final_label final \
#     --species_tree species_tree.forSplit.nwk
# ==============================================================================

use strict;
use warnings;
use Getopt::Long qw(GetOptions);
use FindBin qw($RealBin);
use File::Spec;
use Cwd qw(abs_path);

my ($mode, $species_tree, $gene_trees, $label, $free_matrix, $fix_matrix, $final_label, $help);

GetOptions(
    'mode=s'        => \$mode,
    'species=s'     => \$species_tree,
    'gene=s'        => \$gene_trees,
    'label=s'       => \$label,
    'free=s'        => \$free_matrix,
    'fix=s'         => \$fix_matrix,
    'final_label=s' => \$final_label,
    'species_tree=s' => \$species_tree,
    'help|h'        => \$help,
) or die usage();

if ($help) {
    print usage();
    exit 0;
}

die usage() unless defined $mode;

my $scripts_dir = File::Spec->catdir($RealBin, 'scripts');
for my $required_script (
    qw(label_species_tree.pl tree_to_splits.pl split_branch_label.pl
       generate_branch_matrix.pl extract_na_fuse.pl confirm_na_structure.pl
       classify_fix_missingness.pl)
) {
    my $path = File::Spec->catfile($scripts_dir, $required_script);
    die "[ERROR] Required helper script not found: $path\n" unless -e $path;
}

if ($mode eq 'matrix') {
    run_matrix_mode(
        species => $species_tree,
        gene    => $gene_trees,
        label   => $label,
    );
}
elsif ($mode eq 'finalize') {
    run_finalize_mode(
        free        => $free_matrix,
        fix         => $fix_matrix,
        final_label => $final_label,
        species_tree => $species_tree,
    );
}
elsif ($mode eq 'finalize_fix') {
    run_finalize_fix_mode(
        fix         => $fix_matrix,
        final_label => $final_label,
    );
}
else {
    die "[ERROR] Unsupported --mode: $mode\n" . usage();
}

exit 0;

# ------------------------------------------------------------------------------
# matrix mode
# ------------------------------------------------------------------------------
sub run_matrix_mode {
    my %arg = @_;

    my $species = $arg{species};
    my $gene    = $arg{gene};
    my $label0  = $arg{label};

    die "[ERROR] --species is required in --mode matrix\n" unless defined $species && $species ne '';
    die "[ERROR] --gene is required in --mode matrix\n"    unless defined $gene    && $gene ne '';
    die "[ERROR] --label is required in --mode matrix\n"   unless defined $label0  && $label0 ne '';

    die "[ERROR] Species tree file not found: $species\n" unless -e $species;
    die "[ERROR] Gene tree file not found: $gene\n"       unless -e $gene;
    validate_label($label0, '--label');

    my $species_prefix = 'species_tree';

    print STDERR "[INFO] Running SplitAligner matrix mode for label: $label0\n";
    print STDERR "[INFO] Species tree: $species\n";
    print STDERR "[INFO] Gene trees   : $gene\n";

    run_perl_script(
        'label_species_tree.pl',
        ['-i', $species, '-o', $species_prefix, '-l', $species_prefix],
        'Label species tree',
    );

    run_perl_script(
        'tree_to_splits.pl',
        ['-i', "$species_prefix.forSplit.nwk", '-m', 'species', '-o', $species_prefix],
        'Generate species split axis',
    );

    run_perl_script(
        'tree_to_splits.pl',
        ['-i', $gene, '-m', 'gene', '-o', $label0],
        'Generate gene splits',
    );

    run_perl_script(
        'split_branch_label.pl',
        ['-i', "${label0}_splits", '-j', "$species_prefix.splits.txt", '-o', $label0],
        'Map gene splits onto species axis',
    );

    run_perl_script(
        'generate_branch_matrix.pl',
        ['-i', "${label0}_split_branch_label", '-o', $label0],
        'Generate branch matrices',
    );

    for my $must_exist (
        "$species_prefix.branch_map.txt",
        "$label0.matrix_no_fuse.txt",
        "$label0.matrix_with_fuse.txt",
    ) {
        die "[ERROR] Expected output not found after matrix mode: $must_exist\n" unless -e $must_exist;
    }

    print STDERR "[INFO] Matrix mode completed successfully for label: $label0\n";
    print STDERR "[INFO] Generated: $label0.matrix_no_fuse.txt\n";
    print STDERR "[INFO] Generated: $label0.matrix_with_fuse.txt\n";
}

# ------------------------------------------------------------------------------
# finalize mode
# ------------------------------------------------------------------------------
sub run_finalize_mode {
    my %arg = @_;

    my $free    = $arg{free};
    my $fix     = $arg{fix};
    my $final   = $arg{final_label};
    my $species = $arg{species_tree};

    die "[ERROR] --free is required in --mode finalize\n"        unless defined $free  && $free ne '';
    die "[ERROR] --fix is required in --mode finalize\n"         unless defined $fix   && $fix ne '';
    die "[ERROR] --final_label is required in --mode finalize\n" unless defined $final && $final ne '';

    die "[ERROR] FREE matrix file not found: $free\n" unless -e $free;
    die "[ERROR] FIX matrix file not found: $fix\n"   unless -e $fix;
    if (defined $species && $species ne '') {
        die "[ERROR] Species tree file not found: $species\n" unless -e $species;
    }
    validate_label($final, '--final_label');

    print STDERR "[INFO] Running SplitAligner finalize mode\n";
    print STDERR "[INFO] FREE matrix : $free\n";
    print STDERR "[INFO] FIX matrix  : $fix\n";

    run_perl_script(
        'extract_na_fuse.pl',
        ['-i', $free],
        'Mark NA_fuse in free matrix',
    );

    run_perl_script(
        'extract_na_fuse.pl',
        ['-i', $fix],
        'Mark NA_fuse in fix matrix',
    );

    my $free_na_fuse = default_na_fuse_name($free);
    my $fix_na_fuse  = default_na_fuse_name($fix);

    die "[ERROR] Expected output not found: $free_na_fuse\n" unless -e $free_na_fuse;
    die "[ERROR] Expected output not found: $fix_na_fuse\n"  unless -e $fix_na_fuse;

    run_perl_script(
        'confirm_na_structure.pl',
        [
            '--fix', $fix_na_fuse,
            '--free', $free_na_fuse,
            (defined $species && $species ne '' ? ('--species_tree', $species) : ()),
            '-o', $final,
        ],
        'Classify NA_struct and NA_topo',
    );

    for my $must_exist (
        "$final.fix.na_classified.txt",
        "$final.free.na_classified.txt",
    ) {
        die "[ERROR] Expected output not found after finalize mode: $must_exist\n" unless -e $must_exist;
    }

    print STDERR "[INFO] Finalize mode completed successfully\n";
    print STDERR "[INFO] Generated: $final.fix.na_classified.txt\n";
    print STDERR "[INFO] Generated: $final.free.na_classified.txt\n";
    if (defined $species && $species ne '') {
        print STDERR "[INFO] Generated: $final.support_b.txt\n";
    }
}

# ------------------------------------------------------------------------------
# finalize_fix mode
# ------------------------------------------------------------------------------
sub run_finalize_fix_mode {
    my %arg = @_;

    my $fix   = $arg{fix};
    my $final = $arg{final_label};

    die "[ERROR] --fix is required in --mode finalize_fix\n" unless defined $fix && $fix ne '';
    die "[ERROR] --final_label is required in --mode finalize_fix\n" unless defined $final && $final ne '';

    die "[ERROR] FIX matrix file not found: $fix\n" unless -e $fix;
    validate_label($final, '--final_label');

    print STDERR "[INFO] Running SplitAligner finalize_fix mode\n";
    print STDERR "[INFO] FIX matrix  : $fix\n";

    run_perl_script(
        'extract_na_fuse.pl',
        ['-i', $fix],
        'Mark NA_fuse in fix matrix',
    );

    my $fix_na_fuse = default_na_fuse_name($fix);
    die "[ERROR] Expected output not found: $fix_na_fuse\n" unless -e $fix_na_fuse;

    run_perl_script(
        'classify_fix_missingness.pl',
        [
            '--fix', $fix_na_fuse,
            '--output', "$final.fix.na_classified.txt",
        ],
        'Classify remaining NA as NA_struct in fix-only mode',
    );

    die "[ERROR] Expected output not found after finalize_fix mode: $final.fix.na_classified.txt\n"
        unless -e "$final.fix.na_classified.txt";

    print STDERR "[INFO] Finalize_fix mode completed successfully\n";
    print STDERR "[INFO] Generated: $final.fix.na_classified.txt\n";
}

# ------------------------------------------------------------------------------
# Helpers
# ------------------------------------------------------------------------------
sub run_perl_script {
    my ($script_name, $args_ref, $desc) = @_;
    my $script_path = File::Spec->catfile($scripts_dir, $script_name);

    print STDERR "[INFO] $desc ...\n";
    my @cmd = ('perl', $script_path, @{$args_ref});

    my $exit_code = system(@cmd);
    if ($exit_code != 0) {
        my $status = $exit_code >> 8;
        die "[ERROR] Step failed: $desc\n[ERROR] Command: @cmd\n[ERROR] Exit status: $status\n";
    }
}

sub default_na_fuse_name {
    my ($path) = @_;
    if ($path =~ /\.txt$/) {
        (my $out = $path) =~ s/\.txt$/.na_fuse.txt/;
        return $out;
    }
    return $path . '.na_fuse.txt';
}

sub validate_label {
    my ($value, $arg_name) = @_;
    die "[ERROR] $arg_name cannot be empty\n" unless defined $value && $value ne '';
    if ($value =~ m{[\\/]}) {
        die "[ERROR] $arg_name should be a simple label/prefix, not a path: $value\n";
    }
}

sub usage {
    return <<'USAGE';
Usage:
  SplitAligner.pl --mode matrix --species <species_tree.nwk> --gene <gene_trees.nwk> --label <label>
  SplitAligner.pl --mode finalize --free <free.matrix_with_fuse.txt> --fix <fix.matrix_with_fuse.txt> --final_label <prefix> [--species_tree <species_tree.forSplit.nwk>]
  SplitAligner.pl --mode finalize_fix --fix <fix.matrix_with_fuse.txt> --final_label <prefix>

Modes:
  matrix
    Build branch matrices from a species tree and a gene-tree file.

    Required:
      --species    Species tree in Newick format
      --gene       Gene-tree file in Newick format
      --label      Output label/prefix (e.g. free, fix)

    Outputs:
      species_tree.forSplit.nwk
      species_tree.FigTree.tre
      species_tree.splits.txt
      species_tree.branch_map.txt
      <label>_splits/
      <label>_split_branch_label/
      <label>.matrix_no_fuse.txt
      <label>.matrix_with_fuse.txt

  finalize
    Apply NA_fuse, then classify NA_struct / NA_topo using shared genes between
    fixed-tree and free-tree matrices. If --species_tree is provided, also
    compute branch-wise Support and write an annotated species tree.

    Required:
      --free         FREE matrix_with_fuse file
      --fix          FIX matrix_with_fuse file
      --final_label  Output prefix for final classified matrices

    Optional:
      --species_tree Species tree in forSplit format for Support annotation

    Outputs:
      <free>.na_fuse.txt
      <fix>.na_fuse.txt
      <final_label>.fix.na_classified.txt
      <final_label>.free.na_classified.txt
      <final_label>.support_b.txt
      <species_prefix>.support_b.nwk

  finalize_fix
    Apply NA_fuse to a fixed-topology matrix, then classify all remaining NA
    cells as NA_struct. This mode is intended for fix-only analyses when no
    free-topology comparator is available.

    Required:
      --fix          FIX matrix_with_fuse file
      --final_label  Output prefix for the classified FIX matrix

    Outputs:
      <fix>.na_fuse.txt
      <final_label>.fix.na_classified.txt

Notes:
  - Helper scripts are expected in: scripts/
  - The final NA-classified matrices are generated for genes shared between the
    fixed-tree and free-tree inputs.
  - The program stops with an error if no shared genes are found.
USAGE
}
