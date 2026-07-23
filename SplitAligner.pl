#!/usr/bin/env perl
# ==============================================================================
# Script:      SplitAligner.pl
# Author:      Jiaqi Wu (wujiaqi@hiroshima-u.ac.jp)
# Description: Transactional controller for SplitAligner matrix construction
#              and evidence-based missingness classification.
# ==============================================================================

use strict;
use warnings;
use File::Basename qw(basename dirname);
use File::Copy qw(copy);
use File::Path qw(make_path remove_tree);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($RealBin);
use Getopt::Long qw(GetOptions);
use lib "$RealBin/lib";
use SplitAligner::Newick qw(
    parse_gene_tree_file
    parse_labeled_species_tree_file
    parse_species_tree_file
    retained_primitive_axis
    validate_gene_taxa
);
use SplitAligner::Provenance qw(
    axis_sha256
    file_sha256
    read_json_file
    verify_axis_equal
    write_json_file
);
use SplitAligner::CoordinateState qw(
    read_matrix_table
    read_state_matrix
    state_code_legend
    state_schema_version
    validate_state_for_matrix
);
use SplitAligner::IOSafety qw(
    assert_io_namespaces_disjoint
    build_destination_record
    build_input_record
    preflight_publication_devices
    verify_input_snapshots
    verify_publication_device_guard
);
use SplitAligner::OutputOwnership qw(
    attach_output_ownership
    build_publication_plan
    object_matches_snapshot
    path_exists
    preflight_output_ownership
    verify_publication_plan
    verify_staged_replacement
);
use SplitAligner::TextIO qw(
    abs_path_utf8
    configure_utf8_stdio
    decode_argv_utf8
    decode_filesystem_path_utf8
    getcwd_utf8
    read_utf8_lines
);

decode_argv_utf8();
configure_utf8_stdio();

my ($mode, $species_tree, $gene_trees, $label, $free_matrix, $fix_matrix,
    $final_label, $free_manifest, $fix_manifest, $force, $help);

GetOptions(
    'mode=s'          => \$mode,
    'species=s'       => \$species_tree,
    'gene=s'          => \$gene_trees,
    'label=s'         => \$label,
    'free=s'          => \$free_matrix,
    'fix=s'           => \$fix_matrix,
    'final_label=s'   => \$final_label,
    'species_tree=s'  => \$species_tree,
    'free_manifest=s' => \$free_manifest,
    'fix_manifest=s'  => \$fix_manifest,
    'force!'          => \$force,
    'help|h'          => \$help,
) or die usage();

if ($help) {
    print usage();
    exit 0;
}
die usage() unless defined $mode;

my $program_root = decode_filesystem_path_utf8($RealBin, 'SplitAligner program directory');
my $scripts_dir = File::Spec->catdir($program_root, 'scripts');
for my $required_script (qw(
    label_species_tree.pl tree_to_splits.pl split_branch_label.pl
    generate_branch_matrix.pl extract_na_fuse.pl confirm_na_structure.pl
    classify_fix_missingness.pl
)) {
    my $path = File::Spec->catfile($scripts_dir, $required_script);
    die "[ERROR] Required helper script not found: $path\n" unless -e $path;
}

if ($mode eq 'matrix') {
    run_matrix_mode(
        species => $species_tree,
        gene    => $gene_trees,
        label   => $label,
        force   => $force,
    );
} elsif ($mode eq 'finalize') {
    run_finalize_mode(
        free          => $free_matrix,
        fix           => $fix_matrix,
        final_label   => $final_label,
        species_tree  => $species_tree,
        free_manifest => $free_manifest,
        fix_manifest  => $fix_manifest,
        force         => $force,
    );
} elsif ($mode eq 'finalize_fix') {
    run_finalize_fix_mode(
        fix          => $fix_matrix,
        final_label  => $final_label,
        fix_manifest => $fix_manifest,
        force        => $force,
    );
} else {
    die "[ERROR] Unsupported --mode: $mode\n" . usage();
}

exit 0;

sub run_matrix_mode {
    my %arg = @_;
    my ($species, $gene, $label0) = @arg{qw(species gene label)};

    die "[ERROR] --species is required in --mode matrix\n" unless defined $species && $species ne '';
    die "[ERROR] --gene is required in --mode matrix\n" unless defined $gene && $gene ne '';
    die "[ERROR] --label is required in --mode matrix\n" unless defined $label0 && $label0 ne '';
    validate_matrix_label($label0);

    my $cwd = getcwd_utf8();
    my @inputs = (
        build_input_record(role => 'matrix species tree', path => $species),
        build_input_record(role => 'matrix gene trees', path => $gene),
    );
    my $species_abs = $inputs[0]{resolved_path};
    my $gene_abs = $inputs[1]{resolved_path};
    my @common_names = qw(
        species_tree.forSplit.nwk
        species_tree.FigTree.tre
        species_tree.splits.txt
        species_tree.branch_map.txt
        species_tree.primitive_axis.tsv
    );
    my @label_owned = (
        "${label0}_splits",
        "${label0}_split_branch_label",
        "$label0.gene_id_map.tsv",
        "$label0.primitive_axis.tsv",
        "$label0.primitive_state.tsv",
        "$label0.matrix_no_fuse.txt",
        "$label0.matrix_with_fuse.txt",
        "$label0.run_manifest.json",
    );
    my %common_role = (
        'species_tree.forSplit.nwk'      => 'matrix labeled species tree',
        'species_tree.FigTree.tre'       => 'matrix FigTree species tree',
        'species_tree.splits.txt'        => 'matrix species split axis',
        'species_tree.branch_map.txt'    => 'matrix species branch map',
        'species_tree.primitive_axis.tsv'=> 'matrix species primitive axis',
    );
    my @items;
    push @items, map {
        {
            source      => $_,
            destination => File::Spec->catfile($cwd, $_),
            role        => $common_role{$_},
            type        => 'file',
            shared      => 1,
        }
    } @common_names;
    push @items,
        {
            source      => "${label0}_splits",
            destination => File::Spec->catfile($cwd, "${label0}_splits"),
            role        => 'matrix gene split directory',
            type        => 'directory',
        },
        {
            source      => "${label0}_split_branch_label",
            destination => File::Spec->catfile($cwd, "${label0}_split_branch_label"),
            role        => 'matrix branch-label directory',
            type        => 'directory',
        };
    push @items, map {
        {
            source         => $_,
            destination    => File::Spec->catfile($cwd, $_),
            role           => "matrix output $_",
            type           => 'file',
            owner_manifest => ($_ eq "$label0.run_manifest.json" ? 1 : 0),
        }
    } @label_owned[2 .. $#label_owned];
    my @destinations = map {
        build_destination_record(
            role => $_->{role}, path => $_->{destination}, type => $_->{type})
    } @items;
    assert_io_namespaces_disjoint(
        inputs       => \@inputs,
        destinations => \@destinations,
    );
    my $publication_device_guard = preflight_publication_devices(
        mode             => 'matrix',
        transaction_root => $cwd,
        destinations     => \@destinations,
    );
    my $owner_manifest_path = File::Spec->catfile(
        $cwd, "$label0.run_manifest.json");
    my $ownership_authorization = preflight_output_ownership(
        mode                => 'matrix',
        owner_label         => $label0,
        owner_manifest_path => $owner_manifest_path,
        items               => \@items,
        force               => $arg{force},
    );

    print STDERR "[INFO] Preflighting species and gene trees before creating outputs.\n";
    my $species_root = parse_species_tree_file($species_abs);
    my $gene_records = parse_gene_tree_file($gene_abs);
    validate_gene_taxa($gene_records, $species_root, $gene_abs);

    my $tmp = create_workdir($cwd, $label0);
    my $ok = eval {
        in_directory($tmp, sub {
            run_perl_script('label_species_tree.pl', [
                '-i', $species_abs, '-o', 'species_tree', '-l', 'species_tree',
            ], 'Label species tree');
            run_perl_script('tree_to_splits.pl', [
                '-i', 'species_tree.forSplit.nwk', '-m', 'species', '-o', 'species_tree',
            ], 'Generate species split axis and primitive ledger');
            run_perl_script('tree_to_splits.pl', [
                '-i', $gene_abs, '-m', 'gene', '-o', $label0,
            ], 'Generate gene splits and identity map');
            run_perl_script('split_branch_label.pl', [
                '-i', "${label0}_splits",
                '-j', 'species_tree.splits.txt',
                '-a', 'species_tree.primitive_axis.tsv',
                '-g', "$label0.gene_id_map.tsv",
                '-o', $label0,
            ], 'Map gene splits onto species axis');
            run_perl_script('generate_branch_matrix.pl', [
                '-i', "${label0}_split_branch_label",
                '-o', $label0,
                '-a', 'species_tree.primitive_axis.tsv',
                '-g', "$label0.gene_id_map.tsv",
            ], 'Generate branch matrices from the complete primitive ledger');

            copy('species_tree.primitive_axis.tsv', "$label0.primitive_axis.tsv")
                or die "[ERROR] Cannot create label-specific primitive ledger: $!\n";

            for my $required (
                'species_tree.forSplit.nwk', 'species_tree.FigTree.tre',
                'species_tree.splits.txt', 'species_tree.branch_map.txt',
                'species_tree.primitive_axis.tsv', "$label0.gene_id_map.tsv",
                "$label0.primitive_state.tsv",
                "$label0.matrix_no_fuse.txt", "$label0.matrix_with_fuse.txt",
            ) {
                die "[ERROR] Expected matrix-mode output not found: $required\n" unless -e $required;
            }

            my $axis = read_axis_from_manifest_sidecar("$label0.primitive_axis.tsv");
            my $gene_identity = read_gene_identity_rows("$label0.gene_id_map.tsv");
            my @axis_ids = map { $_->{branch_id} } @{$axis};
            my $state = read_state_matrix("$label0.primitive_state.tsv");
            my $matrix = read_matrix_table("$label0.matrix_with_fuse.txt");
            validate_state_for_matrix(
                state  => $state,
                matrix => $matrix,
                axis   => \@axis_ids,
                role   => "matrix run '$label0'",
            );
            my @manifest_genes = map { $_->{gene_id} } @{$gene_identity};
            verify_text_array_equal(
                $state->{genes}, \@manifest_genes,
                "state genes for '$label0'", "gene identity map for '$label0'",
            );
            my $matrix_with_fuse_sha = file_sha256("$label0.matrix_with_fuse.txt");
            my $manifest = {
                schema_version => 'SplitAligner-run-manifest-v3',
                mode           => 'matrix',
                label          => $label0,
                program        => implementation_identity(),
                axis           => {
                    serialization => 'UTF-8/LF; one row per coordinate: B_alias<TAB>canonical_split_key<LF>',
                    canonical_key_schema => 'SplitAligner-canonical-split-key-v2',
                    coordinate_count => scalar(@{$axis}),
                    sha256          => axis_sha256($axis),
                    coordinates     => $axis,
                    ledger_file     => "$label0.primitive_axis.tsv",
                },
                inputs => {
                    species_tree => {
                        path_as_provided => $species,
                        sha256           => file_sha256($species_abs),
                    },
                    gene_trees => {
                        path_as_provided => $gene,
                        sha256           => file_sha256($gene_abs),
                        record_count     => scalar(@{$gene_records}),
                    },
                },
                gene_identity => $gene_identity,
                coordinate_state => {
                    schema_version       => state_schema_version(),
                    serialization        => 'UTF-8/LF; header gene<TAB>B1...; one S/D/F/U code per primitive cell',
                    code_legend          => state_code_legend(),
                    filename             => "$label0.primitive_state.tsv",
                    sha256               => file_sha256("$label0.primitive_state.tsv"),
                    primitive_axis_sha256 => axis_sha256($axis),
                    gene_record_count    => $state->{gene_count},
                    cell_count           => $state->{cell_count},
                    matrix_with_fuse_sha256 => $matrix_with_fuse_sha,
                },
                species_artifacts => {
                    labeled_tree_sha256 => file_sha256('species_tree.forSplit.nwk'),
                    branch_map_sha256   => file_sha256('species_tree.branch_map.txt'),
                    splits_sha256       => file_sha256('species_tree.splits.txt'),
                },
                outputs => {
                    matrix_no_fuse => {
                        filename => "$label0.matrix_no_fuse.txt",
                        sha256   => file_sha256("$label0.matrix_no_fuse.txt"),
                    },
                    matrix_with_fuse => {
                        filename => "$label0.matrix_with_fuse.txt",
                        sha256   => $matrix_with_fuse_sha,
                    },
                    primitive_state => {
                        filename => "$label0.primitive_state.tsv",
                        sha256   => file_sha256("$label0.primitive_state.tsv"),
                    },
                },
            };
            attach_output_ownership(
                $manifest,
                mode                => 'matrix',
                owner_label         => $label0,
                owner_manifest_path => $owner_manifest_path,
                source_root         => $tmp,
                items               => \@items,
            );
            write_json_file("$label0.run_manifest.json", $manifest);
        });
        publish_transaction(
            $tmp, \@items, $arg{force}, \@inputs, \@destinations,
            $ownership_authorization, $publication_device_guard);
        1;
    };
    my $error = $@;
    remove_tree($tmp) if -d $tmp;
    die $error unless $ok;

    print STDERR "[INFO] Matrix mode completed successfully for label: $label0\n";
    print STDERR "[INFO] Manifest: $label0.run_manifest.json\n";
}

sub run_finalize_mode {
    my %arg = @_;
    my ($free, $fix, $final, $species) = @arg{qw(free fix final_label species_tree)};
    die "[ERROR] --free is required in --mode finalize\n" unless defined $free && $free ne '';
    die "[ERROR] --fix is required in --mode finalize\n" unless defined $fix && $fix ne '';
    die "[ERROR] --final_label is required in --mode finalize\n" unless defined $final && $final ne '';
    validate_label($final, '--final_label');

    my $free_input = build_input_record(role => 'FREE matrix', path => $free);
    my $fix_input = build_input_record(role => 'FIX matrix', path => $fix);
    my $free_abs = $free_input->{resolved_path};
    my $fix_abs = $fix_input->{resolved_path};
    die "[ERROR] FREE and FIX matrix basenames must differ for transactional finalize mode.\n"
        if basename($free_abs) eq basename($fix_abs);

    my ($free_meta, $free_manifest_path, $free_state_abs,
        $free_manifest_input, $free_state_input) = load_and_verify_matrix_manifest(
        $free_abs, $arg{free_manifest}, 'FREE');
    my ($fix_meta, $fix_manifest_path, $fix_state_abs,
        $fix_manifest_input, $fix_state_input) = load_and_verify_matrix_manifest(
        $fix_abs, $arg{fix_manifest}, 'FIX');
    die "[ERROR] FREE and FIX state-sidecar basenames must differ for transactional finalize mode.\n"
        if basename($free_state_abs) eq basename($fix_state_abs);
    verify_axis_equal(
        $fix_meta->{axis}{coordinates}, $free_meta->{axis}{coordinates},
        "FIX manifest $fix_manifest_path", "FREE manifest $free_manifest_path",
    );

    my @inputs = (
        $free_input, $fix_input,
        $free_manifest_input, $fix_manifest_input,
        $free_state_input, $fix_state_input,
    );
    my ($species_abs, $branch_map_abs, $support_tree_name);
    if (defined $species && $species ne '') {
        my $species_input = build_input_record(
            role => 'Support species tree', path => $species);
        $species_abs = $species_input->{resolved_path};
        my (undef, $support_root) = parse_labeled_species_tree_file($species_abs);
        my $support_axis = retained_primitive_axis($support_root);
        verify_axis_equal(
            $fix_meta->{axis}{coordinates}, $support_axis,
            "FIX manifest $fix_manifest_path", "Support tree $species_abs",
        );
        $branch_map_abs = derive_existing_branch_map($species_abs);
        my $branch_map_input = build_input_record(
            role => 'Support branch map', path => $branch_map_abs);
        $branch_map_abs = $branch_map_input->{resolved_path};
        for my $pair (
            [$free_meta, 'FREE'],
            [$fix_meta, 'FIX'],
        ) {
            my ($meta, $role) = @{$pair};
            my $expected_tree = $meta->{species_artifacts}{labeled_tree_sha256} // '';
            my $expected_map  = $meta->{species_artifacts}{branch_map_sha256} // '';
            die "[ERROR] Support tree SHA256 does not match the $role matrix run manifest.\n"
                unless $expected_tree ne '' && file_sha256($species_abs) eq $expected_tree;
            die "[ERROR] Support branch-map SHA256 does not match the $role matrix run manifest.\n"
                unless $expected_map ne '' && file_sha256($branch_map_abs) eq $expected_map;
        }
        $support_tree_name = support_tree_output_name(basename($species_abs));
        push @inputs, $species_input, $branch_map_input;
    }

    my $cwd = getcwd_utf8();
    my $free_na_name = default_na_fuse_name(basename($free_abs));
    my $fix_na_name = default_na_fuse_name(basename($fix_abs));
    my @items = (
        {
            source      => $free_na_name,
            destination => default_na_fuse_name($free_abs),
            role        => 'FREE matrix NA_fuse output',
            type        => 'file',
            shared      => 1,
        },
        {
            source      => $fix_na_name,
            destination => default_na_fuse_name($fix_abs),
            role        => 'FIX matrix NA_fuse output',
            type        => 'file',
            shared      => 1,
        },
        {
            source      => "$final.fix.na_classified.txt",
            destination => File::Spec->catfile($cwd, "$final.fix.na_classified.txt"),
            role        => 'final FIX classified matrix',
            type        => 'file',
        },
        {
            source      => "$final.free.na_classified.txt",
            destination => File::Spec->catfile($cwd, "$final.free.na_classified.txt"),
            role        => 'final FREE classified matrix',
            type        => 'file',
        },
    );
    if ($species_abs) {
        push @items,
            {
                source      => "$final.support_b.txt",
                destination => File::Spec->catfile($cwd, "$final.support_b.txt"),
                role        => 'Support table',
                type        => 'file',
            },
            {
                source      => $support_tree_name,
                destination => File::Spec->catfile($cwd, $support_tree_name),
                role        => 'Support-annotated species tree',
                type        => 'file',
                shared      => 1,
            };
    }
    my $owner_manifest_path = File::Spec->catfile(
        $cwd, "$final.finalize_manifest.json");
    push @items, {
        source         => "$final.finalize_manifest.json",
        destination    => $owner_manifest_path,
        role           => 'finalize manifest',
        type           => 'file',
        owner_manifest => 1,
    };
    my @destinations = map {
        build_destination_record(
            role => $_->{role}, path => $_->{destination}, type => $_->{type},
        )
    } @items;
    assert_io_namespaces_disjoint(
        inputs       => \@inputs,
        destinations => \@destinations,
    );
    my $publication_device_guard = preflight_publication_devices(
        mode             => 'finalize',
        transaction_root => $cwd,
        destinations     => \@destinations,
    );
    my $ownership_authorization = preflight_output_ownership(
        mode                => 'finalize',
        owner_label         => $final,
        owner_manifest_path => $owner_manifest_path,
        items               => \@items,
        force               => $arg{force},
    );

    my $tmp = create_workdir($cwd, $final);
    my $ok = eval {
        copy($free_abs, File::Spec->catfile($tmp, basename($free_abs)))
            or die "[ERROR] Cannot stage FREE matrix: $!\n";
        copy($fix_abs, File::Spec->catfile($tmp, basename($fix_abs)))
            or die "[ERROR] Cannot stage FIX matrix: $!\n";
        copy($free_state_abs, File::Spec->catfile($tmp, basename($free_state_abs)))
            or die "[ERROR] Cannot stage FREE coordinate state: $!\n";
        copy($fix_state_abs, File::Spec->catfile($tmp, basename($fix_state_abs)))
            or die "[ERROR] Cannot stage FIX coordinate state: $!\n";
        if ($species_abs) {
            copy($species_abs, File::Spec->catfile($tmp, basename($species_abs)))
                or die "[ERROR] Cannot stage species tree: $!\n";
            copy($branch_map_abs, File::Spec->catfile($tmp, basename($branch_map_abs)))
                or die "[ERROR] Cannot stage species branch map: $!\n";
        }

        in_directory($tmp, sub {
            run_perl_script('extract_na_fuse.pl', [
                '-i', basename($free_abs), '-o', $free_na_name,
                '--state', basename($free_state_abs),
            ], 'Mark NA_fuse in free matrix');
            run_perl_script('extract_na_fuse.pl', [
                '-i', basename($fix_abs), '-o', $fix_na_name,
                '--state', basename($fix_state_abs),
            ], 'Mark NA_fuse in fix matrix');
            run_perl_script('confirm_na_structure.pl', [
                '--fix', $fix_na_name,
                '--free', $free_na_name,
                '--fix_state', basename($fix_state_abs),
                '--free_state', basename($free_state_abs),
                ($species_abs ? ('--species_tree', basename($species_abs)) : ()),
                '-o', $final,
            ], 'Classify NA_struct and NA_topo');

            my @output_names = (
                $free_na_name, $fix_na_name,
                "$final.fix.na_classified.txt", "$final.free.na_classified.txt",
            );
            push @output_names, "$final.support_b.txt", $support_tree_name if $species_abs;
            my %output_hash = map {
                $_ => { filename => $_, sha256 => file_sha256($_) }
            } @output_names;
            my $manifest = {
                schema_version => 'SplitAligner-finalize-manifest-v2',
                mode           => 'finalize',
                final_label    => $final,
                program        => implementation_identity(),
                primitive_axis_sha256 => $fix_meta->{axis}{sha256},
                inputs => {
                    free_matrix_sha256   => file_sha256(basename($free_abs)),
                    fix_matrix_sha256    => file_sha256(basename($fix_abs)),
                    free_manifest_sha256 => file_sha256($free_manifest_path),
                    fix_manifest_sha256  => file_sha256($fix_manifest_path),
                    free_state_sha256    => file_sha256(basename($free_state_abs)),
                    fix_state_sha256     => file_sha256(basename($fix_state_abs)),
                    state_schema_version => state_schema_version(),
                },
                outputs => \%output_hash,
            };
            attach_output_ownership(
                $manifest,
                mode                => 'finalize',
                owner_label         => $final,
                owner_manifest_path => $owner_manifest_path,
                source_root         => $tmp,
                items               => \@items,
            );
            write_json_file("$final.finalize_manifest.json", $manifest);
        });
        publish_transaction(
            $tmp, \@items, $arg{force}, \@inputs, \@destinations,
            $ownership_authorization, $publication_device_guard);
        1;
    };
    my $error = $@;
    remove_tree($tmp) if -d $tmp;
    die $error unless $ok;
    print STDERR "[INFO] Finalize mode completed successfully.\n";
}

sub run_finalize_fix_mode {
    my %arg = @_;
    my ($fix, $final) = @arg{qw(fix final_label)};
    die "[ERROR] --fix is required in --mode finalize_fix\n" unless defined $fix && $fix ne '';
    die "[ERROR] --final_label is required in --mode finalize_fix\n" unless defined $final && $final ne '';
    validate_label($final, '--final_label');

    my $fix_input = build_input_record(role => 'FIX matrix', path => $fix);
    my $fix_abs = $fix_input->{resolved_path};
    my ($fix_meta, $fix_manifest_path, $fix_state_abs,
        $fix_manifest_input, $fix_state_input) = load_and_verify_matrix_manifest(
        $fix_abs, $arg{fix_manifest}, 'FIX');
    my @inputs = ($fix_input, $fix_manifest_input, $fix_state_input);

    my $cwd = getcwd_utf8();
    my $fix_na_name = default_na_fuse_name(basename($fix_abs));
    my $owner_manifest_path = File::Spec->catfile(
        $cwd, "$final.finalize_manifest.json");
    my @items = (
        {
            source      => $fix_na_name,
            destination => default_na_fuse_name($fix_abs),
            role        => 'FIX matrix NA_fuse output',
            type        => 'file',
            shared      => 1,
        },
        {
            source      => "$final.fix.na_classified.txt",
            destination => File::Spec->catfile($cwd, "$final.fix.na_classified.txt"),
            role        => 'final FIX classified matrix',
            type        => 'file',
        },
        {
            source         => "$final.finalize_manifest.json",
            destination    => $owner_manifest_path,
            role           => 'finalize_fix manifest',
            type           => 'file',
            owner_manifest => 1,
        },
    );
    my @destinations = map {
        build_destination_record(
            role => $_->{role}, path => $_->{destination}, type => $_->{type},
        )
    } @items;
    assert_io_namespaces_disjoint(
        inputs       => \@inputs,
        destinations => \@destinations,
    );
    my $publication_device_guard = preflight_publication_devices(
        mode             => 'finalize_fix',
        transaction_root => $cwd,
        destinations     => \@destinations,
    );
    my $ownership_authorization = preflight_output_ownership(
        mode                => 'finalize_fix',
        owner_label         => $final,
        owner_manifest_path => $owner_manifest_path,
        items               => \@items,
        force               => $arg{force},
    );

    my $tmp = create_workdir($cwd, $final);
    my $ok = eval {
        copy($fix_abs, File::Spec->catfile($tmp, basename($fix_abs)))
            or die "[ERROR] Cannot stage FIX matrix: $!\n";
        copy($fix_state_abs, File::Spec->catfile($tmp, basename($fix_state_abs)))
            or die "[ERROR] Cannot stage FIX coordinate state: $!\n";
        in_directory($tmp, sub {
            run_perl_script('extract_na_fuse.pl', [
                '-i', basename($fix_abs), '-o', $fix_na_name,
                '--state', basename($fix_state_abs),
            ], 'Mark NA_fuse in fix matrix');
            run_perl_script('classify_fix_missingness.pl', [
                '--fix', $fix_na_name,
                '--state', basename($fix_state_abs),
                '--output', "$final.fix.na_classified.txt",
            ], 'Classify fixed missingness from explicit coordinate state');
            my $manifest = {
                schema_version => 'SplitAligner-finalize-manifest-v2',
                mode           => 'finalize_fix',
                final_label    => $final,
                program        => implementation_identity(),
                primitive_axis_sha256 => $fix_meta->{axis}{sha256},
                inputs => {
                    fix_matrix_sha256   => file_sha256(basename($fix_abs)),
                    fix_manifest_sha256 => file_sha256($fix_manifest_path),
                    fix_state_sha256    => file_sha256(basename($fix_state_abs)),
                    state_schema_version => state_schema_version(),
                },
                outputs => {
                    $fix_na_name => { filename => $fix_na_name, sha256 => file_sha256($fix_na_name) },
                    "$final.fix.na_classified.txt" => {
                        filename => "$final.fix.na_classified.txt",
                        sha256 => file_sha256("$final.fix.na_classified.txt"),
                    },
                },
            };
            attach_output_ownership(
                $manifest,
                mode                => 'finalize_fix',
                owner_label         => $final,
                owner_manifest_path => $owner_manifest_path,
                source_root         => $tmp,
                items               => \@items,
            );
            write_json_file("$final.finalize_manifest.json", $manifest);
        });
        publish_transaction(
            $tmp, \@items, $arg{force}, \@inputs, \@destinations,
            $ownership_authorization, $publication_device_guard);
        1;
    };
    my $error = $@;
    remove_tree($tmp) if -d $tmp;
    die $error unless $ok;
    print STDERR "[INFO] Finalize_fix mode completed successfully.\n";
}

sub run_perl_script {
    my ($script_name, $args, $description) = @_;
    my $script_path = File::Spec->catfile($scripts_dir, $script_name);
    print STDERR "[INFO] $description ...\n";
    my @command = ($^X, $script_path, @{$args});
    my $status = system(@command);
    if ($status != 0) {
        my $exit = $status == -1 ? -1 : ($status >> 8);
        die "[ERROR] Step failed: $description\n[ERROR] Exit status: $exit\n";
    }
}

sub load_and_verify_matrix_manifest {
    my ($matrix_path, $explicit_manifest, $role) = @_;
    my $manifest_path = $explicit_manifest;
    if (!defined $manifest_path || $manifest_path eq '') {
        my $file = basename($matrix_path);
        die "[ERROR] Cannot infer $role manifest from matrix name '$file'; use --" . lc($role) . "_manifest.\n"
            unless $file =~ /\A(.+)\.matrix_with_fuse\.txt\z/;
        $manifest_path = File::Spec->catfile(dirname($matrix_path), "$1.run_manifest.json");
    }
    my $manifest_input = build_input_record(
        role => "$role run manifest", path => $manifest_path);
    $manifest_path = $manifest_input->{resolved_path};
    my $manifest = read_json_file($manifest_path);
    my $schema = $manifest->{schema_version} // '';
    if ($schema eq 'SplitAligner-run-manifest-v2') {
        die "[ERROR] $role matrix was generated without coordinate-state provenance. Regenerate it with the repaired SplitAligner matrix mode before finalization.\n";
    }
    die "[ERROR] Unsupported $role manifest schema in $manifest_path.\n"
        unless $schema eq 'SplitAligner-run-manifest-v3';
    die "[ERROR] Unsupported canonical split-key schema in $role manifest $manifest_path.\n"
        unless ($manifest->{axis}{canonical_key_schema} // '')
            eq 'SplitAligner-canonical-split-key-v2';
    die "[ERROR] $role manifest has no ordered primitive coordinate ledger.\n"
        unless ref($manifest->{axis}{coordinates}) eq 'ARRAY' && @{$manifest->{axis}{coordinates}};
    my $computed_axis_sha = axis_sha256($manifest->{axis}{coordinates});
    die "[ERROR] $role manifest axis SHA256 is inconsistent with its coordinate ledger.\n"
        unless $computed_axis_sha eq ($manifest->{axis}{sha256} // '');
    my $expected_matrix_sha = $manifest->{outputs}{matrix_with_fuse}{sha256} // '';
    die "[ERROR] $role matrix SHA256 does not match its run manifest: $matrix_path\n"
        unless $expected_matrix_sha ne '' && file_sha256($matrix_path) eq $expected_matrix_sha;

    my $state_meta = $manifest->{coordinate_state};
    die "[ERROR] $role manifest has no coordinate-state provenance. Regenerate the matrix run.\n"
        unless ref($state_meta) eq 'HASH';
    die "[ERROR] Unsupported $role coordinate-state schema in $manifest_path.\n"
        unless ($state_meta->{schema_version} // '') eq state_schema_version();
    die "[ERROR] $role coordinate-state axis SHA256 does not match the run axis.\n"
        unless ($state_meta->{primitive_axis_sha256} // '') eq $computed_axis_sha;
    die "[ERROR] $role coordinate state is not bound to the supplied matrix.\n"
        unless ($state_meta->{matrix_with_fuse_sha256} // '') eq $expected_matrix_sha;

    my $state_filename = $state_meta->{filename} // '';
    die "[ERROR] Unsafe or missing $role coordinate-state filename in $manifest_path.\n"
        unless $state_filename ne '' && basename($state_filename) eq $state_filename
            && $state_filename !~ /[\\\/\x00-\x1f]/;
    my $state_path = File::Spec->catfile(dirname($manifest_path), $state_filename);
    my $state_input = build_input_record(
        role => "$role coordinate-state sidecar", path => $state_path);
    $state_path = $state_input->{resolved_path};
    die "[ERROR] $role coordinate-state SHA256 does not match its run manifest.\n"
        unless ($state_meta->{sha256} // '') ne ''
            && file_sha256($state_path) eq $state_meta->{sha256};
    die "[ERROR] $role coordinate-state output record is inconsistent with coordinate_state metadata.\n"
        unless ref($manifest->{outputs}{primitive_state}) eq 'HASH'
            && ($manifest->{outputs}{primitive_state}{filename} // '') eq $state_filename
            && ($manifest->{outputs}{primitive_state}{sha256} // '') eq $state_meta->{sha256};

    my @axis_ids = map { $_->{branch_id} } @{$manifest->{axis}{coordinates}};
    my $state = read_state_matrix($state_path);
    my $matrix = read_matrix_table($matrix_path);
    validate_state_for_matrix(
        state  => $state,
        matrix => $matrix,
        axis   => \@axis_ids,
        role   => $role,
    );
    die "[ERROR] $role coordinate-state gene count does not match its run manifest.\n"
        unless 0 + ($state_meta->{gene_record_count} // -1) == $state->{gene_count};
    die "[ERROR] $role coordinate-state cell count does not match its run manifest.\n"
        unless 0 + ($state_meta->{cell_count} // -1) == $state->{cell_count};

    die "[ERROR] $role manifest has no gene identity ledger.\n"
        unless ref($manifest->{gene_identity}) eq 'ARRAY';
    my @manifest_genes = map { $_->{gene_id} // '' } @{$manifest->{gene_identity}};
    verify_text_array_equal(
        $state->{genes}, \@manifest_genes,
        "$role coordinate-state genes", "$role manifest gene ledger",
    );

    return ($manifest, $manifest_path, $state_path, $manifest_input, $state_input);
}

sub verify_text_array_equal {
    my ($left, $right, $left_name, $right_name) = @_;
    die "[ERROR] Invalid comparison between $left_name and $right_name.\n"
        unless ref($left) eq 'ARRAY' && ref($right) eq 'ARRAY';
    my $max = @{$left} > @{$right} ? @{$left} : @{$right};
    for my $i (0 .. $max - 1) {
        my $l = $i < @{$left} ? $left->[$i] : '<missing>';
        my $r = $i < @{$right} ? $right->[$i] : '<missing>';
        next if $l eq $r;
        die "[ERROR] $left_name differs from $right_name at position " . ($i + 1)
            . ": '$l' vs '$r'.\n";
    }
    return 1;
}

sub read_axis_from_manifest_sidecar {
    my ($path) = @_;
    require SplitAligner::Provenance;
    return SplitAligner::Provenance::read_axis_ledger($path);
}

sub read_gene_identity_rows {
    my ($path) = @_;
    my $lines = read_utf8_lines($path);
    die "[ERROR] Gene identity map is empty: $path\n" unless @{$lines};
    die "[ERROR] Unexpected gene identity map header in $path.\n"
        unless $lines->[0] eq "storage_key\tgene_id\tinput_line";
    my @rows;
    for my $index (1 .. $#{$lines}) {
        my $line = $lines->[$index];
        next if $line eq '';
        my ($storage_key, $gene_id, $input_line) = split(/\t/, $line, -1);
        push @rows, {
            storage_key => $storage_key,
            gene_id     => $gene_id,
            input_line  => 0 + $input_line,
        };
    }
    return \@rows;
}

sub implementation_identity {
    my %files;
    for my $relative (
        'SplitAligner.pl',
        'lib/SplitAligner/Newick.pm',
        'lib/SplitAligner/Provenance.pm',
        'lib/SplitAligner/TextIO.pm',
        'lib/SplitAligner/CoordinateState.pm',
        'lib/SplitAligner/IOSafety.pm',
        'lib/SplitAligner/OutputOwnership.pm',
        map { "scripts/$_" } qw(
            label_species_tree.pl tree_to_splits.pl split_branch_label.pl
            generate_branch_matrix.pl extract_na_fuse.pl confirm_na_structure.pl
            classify_fix_missingness.pl
        ),
    ) {
        my $path = File::Spec->catfile($program_root, split('/', $relative));
        $files{$relative} = file_sha256($path);
    }
    return {
        name                => 'SplitAligner',
        frozen_base_version => 'v1.2.0',
        frozen_base_commit  => '4337e2cb62ac24238246893e85e3d2a5a51c8a7a',
        repair_contract     => 'RECERT-BlockerRepair-20260716+UTF8PathClosure-20260717+IOAliasClosure-20260718+OutputOwnershipClosure-20260718+CrossDevicePreflightClosure-20260718+PublicationParentDeviceClosure-20260718+POSIXInputResolutionClosure-20260719+StrictPOSIXInputLookupClosure-20260719',
        implementation_files_sha256 => \%files,
    };
}

sub derive_existing_branch_map {
    my ($tree_path) = @_;
    my $dir = dirname($tree_path);
    my $file = basename($tree_path);
    my @candidate;
    if ($file =~ /\A(.+)\.forSplit\.nwk\z/) {
        push @candidate, File::Spec->catfile($dir, "$1.branch_map.txt");
    }
    push @candidate, File::Spec->catfile($dir, 'species_tree.branch_map.txt');
    for my $path (@candidate) {
        # Return the derived path as provided. build_input_record() performs
        # OS-level resolution and retains this path for retarget detection.
        return $path if -e $path || -l $path;
    }
    die "[ERROR] A branch map matching the Support tree is required but was not found beside $tree_path.\n";
}

sub support_tree_output_name {
    my ($name) = @_;
    $name =~ s/\.forSplit\.nwk\z/.support_b.nwk/ and return $name;
    $name =~ s/\.nwk\z/.support_b.nwk/ and return $name;
    return $name . '.support_b.nwk';
}

sub create_workdir {
    my ($parent, $label0) = @_;
    (my $safe = $label0) =~ s/[^A-Za-z0-9_.-]+/_/g;
    my $path = tempdir(".splitaligner-$safe-XXXXXX", DIR => $parent, CLEANUP => 0);
    return decode_filesystem_path_utf8($path, 'temporary SplitAligner workspace');
}

sub in_directory {
    my ($directory, $code) = @_;
    my $old = getcwd_utf8();
    chdir $directory or die "[ERROR] Cannot enter temporary work directory $directory: $!\n";
    my $ok = eval { $code->(); 1 };
    my $error = $@;
    chdir $old or die "[ERROR] Cannot return to $old: $!\n";
    die $error unless $ok;
}

sub publish_transaction {
    my ($tmp, $items, $force0, $inputs, $destinations, $authorization,
        $publication_device_guard) = @_;
    verify_input_snapshots($inputs);
    assert_io_namespaces_disjoint(
        inputs       => $inputs,
        destinations => $destinations,
    );
    my $plan = build_publication_plan(
        source_root   => $tmp,
        items         => $items,
        force         => $force0,
        authorization => $authorization,
    );
    verify_publication_plan($plan);
    my @publish = grep { $_->{action} ne 'reuse' } @{$plan->{items}};
    verify_publication_device_guard(
        guard             => $publication_device_guard,
        transaction_root => getcwd_utf8(),
        destinations     => $destinations,
    );
    return unless @publish;

    my $backup = tempdir(
        '.splitaligner-publish-backup-XXXXXX',
        DIR => getcwd_utf8(),
        CLEANUP => 0,
    );
    $backup = decode_filesystem_path_utf8($backup, 'transaction backup workspace');
    my (@backed_up, @installed);
    my $fail_after = publication_test_fail_after();
    my $ok = eval {
        my $index = 0;
        for my $item (@publish) {
            my $destination = $item->{destination};
            if ($item->{action} eq 'replace') {
                die "[ERROR][Output ownership] Authorized destination changed "
                    . "during publication for role '$item->{role}': "
                    . "$destination.\n"
                    unless object_matches_snapshot(
                        $destination, $item->{type}, $item->{existing_snapshot});
                my $backup_path = File::Spec->catfile($backup, sprintf('%06d', ++$index));
                rename($destination, $backup_path)
                    or die "[ERROR] Cannot stage existing output for rollback: $destination: $!\n";
                verify_staged_replacement($item, $backup_path);
                push @backed_up, [$backup_path, $destination];
            } elsif (path_exists($destination)) {
                die "[ERROR][Output ownership] An unowned destination appeared "
                    . "during publication for role '$item->{role}': "
                    . "$destination. It was not replaced.\n";
            }
            rename($item->{source_abs}, $destination)
                or die "[ERROR] Cannot publish $item->{source_abs} to $destination: $!\n";
            die "[ERROR] Published output differs from its generated source "
                . "snapshot for role '$item->{role}': $destination.\n"
                unless object_matches_snapshot(
                    $destination, $item->{type}, $item->{source_snapshot});
            push @installed, [$destination, $item->{type}, $item->{source_snapshot}];
            die "[ERROR][Internal test] Injected late publication failure.\n"
                if defined $fail_after && @installed == $fail_after;
        }
        1;
    };
    my $error = $@;
    if (!$ok) {
        my @rollback_error;
        for my $installed (reverse @installed) {
            my ($path, $type, $snapshot) = @{$installed};
            if (object_matches_snapshot($path, $type, $snapshot)) {
                eval { remove_path($path); 1 }
                    or push @rollback_error,
                        "[ERROR] Cannot remove newly installed output during rollback: $path: $@";
            } else {
                push @rollback_error,
                    "[ERROR] Newly installed output changed during rollback and was preserved: $path\n";
            }
        }
        for my $pair (reverse @backed_up) {
            if (path_exists($pair->[1])) {
                push @rollback_error,
                    "[ERROR] Rollback destination is occupied and was preserved: $pair->[1]\n";
                next;
            }
            rename($pair->[0], $pair->[1])
                or push @rollback_error,
                    "[ERROR] Rollback failed for $pair->[1]: $!\n";
        }
        remove_tree($backup) if -d $backup;
        die $error . join('', @rollback_error);
    }
    remove_tree($backup) if -d $backup;
}

sub publication_test_fail_after {
    return undef unless ($ENV{SPLITALIGNER_INTERNAL_TESTING} // '') eq '1';
    my $value = $ENV{SPLITALIGNER_TEST_FAIL_AFTER_PUBLISH};
    return undef unless defined $value && $value ne '';
    die "[ERROR][Internal test] Invalid publication failpoint.\n"
        unless $value =~ /\A[1-9]\d*\z/;
    return 0 + $value;
}

sub remove_path {
    my ($path) = @_;
    if (-d $path) {
        remove_tree($path);
    } elsif (-e $path) {
        unlink $path or die "[ERROR] Cannot remove $path during rollback: $!\n";
    }
}

sub default_na_fuse_name {
    my ($path) = @_;
    if ($path =~ /\.txt\z/) {
        (my $out = $path) =~ s/\.txt\z/.na_fuse.txt/;
        return $out;
    }
    return $path . '.na_fuse.txt';
}

sub validate_label {
    my ($value, $arg_name) = @_;
    die "[ERROR] $arg_name cannot be empty\n" unless defined $value && $value ne '';
    die "[ERROR] $arg_name must be a simple label without path separators or control characters: $value\n"
        if $value =~ m{[\\/\x00-\x1f]};
}

sub validate_matrix_label {
    my ($value) = @_;
    validate_label($value, '--label');

    # Matrix mode always publishes the fixed common species_tree.* artifacts.
    # The only label-owned output that can coincide with that namespace under
    # the current naming scheme is <label>.primitive_axis.tsv.
    if ($value eq 'species_tree') {
        die "[ERROR] --label 'species_tree' is reserved because it collides with common species-tree artifacts, including species_tree.primitive_axis.tsv. Choose a different matrix label.\n";
    }
}

sub usage {
    return <<'USAGE';
Usage:
  SplitAligner.pl --mode matrix --species <species.nwk> --gene <genes.nwk> --label <label> [--force]
  SplitAligner.pl --mode finalize --free <free.matrix_with_fuse.txt> --fix <fix.matrix_with_fuse.txt> --final_label <prefix> [--species_tree <species_tree.forSplit.nwk>] [--free_manifest <json>] [--fix_manifest <json>] [--force]
  SplitAligner.pl --mode finalize_fix --fix <fix.matrix_with_fuse.txt> --final_label <prefix> [--fix_manifest <json>] [--force]

Safety and provenance:
  - Every input and publication destination is registered before work begins;
    input/output aliases are rejected even with --force.
  - Input bytes are snapshotted at preflight and rechecked immediately before
    transactional publication.
  - Every publication destination's resolved publication parent must use the
    same filesystem device as the current working directory. Existing output
    inode devices are not publication authority; --force cannot override this
    requirement.
  - Matrix runs write <label>.primitive_axis.tsv, <label>.gene_id_map.tsv, and
    <label>.run_manifest.json.
  - finalize and finalize_fix require matching run manifests; legacy fallback is
    intentionally disabled.
  - Every successful mode records a complete output ownership inventory.
    Existing label-owned outputs cause an error by default; --force can replace
    only an intact prior run covered by its matching ownership inventory.
  - Shared outputs are reused without mutation only when identical. Differing,
    unowned, legacy, malformed, or tampered destinations fail closed.
USAGE
}
