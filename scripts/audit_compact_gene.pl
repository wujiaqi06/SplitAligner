#!/usr/bin/env perl
use strict;
use warnings;
use FindBin qw($RealBin);
use lib "$RealBin/../lib";
use Getopt::Long qw(GetOptions);
use File::Basename qw(dirname basename);
use File::Spec;
use File::Temp qw(tempdir);
use Digest::SHA;
use JSON::PP;
use SplitAligner::TextIO qw(decode_argv_utf8 configure_utf8_stdio read_utf8_lines encode_utf8_strict);
use SplitAligner::Provenance qw(read_json_file file_sha256 axis_sha256 verify_axis_equal);
use SplitAligner::Newick qw(parse_newick parse_species_tree_file assign_branch_ids retained_primitive_axis tip_labels);
use SplitAligner::CompactKernel qw(prepare_reference compute_gene materialize_gene construct_cell);
use SplitAligner::CompactPipeline qw(visit_input_lines read_ledger_block);

decode_argv_utf8();
configure_utf8_stdio();
my ($manifest_path, $gene_input, $species_input, $gene_list, @wanted, $help);
GetOptions('manifest=s' => \$manifest_path, 'gene-input=s' => \$gene_input,
    'species-input=s' => \$species_input, 'audit-gene=s@' => \@wanted,
    'audit-gene-file=s' => \$gene_list, 'help' => \$help) or die usage();
if ($help) { print usage(); exit 0; }
die usage() unless defined $manifest_path;
push @wanted, grep { length } @{read_utf8_lines($gene_list)} if defined $gene_list;
die "[ERROR] Select at least one gene using --audit-gene or --audit-gene-file.\n" unless @wanted;
my $meta = read_json_file($manifest_path);
die "[ERROR] Manifest does not describe a compact source ledger.\n"
    unless ($meta->{compact_provenance}{schema_version} // '') eq 'SplitAligner-compact-ledger-v1';
my $base = dirname(File::Spec->rel2abs($manifest_path));
my $tmp = tempdir('splitaligner-audit-XXXXXX', TMPDIR => 1, CLEANUP => 1);
my %files;
for my $kind (qw(index dictionary ledger)) {
    my $info = $meta->{compact_provenance}{files}{$kind};
    die "[ERROR] Invalid compact artifact path.\n"
        unless ref($info) eq 'HASH' && defined($info->{filename})
        && $info->{filename} eq basename($info->{filename}) && $info->{filename} !~ /\A\.{1,2}\z/;
    $files{$kind} = freeze_verified(File::Spec->catfile($base, $info->{filename}),
        File::Spec->catfile($tmp, $kind), $info->{sha256});
}
my $index = read_json_file($files{index});
my $dictionary = read_json_file($files{dictionary});
die "[ERROR] Unsupported compact index or dictionary schema.\n"
    unless ($index->{schema_version} // '') eq 'SplitAligner-compact-ledger-v1'
        && ($dictionary->{schema_version} // '') eq 'SplitAligner-taxon-dictionary-v1';
die "[ERROR] Compact artifact fingerprints disagree.\n"
    unless $index->{ledger_sha256} eq $meta->{compact_provenance}{files}{ledger}{sha256}
        && $index->{dictionary_sha256} eq $meta->{compact_provenance}{files}{dictionary}{sha256};
die "[ERROR] Parser identity differs; use the reader matching this result's parser contract.\n"
    unless file_sha256("$RealBin/../lib/SplitAligner/Newick.pm") eq ($index->{parser_sha256} // '');
for my $kind (qw(species_tree gene_trees)) {
    die "[ERROR] Input identity disagrees between manifest and compact index.\n"
        unless ($index->{input_sha256}{$kind} // '') eq ($meta->{inputs}{$kind}{sha256} // '!');
}
$gene_input //= File::Spec->rel2abs($meta->{inputs}{gene_trees}{path_as_provided}, $base);
$species_input //= File::Spec->rel2abs($meta->{inputs}{species_tree}{path_as_provided}, $base);
my $genes = freeze_verified($gene_input, File::Spec->catfile($tmp, 'genes'), $index->{input_sha256}{gene_trees});
my $species = freeze_verified($species_input, File::Spec->catfile($tmp, 'species'), $index->{input_sha256}{species_tree});
my $root = parse_species_tree_file($species);
assign_branch_ids($root);
my $axis = retained_primitive_axis($root);
verify_axis_equal($axis, $meta->{axis}{coordinates}, 'replayed reference', 'saved reference axis');
die "[ERROR] Reference axis identity differs.\n"
    unless axis_sha256($axis) eq $index->{axis_sha256} && $index->{axis_sha256} eq $dictionary->{axis_sha256};
my @taxa = sort(tip_labels($root));
my $json = JSON::PP->new->canonical(1)->utf8(1);
die "[ERROR] Taxon dictionary disagrees with original reference.\n"
    unless $json->encode(\@taxa) eq $json->encode($dictionary->{taxa});
my $ctx = prepare_reference($axis, \@taxa);
my (%by_gene, %selected, %done, %expected);
my $audit_bytes = 0;
my @axis_ids = map { $_->{branch_id} } @{$axis};
my @targets = @{$index->{fusion_catalog}};
my $view = $meta->{fusion_matrix}{view};
die "[ERROR] Invalid saved fusion view.\n" unless $view eq 'whole' || $view eq 'direct';
for my $entry (@{$index->{genes}}) {
    die "[ERROR] Duplicate gene in compact index.\n" if exists $by_gene{$entry->{gene_id}};
    $by_gene{$entry->{gene_id}} = $entry;
}
for my $gene (@wanted) {
    die "[ERROR] Unknown audit gene '$gene'.\n" unless exists $by_gene{$gene};
    $selected{$gene} = 1;
}
open my $ledger, '<:raw', $files{ledger} or die "[ERROR] Cannot read verified compact ledger: $!\n";
my $evidence_path = File::Spec->catfile($tmp, 'selected-evidence.jsonl');
open my $evidence, '>:raw', $evidence_path or die "[ERROR] Cannot stage audit evidence: $!\n";
my $record_id = 0;
visit_input_lines($genes, $index->{limits}{input_line_bytes}, sub {
    my ($line, $number, $offset, $bytes) = @_;
    return if $line =~ /^\s*\z/;
    ++$record_id;
    my $tree_pos = index($line, '(');
    die "[ERROR] Verified input has no tree at line $number.\n" if $tree_pos < 0;
    my $gene = substr($line, 0, $tree_pos);
    $gene =~ s/^\s+|\s+$//g;
    return unless $selected{$gene};
    my $entry = $by_gene{$gene};
    die "[ERROR] Saved input coordinates differ for gene '$gene'.\n"
        unless $entry->{input_record_id} == $record_id && $entry->{input_line} == $number
            && $entry->{input_offset} == $offset && $entry->{input_bytes} == $bytes;
    my $saved = read_ledger_block($ledger, $entry, $index->{limits}{ledger_block_bytes});
    my $tree = parse_newick(substr($line, $tree_pos), source => "audit gene $gene");
    my $actual = compute_gene($tree, $ctx);
    $actual->{gene_id} = $gene;
    $actual->{input_record_id} = $record_id;
    $actual->{input_line} = $number;
    # Compare the evidence that was saved, not a freshly computed replacement.
    for my $field (sort keys %{$actual}) {
        die "[ERROR] Audit mismatch for gene '$gene', field '$field'.\n"
            unless exists($saved->{$field}) && $json->encode($actual->{$field}) eq $json->encode($saved->{$field});
    }
    die "[ERROR] Unexpected compact record fields for gene '$gene'.\n"
        unless keys(%{$saved}) == keys(%{$actual});
    my $data = materialize_gene($saved, $ctx);
    for my $kind (qw(no with state truth groups cells map)) {
        $expected{$kind}{$gene} = { digest => Digest::SHA->new(256), count => 0 };
    }
    add_expected('no', $gene, $gene, @{$data->{primitive}});
    add_expected('state', $gene, $gene, @{$data->{state}});
    add_expected('truth', $gene, $gene, @{$data->{truth}});
    add_expected('map', $gene, $entry->{storage_key}, $gene, $entry->{input_line});
    for my $group (@{$data->{groups}}) {
        add_expected('groups', $gene, $gene, @{$group}{qw(split members status length)});
    }
    my @values;
    for my $target (@targets) {
        my @members = split(/\|/, $target);
        my @cell = construct_cell(\@members, $data->{owners}, $view);
        add_expected('cells', $gene, $gene, $target, @cell);
        push @values, $cell[-1];
    }
    add_expected('with', $gene, $gene, @{$data->{primitive}}, @values);
    my $encoded = $json->encode($actual);
    $audit_bytes += length($encoded);
    die "[ERROR][Capacity] Expanded audit evidence exceeds 32 MiB; select fewer genes per call.\n"
        if $audit_bytes > 32 * 1024 * 1024;
    print {$evidence} $encoded, "\n" or die "[ERROR] Cannot stage audit evidence: $!\n";
    $done{$gene} = 1;
});
close $ledger or die "[ERROR] Cannot close verified ledger: $!\n";
close $evidence or die "[ERROR] Cannot close staged audit evidence: $!\n";
die "[ERROR] An indexed audit gene is absent from verified input.\n" unless keys(%done) == keys(%selected);
my $label = $meta->{label};
my ($map_owner) = grep { $_->{relative_path} eq "$label.gene_id_map.tsv" } @{$meta->{ownership}{outputs}};
die "[ERROR] Missing saved gene-identity file fingerprint.\n" unless $map_owner;
my %tables = (
    no => [$meta->{outputs}{matrix_no_fuse}, ['gene', @axis_ids]],
    with => [$meta->{outputs}{matrix_with_fuse}, ['gene', @axis_ids, @targets]],
    state => [$meta->{outputs}{primitive_state}, ['gene', @axis_ids]],
    truth => [$meta->{projection_graph}{primitive_truth}, ['gene', @axis_ids]],
    groups => [$meta->{projection_graph}{projected_groups}, [qw(gene projected_split members status length)]],
    cells => [$meta->{outputs}{fusion_cells}, [qw(gene column reference_state recovery_status numeric_status construction_mode components deleted_members cut_components length)]],
    map => [{filename => "$label.gene_id_map.tsv", sha256 => $map_owner->{fingerprint}{sha256}}, [qw(storage_key gene_id input_line)]],
);
for my $kind (sort keys %tables) {
    my ($info, $header) = @{$tables{$kind}};
    my $name = $info->{filename};
    die "[ERROR] Invalid scientific artifact path.\n" unless defined($name)
        && !File::Spec->file_name_is_absolute($name) && $name !~ /(?:\A|\/)\.{1,2}(?:\/|\z)/;
    # Only one additional scientific table is staged at a time. Compare rows
    # from exactly the bytes whose saved fingerprint was just checked.
    my $path = freeze_verified(File::Spec->catfile($base, $name),
        File::Spec->catfile($tmp, 'scientific-table'), $info->{sha256});
    my %seen = map { $_ => {digest => Digest::SHA->new(256), count => 0} } keys %selected;
    my $rows = visit_input_lines($path, $index->{limits}{row_bytes}, sub {
        my ($line, $line_no) = @_;
        if ($line_no == 1) {
            die "[ERROR] Audit scientific header mismatch for '$name'.\n" unless $line eq join("\t", @{$header});
            return;
        }
        my @fields = split(/\t/, $line, -1);
        die "[ERROR] Malformed scientific row in '$name' at line $line_no.\n" unless @fields == @{$header};
        my $gene = $fields[$kind eq 'map' ? 1 : 0];
        return unless $selected{$gene};
        $seen{$gene}{digest}->add(encode_utf8_strict($line . "\n", 'audit table row'));
        ++$seen{$gene}{count};
    });
    die "[ERROR] Empty scientific artifact '$name'.\n" unless $rows;
    for my $gene (sort keys %selected) {
        die "[ERROR] Audit scientific row mismatch for gene '$gene' in '$name'.\n"
            unless $seen{$gene}{count} == $expected{$kind}{$gene}{count}
                && $seen{$gene}{digest}->hexdigest eq $expected{$kind}{$gene}{digest}->hexdigest;
    }
    unlink($path) or die "[ERROR] Cannot remove temporary audit table: $!\n";
}
my $report = JSON::PP->new->canonical->encode({
    status => 'PASS', scope => 'selected genes: original input replay compared with saved source ledger and scientific rows',
    scientific_tables => [sort keys %tables],
    input_sha256 => $index->{input_sha256}, gene_order => 'input-record-order',
});
$report =~ s/\}\z/,"genes":[/;
print $report or die "[ERROR] Cannot write audit report: $!\n";
my $separator = '';
visit_input_lines($evidence_path, $index->{limits}{ledger_block_bytes}, sub {
    print $separator, $_[0] or die "[ERROR] Cannot write audit evidence: $!\n";
    $separator = ',';
});
print "]}\n" or die "[ERROR] Cannot finish audit report: $!\n";

sub add_expected {
    my ($kind, $gene, @fields) = @_;
    $expected{$kind}{$gene}{digest}->add(encode_utf8_strict(join("\t", @fields) . "\n", 'audit expected row'));
    ++$expected{$kind}{$gene}{count};
}

sub freeze_verified {
    my ($source, $destination, $expected) = @_;
    die "[ERROR] Invalid expected SHA-256 for audit input.\n" unless defined($expected) && $expected =~ /\A[0-9a-f]{64}\z/;
    open my $in, '<:raw', $source or die "[ERROR] Cannot access original audit input $source: $!\n";
    open my $out, '>:raw', $destination or die "[ERROR] Cannot stage audit input: $!\n";
    my $sha = Digest::SHA->new(256);
    while (1) {
        my $n = read($in, my $chunk, 65536);
        die "[ERROR] Cannot read audit input $source: $!\n" unless defined $n;
        last unless $n;
        $sha->add($chunk);
        print {$out} $chunk or die "[ERROR] Cannot stage audit input: $!\n";
    }
    close $in or die "[ERROR] Cannot close audit input: $!\n";
    close $out or die "[ERROR] Cannot close audit snapshot: $!\n";
    die "[ERROR] Audit input hash differs: $source.\n" unless $sha->hexdigest eq $expected;
    return $destination;
}

sub usage {
    return "Usage: perl $0 --manifest LABEL.run_manifest.json --audit-gene ID [--audit-gene ID ...]\n"
        . "  [--audit-gene-file IDs.txt] [--gene-input original.genes] [--species-input original.nwk]\n";
}
