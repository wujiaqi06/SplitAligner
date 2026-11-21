package SplitAligner::CompactPipeline;

use strict;
use warnings;
use Exporter qw(import);
use Digest::SHA qw(sha256_hex);
use JSON::PP;
use File::Path qw(make_path);
use File::Temp;
use Time::HiRes qw(time);
use SplitAligner::Newick qw(parse_newick tip_labels);
use SplitAligner::Provenance qw(read_axis_ledger axis_sha256 file_sha256 write_json_file);
use SplitAligner::TextIO qw(decode_utf8_strict encode_utf8_strict open_utf8_writer print_utf8 close_utf8_writer);
use SplitAligner::CompactKernel qw(prepare_reference compute_gene_materialized construct_cell compare_patterns);
use SplitAligner::MaterializedSpool qw(write_spool_block read_spool_block);

our @EXPORT_OK = qw(run_compact_matrix visit_input_lines read_ledger_block compact_limits);
my $JSON = JSON::PP->new->canonical(1)->utf8(1);

# Capacity limits bound the objects retained by this implementation; encoded
# payload budgets are not claims about the allocator's total resident memory.
sub compact_limits {
    return {
        genes => 100_000, index_bytes => 32 * 1024 * 1024,
        fusion_patterns => 100_000, catalog_bytes => 32 * 1024 * 1024,
        fusion_members => 1_000_000, columns => 200_000,
        input_line_bytes => 16 * 1024 * 1024,
        ledger_block_bytes => 64 * 1024 * 1024,
        spool_block_bytes => 64 * 1024 * 1024,
        row_bytes => 16 * 1024 * 1024,
        reference_bitset_bytes => 64 * 1024 * 1024,
    };
}

sub _capacity {
    my ($limits, $kind, $actual) = @_;
    die "[ERROR][Capacity] $kind=$actual exceeds limit $limits->{$kind}; no results published.\n"
        if $actual > $limits->{$kind};
}

# Raw byte accounting is independent of Perl's Unicode character positions.
# Reads are bounded, including an unterminated oversized input line.
sub visit_input_lines {
    my ($path, $maximum, $visit, $source) = @_;
    $source //= $path;
    open my $fh, '<:raw', $path or die "[ERROR] Cannot open $path: $!\n";
    my ($buffer, $offset, $line_no) = ('', 0, 0);
    my $emit = sub {
        my ($raw) = @_;
        my $bytes = length($raw);
        die "[ERROR][Capacity] input_line_bytes=$bytes exceeds limit $maximum.\n" if $bytes > $maximum;
        ++$line_no;
        die "[ERROR][$source] UTF-8 BOM is not permitted.\n"
            if $line_no == 1 && substr($raw, 0, 3) eq "\xEF\xBB\xBF";
        $raw =~ s/\n\z//;
        $raw =~ s/\r\z//;
        my $text = decode_utf8_strict($raw, "$source line $line_no");
        $visit->($text, $line_no, $offset, $bytes);
        $offset += $bytes;
    };
    while (1) {
        my $read = read($fh, my $chunk, 65536);
        die "[ERROR] Cannot read $path: $!\n" unless defined $read;
        last unless $read;
        $buffer .= $chunk;
        while ((my $end = index($buffer, "\n")) >= 0) {
            $emit->(substr($buffer, 0, $end + 1, ''));
        }
        die "[ERROR][Capacity] input_line_bytes exceeds limit $maximum.\n" if length($buffer) > $maximum;
    }
    $emit->($buffer) if length($buffer);
    close $fh or die "[ERROR] Cannot close $path: $!\n";
    return $line_no;
}

sub read_ledger_block {
    my ($fh, $entry, $maximum) = @_;
    my $length = $entry->{ledger_bytes};
    die "[ERROR] Invalid compact ledger block boundary.\n"
        unless defined($length) && $length =~ /\A[1-9]\d*\z/ && $length <= $maximum
        && defined($entry->{ledger_offset}) && $entry->{ledger_offset} =~ /\A\d+\z/;
    seek($fh, $entry->{ledger_offset}, 0) or die "[ERROR] Cannot seek compact ledger: $!\n";
    my $raw = '';
    while (length($raw) < $length) {
        my $n = read($fh, my $part, $length - length($raw));
        die "[ERROR] Short or failed read of compact ledger for gene '$entry->{gene_id}'.\n" unless defined($n) && $n > 0;
        $raw .= $part;
    }
    die "[ERROR] Compact ledger block fingerprint differs for gene '$entry->{gene_id}'.\n"
        unless sha256_hex($raw) eq ($entry->{sha256} // '');
    die "[ERROR] Compact ledger block is not one newline-terminated JSON record.\n"
        unless $raw =~ /\n\z/ && $raw !~ /\n./s;
    my $record = eval { $JSON->decode($raw) };
    die "[ERROR] Invalid compact ledger JSON: $@\n" if $@;
    die "[ERROR] Compact ledger identity or record counts differ for gene '$entry->{gene_id}'.\n"
        unless ref($record) eq 'HASH'
        && ($record->{gene_id} // '') eq $entry->{gene_id}
        && ($record->{input_record_id} // 0) == $entry->{input_record_id}
        && ($record->{input_line} // 0) == $entry->{input_line}
        && ref($record->{observations}) eq 'ARRAY' && @{$record->{observations}} == $entry->{observations}
        && ref($record->{components}) eq 'ARRAY' && @{$record->{components}} == $entry->{components};
    return $record;
}

sub _row {
    my ($fh, $limits, @fields) = @_;
    my $bytes = encode_utf8_strict(join("\t", @fields) . "\n", 'scientific output row');
    _capacity($limits, 'row_bytes', length($bytes));
    print {$fh} $bytes or die "[ERROR] Cannot write scientific output: $!\n";
}

sub run_compact_matrix {
    my %arg = @_;
    my $label = $arg{label};
    # Read only the frozen snapshot; use the caller's path spelling solely in
    # diagnostics. Direct API callers may omit this optional display label.
    my $gene_source = $arg{gene_source} // $arg{gene_path};
    my $limits = compact_limits();
    if ($arg{limits}) {
        for my $key (keys %{$arg{limits}}) {
            die "[ERROR] Invalid capacity limit '$key'.\n"
                unless exists($limits->{$key}) && $arg{limits}{$key} =~ /\A[1-9]\d*\z/;
            $limits->{$key} = $arg{limits}{$key};
        }
    }
    my $axis = read_axis_ledger($arg{axis_file});
    my @taxa = sort(tip_labels($arg{species_root}));
    _capacity($limits, 'reference_bitset_bytes', @{$axis} * int((@taxa + 7) / 8));
    _capacity($limits, 'columns', scalar(@{$axis}));
    my $ctx = prepare_reference($axis, \@taxa);
    my @axis_ids = map { $_->{branch_id} } @{$axis};
    my $ledger_path = "$label.compact.jsonl";
    my $index_path = "$label.compact.index.json";
    my $dictionary_path = "$label.taxon_dictionary.json";
    my $dictionary_sha256 = write_json_file($dictionary_path, {
        schema_version => 'SplitAligner-taxon-dictionary-v1', taxa => \@taxa,
        axis_sha256 => axis_sha256($axis), bit_order => 'vec-LSB-first',
        bitset_bytes => int((@taxa + 7) / 8),
    });
    open my $ledger, '>:raw', $ledger_path or die "[ERROR] Cannot create $ledger_path: $!\n";
    my $ledger_digest = Digest::SHA->new(256);
    my $spool = File::Temp->new(TEMPLATE => '.materialized-XXXXXX', DIR => '.', UNLINK => 1);
    binmode $spool, ':raw' or die "[ERROR] Cannot set spool byte mode: $!\n";
    my (@entries, %seen, %catalog, %patterns, %spool_entries);
    my ($index_bytes, $catalog_bytes, $member_count) = (0, 0, 0);
    my %timing = map { $_ => 0 } qw(parse_compute ledger_write spool_write spool_read ledger_verify fusion_construct fusion_write scientific_write);
    visit_input_lines($arg{gene_path}, $limits->{input_line_bytes}, sub {
        my ($line, $line_no, $input_offset, $input_bytes) = @_;
        return if $line =~ /^\s*\z/;
        my $tree_pos = index($line, '(');
        die "[ERROR][$gene_source line $line_no] No '(' found to start a gene tree.\n" if $tree_pos < 0;
        my $gene = substr($line, 0, $tree_pos);
        $gene =~ s/^\s+|\s+$//g;
        die "[ERROR][$gene_source line $line_no] Empty gene ID.\n" unless length($gene);
        die "[ERROR][$gene_source line $line_no gene $gene] Duplicate gene ID; first seen on line $seen{$gene}.\n" if exists $seen{$gene};
        die "[ERROR][$gene_source line $line_no] Gene ID contains a tab and cannot be represented safely: '$gene'.\n" if $gene =~ /\t/;
        _capacity($limits, 'genes', scalar(@entries) + 1);
        my $started = time;
        my ($record, $data);
        my $source = "$gene_source line $line_no gene $gene";
        # Parser errors already contain this context (and possibly a character
        # offset). Do not add a second ERROR/context wrapper to them.
        my $root = parse_newick(substr($line, $tree_pos), source => $source);
        my $ok = eval {
            ($record, $data) = compute_gene_materialized($root, $ctx);
            1;
        };
        unless ($ok) {
            my $error = $@;
            # Preserve the kernel's category and message; add only missing gene
            # context, with one leading ERROR marker.
            $error =~ s/\A\[ERROR\]//;
            die "[ERROR][$source] $error";
        }
        $record->{gene_id} = $gene;
        $record->{input_record_id} = scalar(@entries) + 1;
        $record->{input_line} = $line_no;
        for my $group (@{$data->{groups}}) {
            my $p = $group->{members};
            next if exists $catalog{$p};
            my @members = split(/\|/, $p);
            _capacity($limits, 'fusion_patterns', scalar(keys %catalog) + 1);
            _capacity($limits, 'catalog_bytes', $catalog_bytes += length($p) + 1);
            _capacity($limits, 'fusion_members', $member_count += @members);
            _capacity($limits, 'columns', @axis_ids + scalar(keys %catalog) + 1);
            $catalog{$p} = \@members;
        }
        $patterns{$_} = 1 for @{$data->{patterns}};
        $timing{parse_compute} += time - $started;
        $started = time;
        my $raw = $JSON->encode($record) . "\n";
        _capacity($limits, 'ledger_block_bytes', length($raw));
        my $position = tell($ledger);
        die "[ERROR] Cannot locate compact ledger position: $!\n" if $position < 0;
        my $entry = {
            gene_id => $gene, input_record_id => $record->{input_record_id}, input_line => $line_no,
            input_offset => $input_offset, input_bytes => $input_bytes,
            ledger_offset => $position, ledger_bytes => length($raw), sha256 => sha256_hex($raw),
            observations => scalar(@{$record->{observations}}), components => scalar(@{$record->{components}}),
        };
        _capacity($limits, 'index_bytes', $index_bytes += length($JSON->encode($entry)) + 1);
        print {$ledger} $raw or die "[ERROR] Cannot write compact ledger: $!\n";
        $ledger_digest->add($raw);
        push @entries, $entry;
        $seen{$gene} = $line_no;
        $timing{ledger_write} += time - $started;
        $started = time;
        my $spooled = write_spool_block($spool, $record, $data, $ctx, $entry->{sha256}, $limits);
        $spool_entries{$gene} = $spooled;
        _capacity($limits, 'index_bytes', $index_bytes += length($JSON->encode($spooled)) + 1);
        $timing{spool_write} += time - $started;
    }, $gene_source);
    close $ledger or die "[ERROR] Cannot close compact ledger: $!\n";
    $spool->flush or die "[ERROR] Cannot flush materialized spool: $!\n";
    my $ledger_sha256 = $ledger_digest->hexdigest;
    die "[ERROR][$gene_source] No non-empty gene-tree records found.\n" unless @entries;
    my @ordered = sort { $a->{gene_id} cmp $b->{gene_id} } @entries;
    my @targets = sort { compare_patterns($a, $b) } keys %catalog;
    my $ordinal = 0;
    $_->{storage_key} = sprintf('g%08d', ++$ordinal) for @ordered;
    my $index = {
        schema_version => 'SplitAligner-compact-ledger-v1',
        ledger_file => $ledger_path, ledger_sha256 => $ledger_sha256,
        dictionary_file => $dictionary_path, dictionary_sha256 => $dictionary_sha256,
        axis_sha256 => axis_sha256($axis), parser_sha256 => $arg{parser_sha256},
        source_coordinates => 'zero-based Newick child preorder; root=0; observation_id is 1-based hex-key order within gene',
        input_sha256 => $arg{input_sha256}, input_bytes => $arg{input_bytes},
        ledger_order => 'input-record-order', gene_order => 'Perl-character-cmp',
        limits => $limits, tg_cache => 'disabled', fusion_catalog => \@targets, genes => \@ordered,
        assembly => 'private-materialized-tsv-v1; no JSON decoding or reference projection in phase two',
    };
    my $index_sha256 = write_json_file($index_path, $index);
    _capacity($limits, 'index_bytes', -s $index_path);
    my $graph_dir = "${label}_split_branch_label";
    make_path($graph_dir);
    my %paths = (
        no => "$label.matrix_no_fuse.txt", with => "$label.matrix_with_fuse.txt",
        state => "$label.primitive_state.tsv", map => "$label.gene_id_map.tsv",
        truth => "$graph_dir/$label.primitive_truth.tsv", groups => "$graph_dir/$label.projected_groups.tsv",
        patterns => "$graph_dir/$label.branch_patterns.txt", cells => "$label.fusion_cells.tsv",
    );
    my %fh = map { $_ => open_utf8_writer($paths{$_}) } keys %paths;
    _row($fh{$_}, $limits, 'gene', @axis_ids) for qw(no state truth);
    _row($fh{with}, $limits, 'gene', @axis_ids, @targets);
    _row($fh{map}, $limits, qw(storage_key gene_id input_line));
    _row($fh{groups}, $limits, qw(gene projected_split members status length));
    _row($fh{cells}, $limits, qw(gene column reference_state recovery_status numeric_status construction_mode components deleted_members cut_components length));
    _row($fh{patterns}, $limits, $_) for sort keys %patterns;
    my %truth_counts = map { $_ => 0 } qw(Mapped NA_struct NA_fuse not_displayed);
    my %group_counts = (recovered => 0, not_displayed => 0);
    for my $entry (@ordered) {
        my $started = time;
        my $data = read_spool_block($spool, $spool_entries{$entry->{gene_id}}, $ctx, $limits);
        $timing{spool_read} += time - $started;
        $started = time;
        my $gene = $entry->{gene_id};
        _row($fh{map}, $limits, $entry->{storage_key}, $gene, $entry->{input_line});
        for my $kind (qw(no state truth)) {
            print {$fh{$kind}} $data->{raw_rows}{$kind}
                or die "[ERROR] Cannot copy materialized scientific row: $!\n";
        }
        $truth_counts{$_}++ for @{$data->{truth}};
        for my $line (@{$data->{raw_rows}{groups}}) {
            print {$fh{groups}} $line or die "[ERROR] Cannot copy materialized group row: $!\n";
        }
        for my $g (@{$data->{groups}}) {
            $group_counts{$g->{status}}++;
        }
        $timing{scientific_write} += time - $started;
        my @values;
        for my $target (@targets) {
            $started = time;
            my @cell = construct_cell($catalog{$target}, $data->{owners}, $arg{view});
            $timing{fusion_construct} += time - $started;
            $started = time;
            _row($fh{cells}, $limits, $gene, $target, @cell);
            push @values, $cell[-1];
            $timing{fusion_write} += time - $started;
        }
        $started = time;
        _row($fh{with}, $limits, $gene, @{$data->{primitive}}, @values);
        $timing{scientific_write} += time - $started;
    }
    close $spool or die "[ERROR] Cannot close materialized spool: $!\n";
    # Written bytes, not a later reread, establish the ledger's authority.
    # Reading it here is only a streaming fingerprint check, never JSON decode.
    my $verify_started = time;
    die "[ERROR] Written compact ledger differs from its accumulated fingerprint.\n"
        unless file_sha256($ledger_path) eq $ledger_sha256;
    $timing{ledger_verify} += time - $verify_started;
    close_utf8_writer($fh{$_}, $paths{$_}) for keys %fh;
    print STDERR '[INFO] Compact matrix stage seconds (disjoint instrumented regions): ',
        JSON::PP->new->canonical->encode(\%timing), "\n";
    return {
        gene_identity => [ map { +{ storage_key => $_->{storage_key}, gene_id => $_->{gene_id}, input_line => $_->{input_line} } } @ordered ],
        gene_count => scalar(@ordered), cell_count => scalar(@ordered) * @axis_ids,
        primitive_counts => \%truth_counts, fused_group_counts => \%group_counts,
        ledger => $ledger_path, index => $index_path, dictionary => $dictionary_path,
        limits => $limits, fusion_pattern_count => scalar(@targets),
        artifact_sha256 => { ledger => $ledger_sha256, index => $index_sha256, dictionary => $dictionary_sha256 },
    };
}

1;
