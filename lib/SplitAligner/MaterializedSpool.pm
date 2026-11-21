package SplitAligner::MaterializedSpool;

use 5.010;
use strict;
use warnings;
use Exporter qw(import);
use Digest::SHA qw(sha256_hex);
use POSIX qw(isfinite);
use SplitAligner::TextIO qw(encode_utf8_strict decode_utf8_strict);
use SplitAligner::Newick qw(is_numeric_value branch_id_num);
use SplitAligner::CompactKernel qw(compare_patterns);

our @EXPORT_OK = qw(write_spool_block read_spool_block);
my $SCHEMA = 'SplitAligner-materialized-spool-v1';

# Private execution cache, never a second reference-truth authority. A block is:
# H<TAB>schema<TAB>record_id<TAB>input_line<TAB>UTF8_gene_hex<TAB>reference_id<TAB>JSON_block_sha
# N/S/T<TAB>gene<TAB>primitive/state/truth values in reference-axis order
# G<TAB>gene<TAB>split<TAB>members<TAB>status<TAB>length (split-text order)
# C<TAB>members<TAB>recovered<TAB>length (numeric member-pattern order)
# All lines end in LF. Entries contain offset, bytes, sha256, input_record_id,
# input_line, gene_id, reference_id, ledger_sha256, groups and components.
# The in-memory write-time entry SHA is authority; a read never re-seals bytes.
# read_spool_block also returns raw_rows->{no,state,truth,groups}, removing only
# the two-byte tag prefix. Each raw row includes gene and its final LF.

sub _fail { die "[ERROR][Materialized spool] $_[0]\n"; }
sub _plain {
    my ($value) = @_;
    return defined($value) && !ref($value) && $value !~ /[\t\n]/;
}
sub _positive { return defined($_[0]) && !ref($_[0]) && $_[0] =~ /\A[1-9][0-9]*\z/; }
sub _sha { return defined($_[0]) && !ref($_[0]) && $_[0] =~ /\A[0-9a-f]{64}\z/; }

sub _limits {
    my ($limits) = @_;
    _fail('Invalid spool or scientific-row capacity limit.')
        unless ref($limits) eq 'HASH' && _positive($limits->{spool_block_bytes})
            && _positive($limits->{row_bytes});
}

sub _capacity {
    my ($kind, $value, $limits) = @_;
    die "[ERROR][Capacity] $kind=$value exceeds limit $limits->{$kind}; no results published.\n"
        if $value > $limits->{$kind};
}

sub _length {
    my ($value) = @_;
    _fail('Invalid length token.') unless _plain($value);
    return if $value eq 'NA';
    # Check a copy so numeric coercion cannot change a record's raw token flags.
    my $copy = $value;
    _fail('Invalid length token.') unless is_numeric_value($copy);
    my $number = unpack('d', pack('d', $copy));
    my ($significand) = split(/[eE]/, $copy);
    _fail('Length is outside binary64 range.')
        unless isfinite($number) && !($number == 0 && $significand =~ /[1-9]/);
}

sub _context {
    my ($ctx) = @_;
    _fail('Invalid reference context.') unless ref($ctx) eq 'HASH'
        && ref($ctx->{axis_ids}) eq 'ARRAY' && @{$ctx->{axis_ids}}
        && ref($ctx->{axis_index}) eq 'HASH' && _sha($ctx->{reference_id});
}

sub _members {
    my ($id, $ctx) = @_;
    _fail('Invalid component member pattern.') unless _plain($id)
        && $id =~ /\AB[0-9]+(?:\|B[0-9]+)*\z/;
    my @members = split(/\|/, $id);
    my $last;
    for my $member (@members) {
        my $number = branch_id_num($member);
        _fail('Component member is absent from the axis or is not strictly ordered.')
            unless exists($ctx->{axis_index}{$member}) && (!defined($last) || $number > $last);
        $last = $number;
    }
    return \@members;
}

# Linear cross-field validation of the already computed component partition.
# No projection, split matching, JSON decoding or branch-length reconstruction.
sub _validate_data {
    my ($data, $ctx) = @_;
    _context($ctx);
    _fail('Invalid materialized data.') unless ref($data) eq 'HASH'
        && ref($data->{owners}) eq 'HASH' && ref($data->{groups}) eq 'ARRAY';
    for my $field (qw(primitive state truth)) {
        _fail("Invalid $field row width.") unless ref($data->{$field}) eq 'ARRAY'
            && @{$data->{$field}} == @{$ctx->{axis_ids}};
        _fail("Invalid $field field.") if grep { !_plain($_) } @{$data->{$field}};
    }
    my %components;
    for my $member (keys %{$data->{owners}}) {
        my $owner = $data->{owners}{$member};
        _fail('Invalid component owner.') unless exists($ctx->{axis_index}{$member})
            && ref($owner) eq 'HASH' && _plain($owner->{id});
        if (exists $components{$owner->{id}}) {
            my $previous = $components{$owner->{id}};
            for my $key (qw(size recovered length)) {
                _fail('Inconsistent aliases for the same component.')
                    unless defined($owner->{$key}) && defined($previous->{$key})
                        && !ref($owner->{$key}) && $owner->{$key} eq $previous->{$key};
            }
        } else {
            $components{$owner->{id}} = $owner;
        }
    }
    _fail('A gene must have at least one complete surviving component.') unless keys %components;
    my (%covered, %fused);
    my @ordered;
    for my $id (sort { compare_patterns($a, $b) } keys %components) {
        my $owner = $components{$id};
        my $members = _members($id, $ctx);
        _fail('Invalid component size or recovery flag.')
            unless _positive($owner->{size}) && $owner->{size} == @{$members}
                && defined($owner->{recovered}) && !ref($owner->{recovered})
                && $owner->{recovered} =~ /\A[01]\z/;
        _length($owner->{length});
        _fail('An unrecovered component must have length NA.')
            if !$owner->{recovered} && $owner->{length} ne 'NA';
        my $truth = @{$members} > 1 ? 'NA_fuse' : $owner->{recovered} ? 'Mapped' : 'not_displayed';
        my $state = $owner->{recovered} ? (@{$members} > 1 ? 'F' : 'D') : 'U';
        my $primitive = @{$members} == 1 && $owner->{recovered} ? $owner->{length} : 'NA';
        for my $member (@{$members}) {
            _fail('Repeated, missing or inconsistent component ownership.')
                if $covered{$member}++ || !exists($data->{owners}{$member})
                    || $data->{owners}{$member}{id} ne $id;
            my $i = $ctx->{axis_index}{$member};
            _fail('Component evidence disagrees with primitive/state/truth rows.')
                unless $data->{truth}[$i] eq $truth && $data->{state}[$i] eq $state
                    && $data->{primitive}[$i] eq $primitive;
        }
        $fused{$id} = $owner if @{$members} > 1;
        push @ordered, $owner;
    }
    _fail('Component owner map contains extraneous members.')
        unless keys(%covered) == keys(%{$data->{owners}});
    for my $i (0 .. $#{$ctx->{axis_ids}}) {
        next if $covered{$ctx->{axis_ids}[$i]};
        _fail('Components do not exactly cover the nonstructural axis.')
            unless $data->{truth}[$i] eq 'NA_struct' && $data->{state}[$i] eq 'S'
                && $data->{primitive}[$i] eq 'NA';
    }
    my (%groups_seen, $last_split);
    for my $group (@{$data->{groups}}) {
        _fail('Invalid or incomplete group row.') unless ref($group) eq 'HASH'
            && !grep { !_plain($group->{$_}) } qw(split members status length);
        _fail('Group splits are empty, repeated or unordered.')
            unless length($group->{split}) && (!defined($last_split) || $last_split lt $group->{split});
        $last_split = $group->{split};
        my $owner = $fused{$group->{members}};
        _fail('Group rows do not correspond one-to-one with complete fused components.')
            unless $owner && !$groups_seen{$group->{members}}++;
        my $status = $owner->{recovered} ? 'recovered' : 'not_displayed';
        _fail('Group recovery or length disagrees with its component.')
            unless $group->{status} eq $status && $group->{length} eq $owner->{length};
    }
    _fail('A fused component is missing its group row.') unless keys(%groups_seen) == keys(%fused);
    return \@ordered;
}

sub _scientific_row {
    my ($limits, @fields) = @_;
    _fail('Invalid field in scientific row.') if grep { !_plain($_) } @fields;
    my $raw = encode_utf8_strict(join("\t", @fields) . "\n", 'materialized scientific row');
    _capacity('row_bytes', length($raw), $limits);
    return $raw;
}

sub write_spool_block {
    my ($fh, $record, $data, $ctx, $ledger_sha256, $limits) = @_;
    _limits($limits);
    _context($ctx);
    _fail('Invalid source-record binding.') unless ref($record) eq 'HASH'
        && _positive($record->{input_record_id}) && _positive($record->{input_line})
        && _plain($record->{gene_id}) && length($record->{gene_id})
        && ($record->{reference_id} // '') eq $ctx->{reference_id} && _sha($ledger_sha256);
    my $owners = _validate_data($data, $ctx);
    my $gene = $record->{gene_id};
    my $gene_hex = unpack('H*', encode_utf8_strict($gene, 'materialized gene identity'));
    my $raw = encode_utf8_strict(join("\t", 'H', $SCHEMA, $record->{input_record_id}, $record->{input_line},
        $gene_hex, $ctx->{reference_id}, $ledger_sha256) . "\n", 'spool header');
    for my $pair (['N', 'primitive'], ['S', 'state'], ['T', 'truth']) {
        $raw .= $pair->[0] . "\t" . _scientific_row($limits, $gene, @{$data->{$pair->[1]}});
        _capacity('spool_block_bytes', length($raw), $limits);
    }
    for my $group (@{$data->{groups}}) {
        $raw .= "G\t" . _scientific_row($limits, $gene, @{$group}{qw(split members status length)});
        _capacity('spool_block_bytes', length($raw), $limits);
    }
    for my $owner (@{$owners}) {
        $raw .= encode_utf8_strict('C' . "\t" . join("\t", @{$owner}{qw(id recovered length)}) . "\n",
            'spool component row');
        _capacity('spool_block_bytes', length($raw), $limits);
    }
    # Decoded ASCII B aliases still carry Perl's UTF8 flag. Appending them to
    # already encoded rows upgrades the entire growing buffer; subsequent
    # length() calls can rescan it. Check before SHA's byte coercion hides it.
    _fail('Internal spool buffer must remain an encoded byte string.') if utf8::is_utf8($raw);
    my $offset = tell($fh);
    _fail("Cannot locate spool write position: $!.") unless defined($offset) && $offset >= 0;
    my $entry = { offset => $offset, bytes => length($raw), sha256 => sha256_hex($raw),
        input_record_id => $record->{input_record_id}, input_line => $record->{input_line},
        gene_id => $gene, reference_id => $ctx->{reference_id}, ledger_sha256 => $ledger_sha256,
        groups => scalar(@{$data->{groups}}), components => scalar(@{$owners}) };
    print {$fh} $raw or _fail("Cannot write spool block: $!.");
    return $entry;
}

sub read_spool_block {
    my ($fh, $entry, $ctx, $limits) = @_;
    _limits($limits);
    _context($ctx);
    _fail('Invalid block boundary or source binding.') unless ref($entry) eq 'HASH'
        && _positive($entry->{bytes}) && defined($entry->{offset}) && !ref($entry->{offset})
        && $entry->{offset} =~ /\A[0-9]+\z/ && _sha($entry->{sha256}) && _sha($entry->{ledger_sha256})
        && _positive($entry->{input_record_id}) && _positive($entry->{input_line})
        && _plain($entry->{gene_id}) && length($entry->{gene_id})
        && ($entry->{reference_id} // '') eq $ctx->{reference_id}
        && defined($entry->{groups}) && !ref($entry->{groups}) && $entry->{groups} =~ /\A[0-9]+\z/
        && _positive($entry->{components});
    _capacity('spool_block_bytes', $entry->{bytes}, $limits);
    seek($fh, $entry->{offset}, 0) or _fail("Cannot seek spool block: $!.");
    my $raw = '';
    while (length($raw) < $entry->{bytes}) {
        my $n = read($fh, my $part, $entry->{bytes} - length($raw));
        _fail('Short or failed spool read.') unless defined($n) && $n > 0;
        $raw .= $part;
    }
    _fail('Block differs from its write-time fingerprint.') unless sha256_hex($raw) eq $entry->{sha256};
    _fail('Block lacks its final newline.') unless $raw =~ /\n\z/;
    my @lines = split(/\n/, $raw, -1);
    pop @lines;
    _fail('Incomplete spool block.') unless @lines >= 5;
    my @header = split(/\t/, shift(@lines), -1);
    my $gene_hex = unpack('H*', encode_utf8_strict($entry->{gene_id}, 'spool gene binding'));
    my @expected = ('H', $SCHEMA, @{$entry}{qw(input_record_id input_line)}, $gene_hex,
        $entry->{reference_id}, $entry->{ledger_sha256});
    _fail('Header differs from its source-record binding.') unless @header == @expected
        && join("\t", @header) eq join("\t", @expected);
    my $data = { groups => [], owners => {}, raw_rows => { groups => [] } };
    for my $pair (['N', 'primitive', 'no'], ['S', 'state', 'state'], ['T', 'truth', 'truth']) {
        my $line = shift @lines;
        _fail('Scientific row is missing or out of order.') unless substr($line, 0, 2) eq $pair->[0] . "\t";
        my $row = substr($line, 2) . "\n";
        _capacity('row_bytes', length($row), $limits);
        my @fields = split(/\t/, decode_utf8_strict(substr($line, 2), 'spool scientific row'), -1);
        _fail('Scientific row gene or width differs.') unless @fields == @{$ctx->{axis_ids}} + 1
            && shift(@fields) eq $entry->{gene_id};
        $data->{$pair->[1]} = \@fields;
        $data->{raw_rows}{$pair->[2]} = $row;
    }
    my ($components, $saw_component, $previous) = (0, 0, undef);
    for my $line (@lines) {
        my @fields = split(/\t/, decode_utf8_strict($line, 'spool component or group row'), -1);
        my $tag = shift @fields;
        if ($tag eq 'G') {
            _fail('Group row is out of order or malformed.') if $saw_component || @fields != 5;
            _fail('Group row belongs to a different gene.') unless shift(@fields) eq $entry->{gene_id};
            my $row = substr($line, 2) . "\n";
            _capacity('row_bytes', length($row), $limits);
            my %group;
            @group{qw(split members status length)} = @fields;
            push @{$data->{groups}}, \%group;
            push @{$data->{raw_rows}{groups}}, $row;
        } elsif ($tag eq 'C') {
            $saw_component = 1;
            _fail('Malformed component row.') unless @fields == 3;
            my ($id, $recovered, $length) = @fields;
            _fail('Invalid component recovery flag.') unless $recovered =~ /\A[01]\z/;
            my $members = _members($id, $ctx);
            _fail('Components are repeated or out of order.')
                if defined($previous) && compare_patterns($previous, $id) >= 0;
            $previous = $id;
            my $owner = { id => $id, size => scalar(@{$members}), recovered => 0 + $recovered, length => $length };
            for my $member (@{$members}) {
                _fail('Components overlap.') if exists $data->{owners}{$member};
                $data->{owners}{$member} = $owner;
            }
            ++$components;
        } else {
            _fail('Unknown or misplaced spool row type.');
        }
    }
    _fail('Block component/group counts differ from the index.')
        unless $components == $entry->{components} && @{$data->{groups}} == $entry->{groups};
    _validate_data($data, $ctx);
    return $data;
}

1;
