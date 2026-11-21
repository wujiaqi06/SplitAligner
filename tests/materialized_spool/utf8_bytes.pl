#!/usr/bin/env perl
use strict;
use warnings;
use utf8;
use FindBin qw($RealBin);
use lib "$RealBin/../../lib";
use File::Temp qw(tempfile);
use Digest::SHA qw(sha256_hex);
use Test::More;
use JSON::PP;
use SplitAligner::TextIO qw(encode_utf8_strict);
use SplitAligner::Newick qw(parse_newick assign_branch_ids retained_primitive_axis tip_labels);
use SplitAligner::CompactKernel qw(prepare_reference compute_gene materialize_gene compare_patterns);
use SplitAligner::MaterializedSpool qw(write_spool_block read_spool_block);

my $root = parse_newick('((A,(B,C)),(D,E));');
assign_branch_ids($root);
my $axis = retained_primitive_axis($root);
utf8::upgrade($_->{branch_id}) for @{$axis};
my $ctx = prepare_reference($axis, [tip_labels($root)]);
my $record = compute_gene(parse_newick('((A:1.00,(B:2,C:3):4):5,(D:-0,E:7):8);'), $ctx);
@{$record}{qw(gene_id input_record_id input_line)} = ('猫样本', 1, 7);
my $data = materialize_gene($record, $ctx);
ok(scalar(grep { utf8::is_utf8($_->{id}) } values %{$data->{owners}}),
    'fixture component IDs carry UTF8 flags although IDs are ASCII');
my $before_json = JSON::PP->new->canonical->utf8->encode($record);
my $limits = { spool_block_bytes => 1048576, row_bytes => 1048576 };
my ($fh, $path) = tempfile(UNLINK => 1);
binmode $fh, ':raw';
my $entry = write_spool_block($fh, $record, $data, $ctx, 'a' x 64, $limits);
is(JSON::PP->new->canonical->utf8->encode($record), $before_json,
    'writing byte rows preserves JSON record tokens and scalar representation');
seek($fh, 0, 0) or die $!;
my $written = do { local $/; <$fh> };
my $text = join("\t", 'H', 'SplitAligner-materialized-spool-v1', 1, 7,
    unpack('H*', encode_utf8_strict($record->{gene_id})), $ctx->{reference_id}, 'a' x 64) . "\n";
for my $pair (['N', 'primitive'], ['S', 'state'], ['T', 'truth']) {
    $text .= join("\t", $pair->[0], $record->{gene_id}, @{$data->{$pair->[1]}}) . "\n";
}
for my $group (@{$data->{groups}}) {
    $text .= join("\t", 'G', $record->{gene_id}, @{$group}{qw(split members status length)}) . "\n";
}
my %owners = map { $_->{id} => $_ } values %{$data->{owners}};
for my $id (sort { compare_patterns($a, $b) } keys %owners) {
    $text .= join("\t", 'C', @{$owners{$id}}{qw(id recovered length)}) . "\n";
}
my $expected = encode_utf8_strict($text);
is($written, $expected, 'byte buffer equals one canonical UTF8 encoding of the complete TSV');
is($entry->{bytes}, length($expected), 'entry length counts bytes rather than characters');
is($entry->{sha256}, sha256_hex($expected), 'write-time fingerprint binds exact UTF8 bytes');
my $decoded = read_spool_block($fh, $entry, $ctx, $limits);
is_deeply($decoded->{primitive}, $data->{primitive}, 'raw numeric tokens survive flagged aliases');
is_deeply($decoded->{owners}, $data->{owners}, 'component evidence survives flagged aliases');

# Local helper injection emulates a missing encode call. This regression checks
# the guard before SHA can silently coerce the upgraded buffer back to bytes.
my $encoder = \&SplitAligner::MaterializedSpool::encode_utf8_strict;
for my $target ('spool header', 'spool component row') {
    my ($bad_fh, $bad_path) = tempfile(UNLINK => 1);
    my $success;
    {
        no warnings 'redefine';
        local *SplitAligner::MaterializedSpool::encode_utf8_strict = sub {
            my ($value, $context) = @_;
            if (($context // '') eq $target) {
                my $bad = $value;
                utf8::upgrade($bad);
                return $bad;
            }
            return $encoder->(@_);
        };
        $success = eval { write_spool_block($bad_fh, $record, $data, $ctx, 'a' x 64, $limits); 1 };
    }
    ok(!$success && $@ =~ /Internal spool buffer must remain an encoded byte string/,
        "$target flag regression fails before SHA can hide promotion");
    is(-s $bad_path, 0, 'no flagged block is written after invariant failure');
}
done_testing();
