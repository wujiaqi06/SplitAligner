#!/usr/bin/env perl
use 5.010;
use strict;
use warnings;
use utf8;
use FindBin qw($RealBin);
use lib "$RealBin/../../lib";
use Test::More;
use File::Temp qw(tempfile);
use Digest::SHA qw(sha256_hex);
use JSON::PP;
use SplitAligner::Newick qw(parse_newick assign_branch_ids retained_primitive_axis tip_labels);
use SplitAligner::TextIO qw(encode_utf8_strict);
use SplitAligner::CompactKernel qw(prepare_reference compute_gene materialize_gene construct_cell);
use SplitAligner::CompactPipeline qw(compact_limits);
use SplitAligner::MaterializedSpool qw(write_spool_block read_spool_block);

# The source ledger and materializer are the expected-value oracle here. The
# execution spool is only a sealed cache. Artificially rebinding its SHA tests
# syntax/cross-field checks, not an impossible proof that every internally
# consistent replacement still agrees with a JSON record the reader never reads.
my $root = parse_newick('(((A,B),C),(D,E));');
assign_branch_ids($root);
my $ctx = prepare_reference(retained_primitive_axis($root), [tip_labels($root)]);
my $limits = compact_limits();
$limits->{spool_block_bytes} = 64 * 1024 * 1024;
my $json = JSON::PP->new->canonical(1)->utf8(1);

sub fixture {
    my ($tree, $gene, $number) = @_;
    my $record = compute_gene(parse_newick($tree), $ctx);
    $record->{gene_id} = $gene;
    $record->{input_record_id} = $number;
    $record->{input_line} = $number + 3;
    my $data = materialize_gene($record, $ctx);
    return ($record, $data, sha256_hex($json->encode($record) . "\n"));
}

sub file_with {
    my ($bytes) = @_;
    my ($fh, $path) = tempfile(UNLINK => 1);
    binmode($fh, ':raw') or die $!;
    print {$fh} $bytes or die $!;
    seek($fh, 0, 0) or die $!;
    return ($fh, $path);
}

sub rejected {
    my ($code, $name) = @_;
    my @warnings;
    my $ok;
    my $error;
    {
        local $SIG{__WARN__} = sub { push @warnings, @_ };
        $ok = eval { $code->(); 1 };
        $error = $@;
    }
    ok(!$ok && length($error), $name) or diag('unexpected success');
    unlike($error, qr/Undefined subroutine|Can't locate .*\.pm/, "$name uses the actual reader");
    is(scalar(@warnings), 0, "$name has no incidental warnings") or diag(join('', @warnings));
}

my @cases = (
    ['full raw tokens', '(((A:1.00,B:0):4,C):5,(D:-0,E:7):8);', "基因/α\rinside", 7],
    ['fused recovered', '((A:1.00,B:0):5,(D:6,E:7):8);', 'g/0', 8],
    ['fused unrecovered', '((A:1,D:2):5,(B:6,E:7):8);', 'g-topology', 9],
    ['two retained taxa', '(A:1.00,E:2);', '0', 10],
    ['unicode tab-free gene', '(((A:1e+0,B:.5):4,C:+0):5,(D:6,E:7):8);', '猫\literal_.._%', 11],
);
my ($base_raw, $base_entry, $base_data, $base_record, $base_digest);
my ($fused_raw, $fused_entry, $fused_data);
my ($unrecovered_raw, $unrecovered_entry, $unrecovered_data);
for my $case (@cases) {
    my ($name, $tree, $gene, $number) = @{$case};
    my ($record, $data, $digest) = fixture($tree, $gene, $number);
    my $source_json = $json->encode($record);
    my ($fh, $path) = file_with('prefix');
    seek($fh, 0, 2) or die $!;
    my $entry = write_spool_block($fh, $record, $data, $ctx, $digest, $limits);
    is($json->encode($record), $source_json, "$name writer preserves source values and JSON scalar types");
    my $actual = read_spool_block($fh, $entry, $ctx, $limits);
    is($entry->{offset}, 6, "$name exact nonzero byte offset");
    is($entry->{gene_id}, $gene, "$name original gene identity");
    is($entry->{ledger_sha256}, $digest, "$name source block binding");
    for my $field (qw(primitive state truth groups owners)) {
        is_deeply($actual->{$field}, $data->{$field}, "$name $field roundtrip");
    }
    for my $pair ([no => 'primitive'], [state => 'state'], [truth => 'truth']) {
        my ($file, $field) = @{$pair};
        is($actual->{raw_rows}{$file}, encode_utf8_strict(join("\t", $gene, @{$data->{$field}}) . "\n"),
            "$name $file exact UTF8/LF row");
    }
    my @groups = map { encode_utf8_strict(join("\t", $gene, @{$_}{qw(split members status length)}) . "\n") }
        @{$data->{groups}};
    is_deeply($actual->{raw_rows}{groups}, \@groups, "$name exact sorted group rows");
    for my $view (qw(whole direct)) {
        my @target = @{$ctx->{axis_ids}};
        is_deeply([construct_cell(\@target, $actual->{owners}, $view)],
            [construct_cell(\@target, $data->{owners}, $view)], "$name $view full-target construction unchanged");
        for my $group (@{$data->{groups}}) {
            my @members = split(/\|/, $group->{members});
            is_deeply([construct_cell(\@members, $actual->{owners}, $view)],
                [construct_cell(\@members, $data->{owners}, $view)], "$name $view active construction unchanged");
        }
    }
    seek($fh, $entry->{offset}, 0) or die $!;
    read($fh, my $raw, $entry->{bytes}) == $entry->{bytes} or die $!;
    is(sha256_hex($raw), $entry->{sha256}, "$name seals emitted bytes");
    if ($name eq 'full raw tokens') {
        ($base_raw, $base_entry, $base_data, $base_record, $base_digest) = ($raw, $entry, $data, $record, $digest);
        ok(grep({ $_ eq '1.00' } @{$actual->{primitive}}), 'raw 1.00 remains unchanged');
        ok(grep({ $_ eq '-0' } @{$actual->{primitive}}), 'signed raw zero remains unchanged');
    }
    ($fused_raw, $fused_entry, $fused_data) = ($raw, $entry, $data) if $name eq 'fused recovered';
    ($unrecovered_raw, $unrecovered_entry, $unrecovered_data) = ($raw, $entry, $data) if $name eq 'fused unrecovered';
    close $fh or die $!;
}
ok(grep({ $_->{status} eq 'not_displayed' } @{$unrecovered_data->{groups}}), 'fixture contains an unrecovered fused group');
ok(grep({ $_ eq 'NA_fuse' } @{$unrecovered_data->{truth}}), 'unrecovered group retains NA_fuse truth');
ok(grep({ $_ eq 'U' } @{$unrecovered_data->{state}}), 'unrecovered group carries U state');

sub mutate {
    my ($raw, $entry, $change, $rebind, $name) = @_;
    my $changed = $raw;
    $change->(\$changed);
    my %replacement = %{$entry};
    if ($rebind) {
        $replacement{sha256} = sha256_hex($changed);
        $replacement{bytes} = length($changed);
    }
    my ($fh, $path) = file_with(('x' x $entry->{offset}) . $changed);
    rejected(sub { read_spool_block($fh, \%replacement, $ctx, $limits) }, $name);
    close $fh or die $!;
}

mutate($base_raw, $base_entry, sub { ${$_[0]} =~ s/1\.00/9.00/ or die 'mutation absent' }, 0, 'unsealed byte change');
mutate($base_raw, $base_entry, sub { chop ${$_[0]} }, 0, 'short final block');
mutate($base_raw, $base_entry, sub { chop ${$_[0]} }, 1, 'missing final newline despite rebound SHA');
mutate($base_raw, $base_entry, sub { ${$_[0]} .= "Z\textra\n" }, 1, 'unknown row kind');
mutate($base_raw, $base_entry, sub { ${$_[0]} =~ s/\AH\t/H\twrong\t/ }, 1, 'header wrong field count');
mutate($base_raw, $base_entry, sub { ${$_[0]} =~ s/SplitAligner-materialized-spool-v1/SplitAligner-materialized-spool-v999/ }, 1, 'unknown spool version');
mutate($base_raw, $base_entry, sub { ${$_[0]} =~ s/\AH([^\n]*)\n/H$1\nH$1\n/ }, 1, 'duplicate header');
mutate($base_raw, $base_entry, sub { ${$_[0]} =~ s/^N[^\n]*\n//m or die 'mutation absent' }, 1, 'missing primitive row');
mutate($base_raw, $base_entry, sub { ${$_[0]} =~ s/^(N[^\n]*\n)/$1$1/m or die 'mutation absent' }, 1, 'duplicate primitive row');
mutate($base_raw, $base_entry, sub { ${$_[0]} =~ s/^(N[^\n]*)\n/$1\textra\n/m or die 'mutation absent' }, 1, 'wrong primitive width');
mutate($base_raw, $base_entry, sub { ${$_[0]} =~ s/^(N[^\n]*?)1\.00/$1not-a-number/m or die 'mutation absent' }, 1, 'invalid primitive numeric token');
mutate($base_raw, $base_entry, sub { ${$_[0]} =~ s/^(S\t[^\t]*\t)[^\t\n]+/$1X/m or die 'mutation absent' }, 1, 'invalid state code');
mutate($base_raw, $base_entry, sub { ${$_[0]} =~ s/^(T\t[^\t]*\t)[^\t\n]+/$1MappedWrong/m or die 'mutation absent' }, 1, 'invalid truth token');
mutate($base_raw, $base_entry, sub { ${$_[0]} =~ s/^(C[^\n]*\n)/$1$1/m or die 'mutation absent' }, 1, 'duplicate owner row');
mutate($base_raw, $base_entry, sub { ${$_[0]} =~ s/^C[^\n]*\n//m or die 'mutation absent' }, 1, 'missing owner for non-S branch');
mutate($base_raw, $base_entry, sub { ${$_[0]} =~ s/^(C\t)[^\t]+/$1B99999/m or die 'mutation absent' }, 1, 'owner member absent from axis');
mutate($base_raw, $base_entry, sub { ${$_[0]} =~ s/^C\tB1\t/C\tB2\t/m or die 'mutation absent' }, 1, 'owners overlap on a reference member');
mutate($base_raw, $base_entry, sub { ${$_[0]} =~ s/^(C\t[^\t]+\t)[01]/${1}2/m or die 'mutation absent' }, 1, 'owner recovery is not Boolean');
mutate($base_raw, $base_entry, sub { ${$_[0]} =~ s/^(C\t[^\t]+\t[01]\t)[^\n]+/$1Infinity/m or die 'mutation absent' }, 1, 'nonfinite owner length');
mutate($base_raw, $base_entry, sub { ${$_[0]} =~ s/^(C\t[^\t]+\t1\t)1\.00/${1}9.00/m or die 'mutation absent' }, 1, 'owner length disagrees with primitive row');
mutate($fused_raw, $fused_entry, sub { ${$_[0]} =~ s/^G[^\n]*\n//m or die 'mutation absent' }, 1, 'missing group for fused owner');
mutate($fused_raw, $fused_entry, sub { ${$_[0]} =~ s/^(G[^\n]*\n)/$1$1/m or die 'mutation absent' }, 1, 'duplicate group row');
mutate($fused_raw, $fused_entry, sub { ${$_[0]} =~ s/^(G[^\n]*\t)([^\t\n]+)\n/${1}999\n/m or die 'mutation absent' }, 1, 'group length disagrees with owner');
mutate($fused_raw, $fused_entry, sub { ${$_[0]} =~ s/^(G[^\n]*)\trecovered\t/$1\tnot_displayed\t/m or die 'mutation absent' }, 1, 'group recovery disagrees with owner');
mutate($fused_raw, $fused_entry, sub { ${$_[0]} =~ s/^(G\t[^\t]*\t[^\t]*\t)[^\t]+/${1}B1/m or die 'mutation absent' }, 1, 'group members disagree with complete owner');
mutate($unrecovered_raw, $unrecovered_entry, sub { ${$_[0]} =~ s/^(C\t[^\t]+\t0\t)NA/${1}9/m or die 'mutation absent' }, 1, 'unrecovered owner cannot claim numeric length');
mutate($base_raw, $base_entry, sub { ${$_[0]} =~ s/^(N\t)/$1\xFF/m or die 'mutation absent' }, 1, 'malformed UTF8 scientific row');

for my $change (
    ['offset-negative', offset => -1], ['offset-fraction', offset => '1.5'],
    ['bytes-zero', bytes => 0], ['bytes-negative', bytes => -1], ['bytes-fraction', bytes => '1.5'],
    ['gene binding', gene_id => 'other-gene'], ['record binding', input_record_id => 100],
    ['line binding', input_line => 100], ['reference binding', reference_id => ('0' x 64)],
    ['ledger binding', ledger_sha256 => ('1' x 64)], ['group count binding', groups => 99],
    ['component count binding', components => 99],
) {
    my ($name, $field, $value) = @{$change};
    my %entry = (%{$base_entry}, $field => $value);
    my ($fh, $path) = file_with(('x' x $base_entry->{offset}) . $base_raw);
    rejected(sub { read_spool_block($fh, \%entry, $ctx, $limits) }, $name);
    close $fh or die $!;
}
{
    my %small = (%{$limits}, spool_block_bytes => $base_entry->{bytes} - 1);
    my ($fh, $path) = file_with(('x' x $base_entry->{offset}) . $base_raw);
    rejected(sub { read_spool_block($fh, $base_entry, $ctx, \%small) }, 'reader block capacity');
    close $fh or die $!;
    ($fh, $path) = file_with('');
    rejected(sub { write_spool_block($fh, $base_record, $base_data, $ctx, $base_digest, \%small) }, 'writer block capacity');
    close $fh or die $!;
}
done_testing();
