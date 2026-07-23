#!/usr/bin/env perl

use strict;
use warnings;
use utf8;
use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($RealBin);
use Test::More;
use lib File::Spec->catdir($RealBin, '..', '..', 'lib');
use SplitAligner::IOSafety qw(
    assert_publication_device_records
    build_destination_record
    preflight_publication_devices
    verify_publication_device_guard
    verify_publication_device_records
);
use SplitAligner::TextIO qw(abs_path_utf8);

my $transaction = synthetic_record(
    role => 'transaction parent', path => '/tmp/run',
    publication_parent => '/tmp/run', parent_device => 19,
    existing_destination_path => '/tmp/run');
my $object_diverges = synthetic_record(
    role => 'existing output', path => '/tmp/run/result.tsv',
    publication_parent => '/tmp/run', parent_device => 19,
    existing_destination_path => '/tmp/run/result.tsv', object_device => 17);

lives(
    'matching publication parents accept an existing object on a divergent device',
    sub {
        assert_publication_device_records(
            mode => 'matrix', transaction => $transaction,
            destinations => [$object_diverges]);
    },
);
is($object_diverges->{existing_object_device}, 17,
    'synthetic divergent existing-object device is retained as non-authoritative context');

my $parent_diverges = synthetic_record(
    role => 'remote output', path => '/dev/shm/result.tsv',
    publication_parent => '/dev/shm', parent_device => 24,
    existing_destination_path => '/dev/shm/result.tsv', object_device => 19);
my $error = capture_error(sub {
    assert_publication_device_records(
        mode => 'finalize', transaction => $transaction,
        destinations => [$parent_diverges]);
});
like($error, qr/\[ERROR\]\[Publication device\]/,
    'publication-parent mismatch is rejected');
like($error, qr/Transaction parent '\/tmp\/run'.*device 19/s,
    'diagnostic names the transaction parent and device');
like($error, qr/publication parent '\/dev\/shm'.*device 24/s,
    'diagnostic names the destination publication parent and device');
like($error, qr/existing destination path: '\/dev\/shm\/result\.tsv'/,
    'diagnostic distinguishes the existing destination path');
unlike($error, qr/device 19.*existing destination.*device 19/s,
    'existing-object device is not presented as publication authority');

my $initial = {
    mode => 'finalize', transaction => $transaction,
    destinations => [$object_diverges],
};
my $retargeted_parent = synthetic_record(
    role => 'existing output', path => '/tmp/run/result.tsv',
    publication_parent => '/tmp/run', resolved_parent => '/other/run',
    parent_device => 19, existing_destination_path => '/tmp/run/result.tsv');
rejects(
    'same-device publication-parent retargeting is rejected at recheck',
    sub {
        verify_publication_device_records(
            initial => $initial,
            current => {
                mode => 'finalize', transaction => $transaction,
                destinations => [$retargeted_parent],
            });
    },
    qr/reason: destination-publication-parent-changed-after-preflight/,
);

my $tmp = tempdir(
    'splitaligner-recert008-property-XXXXXX', TMPDIR => 1, CLEANUP => 1);
my $unicode_root = File::Spec->catdir($tmp, '事务_父目录');
my $real_parent = File::Spec->catdir($unicode_root, '真实_父目录');
my $linked_parent = File::Spec->catdir($unicode_root, '链接_父目录');
make_path($real_parent);
ok(symlink($real_parent, $linked_parent),
    'create symlinked Unicode publication parent');
my $existing_destination = File::Spec->catfile($linked_parent, '已有结果.tsv');
open my $fh, '>:raw', $existing_destination
    or die "Cannot create $existing_destination: $!\n";
print {$fh} "existing\n";
close $fh or die "Cannot close $existing_destination: $!\n";
my $destination = build_destination_record(
    role => 'Unicode existing destination', path => $existing_destination,
    type => 'file');
my $resolved_real_parent = abs_path_utf8(
    $real_parent, 'RECERT-008 expected real publication parent');
my $guard;
lives(
    'production preflight derives device from an existing destination parent',
    sub {
        $guard = preflight_publication_devices(
            mode => 'matrix', transaction_root => $unicode_root,
            destinations => [$destination]);
    },
);
is($guard->{destinations}[0]{publication_parent_path}, $linked_parent,
    'canonical publication parent excludes the destination object');
is($guard->{destinations}[0]{resolved_publication_parent_path}, $resolved_real_parent,
    'symlinked publication parent resolves to the real directory');
is($guard->{destinations}[0]{existing_destination_path}, $existing_destination,
    'existing destination is recorded separately');
is($guard->{destinations}[0]{device}, $guard->{transaction}{device},
    'publication-parent device matches the transaction parent device');

my $nested_destination = build_destination_record(
    role => 'Nested missing destination',
    path => File::Spec->catfile($linked_parent, '尚未创建', '结果.tsv'),
    type => 'file');
my $nested_guard;
lives(
    'missing destination hierarchy uses the nearest existing publication parent',
    sub {
        $nested_guard = preflight_publication_devices(
            mode => 'matrix', transaction_root => $unicode_root,
            destinations => [$nested_destination]);
    },
);
is($nested_guard->{destinations}[0]{nearest_existing_parent_path}, $linked_parent,
    'nearest existing publication parent is recorded');
is($nested_guard->{destinations}[0]{resolved_publication_parent_path}, $resolved_real_parent,
    'nearest existing symlink parent is resolved without Unicode normalization');
lives(
    'unchanged publication-parent guard passes defense-in-depth recheck',
    sub {
        verify_publication_device_guard(
            guard => $nested_guard, transaction_root => $unicode_root,
            destinations => [$nested_destination]);
    },
);

done_testing();

sub synthetic_record {
    my (%arg) = @_;
    my $parent = $arg{publication_parent};
    my $resolved = $arg{resolved_parent} // $parent;
    return {
        role                             => $arg{role},
        path                             => $arg{path},
        canonical_path                   => $arg{path},
        publication_parent_path          => $parent,
        resolved_publication_parent_path => $resolved,
        resolved_path                    => $resolved,
        existing_destination_path        =>
            ($arg{existing_destination_path} // ''),
        existing_object_device           => $arg{object_device},
        device                           => $arg{parent_device},
    };
}

sub capture_error {
    my ($code) = @_;
    my $ok = eval { $code->(); 1 };
    return '' if $ok;
    return $@;
}

sub lives {
    my ($name, $code) = @_;
    my $error = capture_error($code);
    ok($error eq '', $name) or diag($error);
}

sub rejects {
    my ($name, $code, $pattern) = @_;
    my $error = capture_error($code);
    ok($error ne '', "$name is rejected");
    like($error, $pattern, "$name reports the expected diagnostic");
}
