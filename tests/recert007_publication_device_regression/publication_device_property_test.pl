#!/usr/bin/env perl

use strict;
use warnings;
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

my $transaction = synthetic_record(
    role => 'transaction root', path => '/local/工作', device => 100);
my $local_free = synthetic_record(
    role => 'FREE matrix NA_fuse output',
    path => '/local/工作/free.matrix_with_fuse.na_fuse.txt', device => 100);
my $local_fix = synthetic_record(
    role => 'FIX matrix NA_fuse output',
    path => '/local/工作/fix.matrix_with_fuse.na_fuse.txt', device => 100);
my $remote_free = synthetic_record(
    role => 'FREE matrix NA_fuse output',
    path => '/remote/矩阵/free.matrix_with_fuse.na_fuse.txt', device => 200);
my $remote_fix = synthetic_record(
    role => 'FIX matrix NA_fuse output',
    path => '/remote/矩阵/fix.matrix_with_fuse.na_fuse.txt', device => 200);

lives(
    'same-device synthetic destinations are accepted',
    sub {
        assert_publication_device_records(
            mode => 'finalize', transaction => $transaction,
            destinations => [$local_free, $local_fix]);
    },
);

rejects(
    'paired cross-device layout',
    sub {
        assert_publication_device_records(
            mode => 'finalize', transaction => $transaction,
            destinations => [$remote_free, $remote_fix]);
    },
    qr/\[ERROR\]\[Publication device\].*Mode 'finalize'.*FREE matrix NA_fuse output/s,
);

rejects(
    'fix-only cross-device layout',
    sub {
        assert_publication_device_records(
            mode => 'finalize_fix', transaction => $transaction,
            destinations => [$remote_fix]);
    },
    qr/Mode 'finalize_fix'.*FIX matrix NA_fuse output/s,
);

rejects(
    'mixed-device layout identifies the remote role',
    sub {
        assert_publication_device_records(
            mode => 'finalize', transaction => $transaction,
            destinations => [$local_fix, $remote_free]);
    },
    qr/destination role 'FREE matrix NA_fuse output'/,
);

my $unicode_error = capture_error(sub {
    assert_publication_device_records(
        mode => 'finalize', transaction => $transaction,
        destinations => [$remote_free]);
});
like($unicode_error, qr{/local/工作}, 'diagnostic preserves Unicode transaction path');
like($unicode_error, qr{/remote/矩阵}, 'diagnostic preserves Unicode destination path');
like($unicode_error, qr/device 100.*device 200/s, 'diagnostic reports both devices');
like($unicode_error, qr/--force cannot override/, 'device mismatch cannot be bypassed by --force');
like(
    $unicode_error,
    qr/run on the matrix filesystem or copy the matrix run into the transaction filesystem/,
    'diagnostic reports the required remediation',
);

my $initial = {
    mode => 'finalize', transaction => $transaction,
    destinations => [$local_free, $local_fix],
};
my $changed_free = synthetic_record(
    role => 'FREE matrix NA_fuse output',
    path => '/local/工作/free.matrix_with_fuse.na_fuse.txt', device => 200);
my $current = {
    mode => 'finalize', transaction => $transaction,
    destinations => [$changed_free, $local_fix],
};
rejects(
    'destination device change is rejected by publication recheck',
    sub { verify_publication_device_records(initial => $initial, current => $current) },
    qr/reason: destination-publication-parent-device-changed-after-preflight/,
);

my $changed_transaction = synthetic_record(
    role => 'transaction root', path => '/local/工作', device => 300);
rejects(
    'transaction device change is rejected by publication recheck',
    sub {
        verify_publication_device_records(
            initial => $initial,
            current => {
                mode => 'finalize', transaction => $changed_transaction,
                destinations => [$local_free, $local_fix],
            },
        );
    },
    qr/reason: transaction-device-changed-after-preflight/,
);

my $tmp = tempdir(
    'splitaligner-recert007-property-XXXXXX', TMPDIR => 1, CLEANUP => 1);
my $unicode_root = File::Spec->catdir($tmp, '事务_路径');
my $real_destination_root = File::Spec->catdir($unicode_root, '真实_输出');
my $linked_destination_root = File::Spec->catdir($unicode_root, '链接_输出');
make_path($real_destination_root);
ok(
    symlink($real_destination_root, $linked_destination_root),
    'create symlinked same-device destination parent',
);
my $destination_path = File::Spec->catfile(
    $linked_destination_root, '尚未创建', '结果.tsv');
my $destination = build_destination_record(
    role => 'Unicode symlink destination', path => $destination_path, type => 'file');
my $guard;
lives(
    'real same-device Unicode/symlink destination preflight succeeds',
    sub {
        $guard = preflight_publication_devices(
            mode => 'matrix', transaction_root => $unicode_root,
            destinations => [$destination]);
    },
);
is(
    $guard->{transaction}{device},
    $guard->{destinations}[0]{device},
    'real destination device is derived from the nearest resolved parent',
);
lives(
    'unchanged real publication-device guard passes recheck',
    sub {
        verify_publication_device_guard(
            guard => $guard, transaction_root => $unicode_root,
            destinations => [$destination]);
    },
);

done_testing();

sub synthetic_record {
    my (%arg) = @_;
    my $parent = $arg{publication_parent} // $arg{path};
    my $resolved_parent = $arg{resolved_parent} // $parent;
    return {
        role                             => $arg{role},
        path                             => $arg{path},
        canonical_path                   => $arg{path},
        publication_parent_path          => $parent,
        resolved_publication_parent_path => $resolved_parent,
        resolved_path                    => $resolved_parent,
        existing_destination_path        =>
            ($arg{existing_destination_path} // ''),
        device                           => $arg{device},
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
