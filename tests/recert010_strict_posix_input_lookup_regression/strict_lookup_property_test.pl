#!/usr/bin/env perl

use strict;
use warnings;
use Digest::SHA ();
use Fcntl qw(O_NONBLOCK O_RDONLY);
use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($RealBin);
use IO::Socket::UNIX;
use POSIX qw(mkfifo);
use Socket qw(SOCK_STREAM);
use Test::More;
use lib File::Spec->catdir($RealBin, '..', '..', 'lib');
use SplitAligner::IOSafety qw(build_input_record verify_input_snapshots);
use SplitAligner::TextIO qw(abs_path_utf8 write_utf8_file);

my $tmp = tempdir('sa10-XXXXXX', DIR => '/tmp', CLEANUP => 1);

my $strict = File::Spec->catdir($tmp, 'strict-invalid');
make_path(File::Spec->catdir($strict, 'real'));
write_utf8_file(File::Spec->catfile($strict, 'real', 'notdir'), "not-a-directory\n");
write_utf8_file(File::Spec->catfile($strict, 'real', 'input.nwk'), "intended\n");
ok(symlink(File::Spec->catfile('real', 'notdir'), File::Spec->catfile($strict, 'linkfile')),
    'STRICT-01 symlink-to-file fixture created');

for my $case (
    ['STRICT-01 symlink-to-file intermediate', "$strict/linkfile/../input.nwk"],
    ['STRICT-02 regular-file intermediate', "$strict/real/notdir/../input.nwk"],
    ['STRICT-03 file dot', "$strict/real/input.nwk/."],
    ['STRICT-04 file trailing slash', "$strict/real/input.nwk/"],
) {
    my ($name, $path) = @{$case};
    direct_open_rejected($name, $path);
    rejects_build($name, $path, qr/Cannot open input path.*strict POSIX lookup/s);
}

my $valid = File::Spec->catdir($tmp, 'valid-directory-symlink');
make_path(File::Spec->catdir($valid, 'real', 'sub'));
write_utf8_file(File::Spec->catfile($valid, 'real', 'input.nwk'), "directory-symlink-target\n");
ok(symlink(File::Spec->catfile('real', 'sub'), File::Spec->catfile($valid, 'linkdir')),
    'STRICT-08 directory symlink fixture created');
assert_record_matches_direct_open(
    'STRICT-08 directory symlink followed by parent',
    "$valid/linkdir/../input.nwk",
    "directory-symlink-target\n",
);

my $unicode = File::Spec->catdir($tmp, 'Unicode-猫');
make_path(File::Spec->catdir($unicode, '真实', '内层'));
write_utf8_file(File::Spec->catfile($unicode, '真实', '输入.nwk'), "Unicode-target\n");
ok(symlink(File::Spec->catfile('真实', '内层'), File::Spec->catfile($unicode, '第二跳')),
    'STRICT-09 Unicode second hop created');
ok(symlink('第二跳', File::Spec->catfile($unicode, '第一跳')),
    'STRICT-09 Unicode first hop created');
assert_record_matches_direct_open(
    'STRICT-09 Unicode multiple-hop input',
    "$unicode/第一跳/../../真实/输入.nwk",
    "Unicode-target\n",
);

my $ordinary = File::Spec->catdir($tmp, 'ordinary-parent');
make_path(File::Spec->catdir($ordinary, 'a'));
write_utf8_file(File::Spec->catfile($ordinary, 'input.nwk'), "ordinary\n");
assert_record_matches_direct_open(
    'STRICT valid ordinary a/../b', "$ordinary/a/../input.nwk", "ordinary\n");

my $nonregular = File::Spec->catdir($tmp, 'nonregular');
make_path(File::Spec->catdir($nonregular, 'directory'));
rejects_build('STRICT directory input', File::Spec->catdir($nonregular, 'directory'),
    qr/not-regular-file/);

my $fifo = File::Spec->catfile($nonregular, 'fifo');
ok(mkfifo($fifo, 0600), 'STRICT FIFO fixture created');
rejects_build('STRICT FIFO input', $fifo, qr/not-regular-file/);

rejects_build('STRICT device input', '/dev/null', qr/not-regular-file/);

my $socket_path = File::Spec->catfile($nonregular, 'socket');
SKIP: {
    my $socket = IO::Socket::UNIX->new(
        Type   => SOCK_STREAM,
        Local  => $socket_path,
        Listen => 1,
    );
    skip "Unix socket fixture unavailable: $!", 5 unless defined $socket;
    pass('STRICT Unix socket fixture created');
    rejects_build('STRICT socket input', $socket_path,
        qr/(?:direct-open-failed|not-regular-file)/);
}

my $permission_root = File::Spec->catdir($tmp, 'permission-denied');
make_path($permission_root);
write_utf8_file(File::Spec->catfile($permission_root, 'input.nwk'), "denied\n");
SKIP: {
    skip 'permission denial is not meaningful for root', 2 if $> == 0;
    chmod 0000, $permission_root or die "cannot restrict permission fixture: $!\n";
    direct_open_rejected('STRICT permission-denied intermediate',
        File::Spec->catfile($permission_root, 'input.nwk'));
    rejects_build('STRICT permission-denied input',
        File::Spec->catfile($permission_root, 'input.nwk'), qr/direct-open-failed/);
    chmod 0700, $permission_root or die "cannot restore permission fixture: $!\n";
}

my $retarget = File::Spec->catdir($tmp, 'retarget');
make_path(File::Spec->catdir($retarget, 'first', 'sub'));
make_path(File::Spec->catdir($retarget, 'second', 'sub'));
write_utf8_file(File::Spec->catfile($retarget, 'first', 'input.nwk'), "first\n");
write_utf8_file(File::Spec->catfile($retarget, 'second', 'input.nwk'), "second\n");
my $link = File::Spec->catfile($retarget, 'link');
ok(symlink(File::Spec->catfile('first', 'sub'), $link), 'STRICT-10 retarget link created');
my $provided = "$retarget/link/../input.nwk";
my $record = build_input_record(role => 'STRICT-10 retarget input', path => $provided);
ok(unlink($link), 'STRICT-10 initial link removed');
ok(symlink(File::Spec->catfile('second', 'sub'), $link), 'STRICT-10 link retargeted');
my $snapshot_ok = eval { verify_input_snapshots([$record]); 1 };
ok(!$snapshot_ok, 'STRICT-10 retargeted path rejected by strict revalidation');
like($@, qr/resolved-path-changed/, 'STRICT-10 retarget reports resolved path change');

done_testing();

sub direct_open_rejected {
    my ($name, $path) = @_;
    my $fh;
    my $ok = sysopen($fh, $path, O_RDONLY | O_NONBLOCK);
    close($fh) if $ok;
    ok(!$ok, "$name direct kernel lookup rejects exact path");
}

sub rejects_build {
    my ($name, $path, $pattern) = @_;
    my $ok = eval { build_input_record(role => $name, path => $path); 1 };
    ok(!$ok, "$name fails closed");
    like($@, qr/\Q$name\E/s, "$name diagnostic preserves the role");
    like($@, qr/\Q$path\E/s, "$name diagnostic preserves the provided path");
    like($@, $pattern, "$name diagnostic reports the lookup reason");
}

sub direct_identity {
    my ($path) = @_;
    my $fh;
    sysopen($fh, $path, O_RDONLY | O_NONBLOCK)
        or die "direct open failed for valid fixture $path: $!\n";
    binmode($fh, ':raw');
    my @stat = stat($fh);
    my $sha = Digest::SHA->new(256);
    $sha->addfile($fh);
    close($fh);
    return {
        device => 0 + $stat[0],
        inode  => 0 + $stat[1],
        size   => 0 + $stat[7],
        sha256 => $sha->hexdigest,
    };
}

sub assert_record_matches_direct_open {
    my ($name, $path, $expected_content) = @_;
    my $direct = direct_identity($path);
    my $resolved = abs_path_utf8($path, "$name realpath");
    ok(defined $resolved, "$name canonical path available after strict open");
    my $record = build_input_record(role => $name, path => $path);
    is($record->{path_as_provided}, $path, "$name preserves path_as_provided");
    is($record->{resolved_path}, $resolved, "$name resolved path matches OS realpath");
    is($record->{canonical_path}, $resolved, "$name canonical path identifies resolved object");
    is($record->{sha256}, $direct->{sha256}, "$name SHA binds direct-open handle");
    is($record->{size}, $direct->{size}, "$name size binds direct-open handle");
    is($record->{device}, $direct->{device}, "$name device binds direct-open handle");
    is($record->{inode}, $direct->{inode}, "$name inode binds direct-open handle");
    open my $fh, '<:raw', $record->{resolved_path} or die "cannot read $resolved: $!\n";
    local $/;
    my $content = <$fh>;
    close $fh;
    is($content, $expected_content, "$name parser-compatible path reads intended object");
    ok(eval { verify_input_snapshots([$record]); 1 }, "$name strict snapshot succeeds")
        or diag($@);
}
