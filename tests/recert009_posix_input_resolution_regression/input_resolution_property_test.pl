#!/usr/bin/env perl

use strict;
use warnings;
use File::Path qw(make_path);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($RealBin);
use Test::More;
use lib File::Spec->catdir($RealBin, '..', '..', 'lib');
use SplitAligner::IOSafety qw(build_input_record verify_input_snapshots);
use SplitAligner::Provenance qw(file_sha256);
use SplitAligner::TextIO qw(abs_path_utf8 write_utf8_file);

my $tmp = tempdir('splitaligner-recert009-property-XXXXXX', TMPDIR => 1, CLEANUP => 1);

my $decoy = make_symlink_parent_fixture(
    File::Spec->catdir($tmp, 'decoy-present'),
    "intended\n", "decoy\n", 1,
);
assert_record_matches_os('PATH-009-01/02 decoy present', $decoy, "intended\n");

my $no_decoy = make_symlink_parent_fixture(
    File::Spec->catdir($tmp, 'decoy-absent'),
    "intended-only\n", undef, 0,
);
assert_record_matches_os('PATH-009-03/04 decoy absent', $no_decoy, "intended-only\n");

my $unicode_root = File::Spec->catdir($tmp, 'Unicode-猫');
make_path(File::Spec->catdir($unicode_root, '真实', '子'));
write_utf8_file(File::Spec->catfile($unicode_root, '真实', '输入.nwk'), "Unicode-target\n");
ok(symlink(File::Spec->catfile('真实', '子'), File::Spec->catfile($unicode_root, '链接')),
    'PATH-009-11 Unicode symlink fixture created');
my $unicode_provided = "$unicode_root/链接/../输入.nwk";
assert_record_matches_os('PATH-009-11 Unicode components', $unicode_provided, "Unicode-target\n");

my $multi = File::Spec->catdir($tmp, 'multiple-hops');
make_path(File::Spec->catdir($multi, 'target', 'inner'));
write_utf8_file(File::Spec->catfile($multi, 'target', 'input.nwk'), "multiple-hop-target\n");
ok(symlink(File::Spec->catfile('target', 'inner'), File::Spec->catfile($multi, 'hop2')),
    'PATH-009-12 second hop created');
ok(symlink('hop2', File::Spec->catfile($multi, 'hop1')),
    'PATH-009-12 first hop created');
my $multi_provided = "$multi/hop1/../../target/input.nwk";
assert_record_matches_os('PATH-009-12 multiple hops and parents', $multi_provided,
    "multiple-hop-target\n");

my $normal = File::Spec->catdir($tmp, 'normal-parent');
make_path(File::Spec->catdir($normal, 'a'));
write_utf8_file(File::Spec->catfile($normal, 'b.nwk'), "normal-parent\n");
assert_record_matches_os('PATH-009-13 ordinary a/../b', "$normal/a/../b.nwk",
    "normal-parent\n");

my $dangling = File::Spec->catfile($tmp, 'dangling-input');
ok(symlink('missing-target', $dangling), 'PATH-009-14 dangling fixture created');
rejects_build('PATH-009-14 dangling symlink', $dangling,
    qr/(?:Cannot resolve input path|Resolved input is not a regular file)/);

my $loop_a = File::Spec->catfile($tmp, 'loop-a');
my $loop_b = File::Spec->catfile($tmp, 'loop-b');
ok(symlink('loop-b', $loop_a), 'PATH-009-14 loop first link created');
ok(symlink('loop-a', $loop_b), 'PATH-009-14 loop second link created');
rejects_build('PATH-009-14 symlink loop', $loop_a, qr/Cannot resolve input path/);

my $retarget = File::Spec->catdir($tmp, 'retarget');
make_path(File::Spec->catdir($retarget, 'first', 'sub'));
make_path(File::Spec->catdir($retarget, 'second', 'sub'));
write_utf8_file(File::Spec->catfile($retarget, 'first', 'input.nwk'), "first\n");
write_utf8_file(File::Spec->catfile($retarget, 'second', 'input.nwk'), "second\n");
my $retarget_link = File::Spec->catfile($retarget, 'link');
ok(symlink(File::Spec->catfile('first', 'sub'), $retarget_link),
    'PATH-009-15 initial symlink created');
my $retarget_provided = "$retarget/link/../input.nwk";
my $retarget_record = build_input_record(role => 'retargeted input', path => $retarget_provided);
ok(unlink($retarget_link), 'PATH-009-15 initial symlink removed');
ok(symlink(File::Spec->catfile('second', 'sub'), $retarget_link),
    'PATH-009-15 symlink retargeted');
my $snapshot_ok = eval { verify_input_snapshots([$retarget_record]); 1 };
ok(!$snapshot_ok, 'PATH-009-15 retargeted provided path is rejected');
like($@, qr/reason: resolved-path-changed/,
    'PATH-009-15 retarget reports resolved-path change');

rejects_build('POSIX input NUL rejection', "bad\x00path",
    qr/Filesystem path contains a NUL byte/);

done_testing();

sub make_symlink_parent_fixture {
    my ($root, $intended, $decoy_content, $with_decoy) = @_;
    make_path(File::Spec->catdir($root, 'real', 'sub'));
    write_utf8_file(File::Spec->catfile($root, 'real', 'input.nwk'), $intended);
    write_utf8_file(File::Spec->catfile($root, 'input.nwk'), $decoy_content)
        if $with_decoy;
    die "cannot create symlink fixture: $!\n"
        unless symlink(File::Spec->catfile('real', 'sub'), File::Spec->catfile($root, 'link'));
    return "$root/link/../input.nwk";
}

sub assert_record_matches_os {
    my ($name, $provided, $expected_content) = @_;
    my $os_resolved = abs_path_utf8($provided, "$name OS realpath");
    ok(defined $os_resolved, "$name resolves under OS semantics");
    my $record = build_input_record(role => $name, path => $provided);
    is($record->{path_as_provided}, $provided, "$name preserves path_as_provided");
    is($record->{resolved_path}, $os_resolved, "$name resolved_path matches OS realpath");
    is($record->{canonical_path}, $os_resolved, "$name canonical input path is resolved object");
    is($record->{sha256}, file_sha256($os_resolved), "$name SHA binds resolved object");
    is($record->{size}, 0 + (-s $os_resolved), "$name size binds resolved object");
    my @stat = stat($os_resolved);
    is($record->{device}, 0 + $stat[0], "$name device binds resolved object");
    is($record->{inode}, 0 + $stat[1], "$name inode binds resolved object");
    open my $fh, '<:raw', $record->{resolved_path} or die "cannot read $record->{resolved_path}: $!\n";
    local $/;
    my $content = <$fh>;
    close $fh;
    is($content, $expected_content, "$name reads intended object");
    ok(eval { verify_input_snapshots([$record]); 1 }, "$name snapshot revalidation succeeds")
        or diag($@);
}

sub rejects_build {
    my ($name, $path, $pattern) = @_;
    my $ok = eval { build_input_record(role => $name, path => $path); 1 };
    ok(!$ok, "$name fails closed");
    like($@, $pattern, "$name reports resolution failure");
}
