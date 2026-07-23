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
    assert_io_namespaces_disjoint
    build_destination_record
    build_input_record
    verify_input_snapshots
);
use SplitAligner::TextIO qw(write_utf8_file);

my $tmp = tempdir('splitaligner-recert005-property-XXXXXX', TMPDIR => 1, CLEANUP => 1);
my $data = File::Spec->catdir($tmp, '数据');
my $out = File::Spec->catdir($tmp, '输出');
make_path($data, $out);

my $input_path = File::Spec->catfile($data, '输入.nwk');
write_utf8_file($input_path, "((A:1,B:1):2,(C:1,D:1):3);\n");
my $input = build_input_record(role => 'property input', path => $input_path);

rejects(
    'exact canonical path',
    [$input],
    [destination('exact output', $input_path, 'file')],
    qr/reason: exact-path/,
);

my $dot_alias = File::Spec->catfile($data, '.', '..', '数据', '输入.nwk');
rejects(
    'dot and dot-dot lexical alias',
    [$input],
    [destination('dot alias output', $dot_alias, 'file')],
    qr/reason: exact-path/,
);

my $symlink_path = File::Spec->catfile($out, 'symlink-output.nwk');
if (symlink($input_path, $symlink_path)) {
    rejects(
        'symlink-resolved alias',
        [$input],
        [destination('symlink output', $symlink_path, 'file')],
        qr/reason: resolved-path/,
    );
} else {
    fail("create symlink fixture: $!");
}

my $hardlink_path = File::Spec->catfile($out, 'hardlink-output.nwk');
if (link($input_path, $hardlink_path)) {
    rejects(
        'hardlink inode alias',
        [$input],
        [destination('hardlink output', $hardlink_path, 'file')],
        qr/reason: same-inode/,
    );
} else {
    fail("create hardlink fixture: $!");
}

my $nested_dir = File::Spec->catdir($tmp, 'nested-output');
make_path($nested_dir);
my $nested_input_path = File::Spec->catfile($nested_dir, 'inside.nwk');
write_utf8_file($nested_input_path, "(A:1,B:1);\n");
my $nested_input = build_input_record(role => 'nested input', path => $nested_input_path);
rejects(
    'input nested inside output directory',
    [$nested_input],
    [destination('published directory', $nested_dir, 'directory')],
    qr/reason: inside-output-directory/,
);

my $dir_symlink = File::Spec->catdir($tmp, 'directory-link');
if (symlink($nested_dir, $dir_symlink)) {
    rejects(
        'resolved input nested inside symlinked output directory',
        [$nested_input],
        [destination('symlinked published directory', $dir_symlink, 'directory')],
        qr/reason: inside-output-directory/,
    );
} else {
    fail("create directory symlink fixture: $!");
}

accepts(
    'component-aware comparison does not confuse prefix siblings',
    [$input],
    [destination('safe prefix sibling', File::Spec->catfile($tmp, '数据-extra'), 'directory')],
);

rejects(
    'duplicate file destinations',
    [],
    [
        destination('first output', File::Spec->catfile($out, 'same.txt'), 'file'),
        destination('second output', File::Spec->catfile($out, '.', 'same.txt'), 'file'),
    ],
    qr/reason: duplicate-destination/,
);

rejects(
    'file destination nested under directory destination',
    [],
    [
        destination('output directory', File::Spec->catdir($out, 'tree'), 'directory'),
        destination('nested output file', File::Spec->catfile($out, 'tree', 'value.txt'), 'file'),
    ],
    qr/reason: duplicate-destination/,
);

rejects(
    'overlapping directory destinations',
    [],
    [
        destination('outer output directory', File::Spec->catdir($out, 'outer'), 'directory'),
        destination('inner output directory', File::Spec->catdir($out, 'outer', 'inner'), 'directory'),
    ],
    qr/reason: duplicate-destination/,
);

accepts(
    'disjoint Unicode paths',
    [$input],
    [
        destination('Unicode output file', File::Spec->catfile($out, '结果.tsv'), 'file'),
        destination('Unicode output directory', File::Spec->catdir($out, '分支目录'), 'directory'),
    ],
);

my $mutation_path = File::Spec->catfile($data, 'mutation.txt');
write_utf8_file($mutation_path, "before\n");
my $mutation_input = build_input_record(role => 'mutable input', path => $mutation_path);
write_utf8_file($mutation_path, "after\n");
my $mutation_ok = eval { verify_input_snapshots([$mutation_input]); 1 };
ok(!$mutation_ok, 'pre-publication input mutation is rejected');
like($@, qr/Input role 'mutable input'.*reason: content-sha256/s,
    'mutation diagnostic names role and SHA-256 reason');

done_testing();

sub destination {
    my ($role, $path, $type) = @_;
    return build_destination_record(role => $role, path => $path, type => $type);
}

sub rejects {
    my ($name, $inputs, $destinations, $pattern) = @_;
    my $ok = eval {
        assert_io_namespaces_disjoint(
            inputs => $inputs,
            destinations => $destinations,
        );
        1;
    };
    ok(!$ok, "$name is rejected");
    like($@, $pattern, "$name reports the expected reason");
    like($@, qr/--force cannot override/, "$name cannot be bypassed by --force");
}

sub accepts {
    my ($name, $inputs, $destinations) = @_;
    my $ok = eval {
        assert_io_namespaces_disjoint(
            inputs => $inputs,
            destinations => $destinations,
        );
        1;
    };
    ok($ok, $name) or diag($@);
}
