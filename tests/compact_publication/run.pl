#!/usr/bin/env perl
use strict;
use warnings;
use Test::More;
use Cwd qw(abs_path);
use Digest::SHA qw(sha256_hex);
use File::Path qw(make_path remove_tree);
use File::Spec;
use File::Temp qw(tempdir);
use FindBin qw($RealBin);
use lib "$RealBin/../../lib";
use SplitAligner::OutputOwnership qw(
    attach_output_ownership build_publication_plan object_matches_snapshot
    path_exists preflight_output_ownership verify_publication_plan
    verify_staged_replacement
);
use SplitAligner::IOSafety qw(
    assert_io_namespaces_disjoint build_destination_record
    preflight_publication_devices verify_input_snapshots
    verify_publication_device_guard
);
use SplitAligner::Provenance qw(write_json_file);
use SplitAligner::TextIO qw(decode_filesystem_path_utf8 getcwd_utf8);

# Exercise the actual controller transaction bodies without starting matrix
# construction. No production test hooks or copies of transaction logic.
my $root = abs_path("$RealBin/../..");
open my $source_fh, '<:raw', "$root/SplitAligner.pl" or die $!;
my $controller = do { local $/; <$source_fh> };
close $source_fh or die $!;
for my $name (qw(publish_transaction publication_test_fail_after remove_path)) {
    my ($body) = $controller =~ /^(sub \Q$name\E\b.*?)(?=^sub |\z)/ms;
    die "Cannot locate production sub $name" unless defined $body;
    eval $body;
    die "Cannot load production sub $name: $@" if $@;
}
my $work = @ARGV ? File::Spec->rel2abs(shift @ARGV) : tempdir(CLEANUP => 1);
if (-d $work) {
    opendir my $existing, $work or die $!;
    my @contents = grep { $_ ne '.' && $_ ne '..' } readdir $existing;
    closedir $existing or die $!;
    die "Refusing to overwrite evidence directory $work" if @contents;
}
make_path($work);
diag("Fixture directory: $work");
write_json_file("$work/source_binding.json", {
    controller_sha256 => sha256_hex($controller),
    output_ownership_sha256 => sha256_hex(read_bytes("$root/lib/SplitAligner/OutputOwnership.pm")),
    provenance_sha256 => sha256_hex(read_bytes("$root/lib/SplitAligner/Provenance.pm")),
    driver_sha256 => sha256_hex(read_bytes("$RealBin/run.pl")),
    scope => 'Actual publication functions and module; no matrix or benchmark run',
});

sub read_bytes {
    my ($p) = @_;
    open my $fh, '<:raw', $p or die "$p: $!";
    my $s = do { local $/; <$fh> };
    close $fh or die $!;
    return $s;
}
sub put {
    my ($p, $s) = @_;
    open my $fh, '>:raw', $p or die "$p: $!";
    print {$fh} $s or die $!;
    close $fh or die $!;
}
sub prepare {
    my ($f, $generation, $stage_name) = @_;
    my $stage = "$f->{base}/$stage_name";
    make_path("$stage/owned");
    put("$stage/value.txt", "$generation value\n");
    put("$stage/owned/member.txt", "$generation member\n");
    put("$stage/shared.txt", "immutable shared\n");
    my $manifest = { schema_version => 'publication-test-v1', generation => $generation };
    attach_output_ownership($manifest,
        mode => 'matrix', owner_label => 'p', source_root => $stage,
        owner_manifest_path => "$f->{out}/p.run_manifest.json", items => $f->{items});
    my $seal = write_json_file("$stage/p.run_manifest.json", $manifest);
    return { %$f, stage => $stage, inventory => $manifest->{ownership}{outputs}, manifest_seal => $seal };
}
sub transaction {
    my ($f, $force) = @_;
    my $auth = preflight_output_ownership(
        mode => 'matrix', owner_label => 'p',
        owner_manifest_path => "$f->{out}/p.run_manifest.json",
        items => $f->{items}, force => $force);
    my @dest = map { build_destination_record(
        role => $_->{role}, path => $_->{destination}, type => $_->{type}) } @{$f->{items}};
    my $guard = preflight_publication_devices(
        mode => 'matrix', transaction_root => $f->{out}, destinations => \@dest);
    my $cwd = getcwd_utf8();
    chdir $f->{out} or die $!;
    my $ok = eval {
        publish_transaction($f->{stage}, $f->{items}, $force, [], \@dest, $auth, $guard, $f->{inventory}, $f->{manifest_seal});
        1;
    };
    my $err = $@;
    chdir $cwd or die $!;
    die $err unless $ok;
    return 1;
}
sub fixture {
    my ($name, $generation, $clean) = @_;
    my $base = "$work/$name";
    die "Refusing to overwrite fixture $base" if -e $base;
    make_path("$base/out");
    my $f = { base => $base, out => "$base/out" };
    $f->{items} = [
        { source => 'value.txt', destination => "$base/out/value.txt", type => 'file', role => 'value' },
        { source => 'owned', destination => "$base/out/owned", type => 'directory', role => 'owned records' },
        { source => 'shared.txt', destination => "$base/out/shared.txt", type => 'file', role => 'shared', shared => 1 },
        { source => 'p.run_manifest.json', destination => "$base/out/p.run_manifest.json",
            type => 'file', role => 'manifest', owner_manifest => 1 },
    ];
    transaction(prepare($f, 'old', 'old-stage'), 0) unless $clean;
    return prepare($f, $generation // 'new', 'next-stage');
}
sub snapshot {
    my ($f) = @_;
    return { map {
        my $p = "$f->{out}/$_";
        my @st = lstat $p;
        $_ => @st ? [sha256_hex(read_bytes($p)), @st[0,1,2,3]] : undef;
    } qw(value.txt owned/member.txt shared.txt p.run_manifest.json) };
}
sub backups {
    my ($f) = @_;
    opendir my $dh, $f->{out} or die $!;
    my @d = map { "$f->{out}/$_" } grep { /^\.splitaligner-publish-backup-/ } readdir $dh;
    closedir $dh or die $!;
    return @d;
}
sub attempt {
    my ($f) = @_;
    my $ok = eval { transaction($f, 1); 1 };
    my $err = $@;
    put("$f->{base}/error.txt", $err);
    return ($ok, $err);
}

{
    my $f = fixture('normal');
    ok(transaction($f, 1), 'normal authorized replacement succeeds');
    is(read_bytes("$f->{out}/value.txt"), "new value\n", 'normal replacement publishes new bytes');
    is(scalar(backups($f)), 0, 'normal replacement cleans obsolete backups');
}
for my $target ('value.txt', 'owned/member.txt') {
    my $case = $target eq 'value.txt' ? 'P1_file' : 'P1_directory';
    my $f = fixture($case);
    my $before = snapshot($f);
    put("$f->{stage}/$target", "changed after manifest\n");
    my ($ok, $err) = attempt($f);
    ok(!$ok, "$case fails");
    like($err, qr/differs from its manifest inventory/, "$case compares manifest authority");
    is_deeply(snapshot($f), $before, "$case keeps every prior output unchanged");
    is(scalar(backups($f)), 0, "$case makes no backup or publication move");
}
{
    my $f = fixture('P1_clean', 'new', 1);
    my $before = snapshot($f);
    put("$f->{stage}/value.txt", "changed after manifest\n");
    my ($ok, $err) = attempt($f);
    ok(!$ok, 'P1 clean run rejects post-manifest generated change');
    like($err, qr/differs from its manifest inventory/, 'P1 clean run uses manifest authority');
    is_deeply(snapshot($f), $before, 'P1 clean run publishes nothing');
}
for my $after_plan (0, 1) {
    my $name = $after_plan ? 'manifest_after_plan' : 'manifest_before_plan';
    my $f = fixture($name);
    my $before = snapshot($f);
    my $manifest_path = "$f->{stage}/p.run_manifest.json";
    is($f->{manifest_seal}, sha256_hex(read_bytes($manifest_path)), "$name writer seal matches emitted bytes");
    my $mutate = sub { put($manifest_path, read_bytes($manifest_path) . " \n"); };
    my $real = \&build_publication_plan;
    no warnings 'redefine';
    local *build_publication_plan = sub {
        my $plan = $real->(@_);
        $mutate->() if $after_plan;
        return $plan;
    };
    $mutate->() unless $after_plan;
    my ($ok, $err) = attempt($f);
    ok(!$ok, "$name rejects changed manifest bytes even when JSON remains equivalent");
    like($err, $after_plan ? qr/Generated transaction source changed/ : qr/serialization seal/,
        "$name identifies the sealed-object failure");
    is_deeply(snapshot($f), $before, "$name leaves all prior output unchanged");
    is(scalar(backups($f)), 0, "$name fails before publication moves");
}
{
    my $f = fixture('inventory_extra');
    push @{$f->{inventory}}, { %{$f->{inventory}[0]}, relative_path => 'unexpected' };
    my $before = snapshot($f);
    my ($ok, $err) = attempt($f);
    ok(!$ok, 'extra inventory entry fails');
    like($err, qr/unexpected publication paths/, 'extra inventory entry identified');
    is_deeply(snapshot($f), $before, 'extra inventory leaves old outputs unchanged');
}
for my $clean (0, 1) {
    my $name = $clean ? 'P2_clean_source_after_plan' : 'P4_replace_source_after_plan';
    my $f = fixture($name, 'new', $clean);
    my $before = snapshot($f);
    my $real = \&build_publication_plan;
    no warnings 'redefine';
    local *build_publication_plan = sub {
        my $plan = $real->(@_);
        put("$f->{stage}/value.txt", "changed after plan\n");
        return $plan;
    };
    my ($ok, $err) = attempt($f);
    ok(!$ok, "$name post-plan source mutation fails");
    like($err, qr/Generated transaction source changed/, "$name still checks content before publication");
    is_deeply(snapshot($f), $before, "$name preserves destination state");
}
for my $case (['P6_owned_reuse', 'old', 'value.txt'], ['P7_shared_reuse', 'old', 'shared.txt'], ['P3_early_replace', 'new', 'value.txt'], ['P8_late_replace', 'new', 'owned/member.txt']) {
    my ($name, $generation, $target) = @$case;
    my $f = fixture($name, $generation);
    my $expected;
    my $real = \&build_publication_plan;
    no warnings 'redefine';
    local *build_publication_plan = sub {
        my $plan = $real->(@_);
        put("$f->{out}/$target", "injected destination change\n");
        $expected = snapshot($f);
        return $plan;
    };
    my ($ok, $err) = attempt($f);
    ok(!$ok, "$name rejects target change");
    is_deeply(snapshot($f), $expected, "$name changes only the injected target");
    is(scalar(backups($f)), 0, "$name rejects before any old output moves");
}
for my $target ('value.txt', 'owned/member.txt') {
    my $case = $target eq 'value.txt' ? 'P5_file' : 'P5_directory';
    my $f = fixture($case);
    my $real = \&verify_publication_device_guard;
    no warnings 'redefine';
    local *verify_publication_device_guard = sub {
        $real->(@_);
        put("$f->{stage}/$target", "changed after all pre-publication checks\n");
        return 1;
    };
    my ($ok, $err) = attempt($f);
    ok(!$ok, "$case late source change fails");
    like($err, qr/Published output differs/, "$case detected after install");
    like($err, qr/prior output retained for recovery at/, "$case locates preserved old output");
    my @backup = backups($f);
    is(scalar @backup, 1, "$case retains recovery directory");
    if (@backup == 1) {
        my $relative = $target eq 'value.txt' ? '000001' : '000002/member.txt';
        is(read_bytes("$backup[0]/$relative"), $target eq 'value.txt' ? "old value\n" : "old member\n",
            "$case old bytes remain recoverable");
        like($err, qr/\Q$backup[0]\E/, "$case reports actual backup path");
    }
    is(read_bytes("$f->{out}/$target"), "changed after all pre-publication checks\n", "$case preserves untrusted installed object");
    is(read_bytes("$f->{out}/shared.txt"), "immutable shared\n", "$case leaves shared output untouched");
}
{
    my $f = fixture('E10_verifier_exception');
    my $before = snapshot($f);
    my $real = \&verify_staged_replacement;
    no warnings 'redefine';
    local *verify_staged_replacement = sub { $real->(@_); die "[INJECTED E10] post-backup verifier exception\n"; };
    my ($ok, $err) = attempt($f);
    ok(!$ok, 'E10 verifier exception fails');
    like($err, qr/INJECTED E10/, 'E10 failure retained');
    is_deeply(snapshot($f), $before, 'E10 exception restores exact prior files and inodes');
    is(scalar(backups($f)), 0, 'E10 complete restoration leaves no backup residue');
}
{
    my $f = fixture('E10_backup_tamper');
    my $real = \&verify_staged_replacement;
    no warnings 'redefine';
    local *verify_staged_replacement = sub {
        my ($item, $path) = @_;
        put($path, "tampered staged backup\n");
        return $real->(@_);
    };
    my ($ok, $err) = attempt($f);
    ok(!$ok, 'E10 actual backup mutation fails');
    like($err, qr/backup differs from its authorized snapshot/, 'E10 rejects automatic restoration of changed backup');
    my @backup = backups($f);
    is(scalar @backup, 1, 'E10 changed backup is retained');
    is(read_bytes("$backup[0]/000001"), "tampered staged backup\n", 'E10 does not delete changed recovery object') if @backup == 1;
    ok(!-e "$f->{out}/value.txt", 'E10 does not claim a changed backup is the original result');
    like($err, qr/Transaction recovery directory retained/, 'E10 reports incomplete recovery');
}
{
    my $f = fixture('ordinary_late_failure');
    my $before = snapshot($f);
    local $ENV{SPLITALIGNER_INTERNAL_TESTING} = '1';
    local $ENV{SPLITALIGNER_TEST_FAIL_AFTER_PUBLISH} = '2';
    my ($ok, $err) = attempt($f);
    ok(!$ok, 'ordinary late failure injected');
    is_deeply(snapshot($f), $before, 'ordinary late failure restores all old outputs');
    is(scalar(backups($f)), 0, 'ordinary late failure has no backup residue');
}
done_testing();
