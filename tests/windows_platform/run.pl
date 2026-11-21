#!/usr/bin/env perl
use strict;
use warnings;
use utf8;
use FindBin qw($RealBin);
use lib "$RealBin/../../lib";
use Test::More;
binmode Test::More->builder->output, ':utf8';
binmode Test::More->builder->failure_output, ':utf8';
use File::Temp qw(tempdir);
use SplitAligner::Platform qw(windows_absolute_path path_identity_key native_windows assert_native_runtime assert_windows_stat assert_windows_path);
use SplitAligner::IOSafety ();
use SplitAligner::OutputOwnership ();
use File::Spec;

for my $case (
    ['C:\\work\\input.nwk', 'C:/work/input.nwk'],
    ['c:/work/./sub/../input.nwk', 'C:/work/input.nwk'],
    ['input.nwk', 'C:/Work Dir/input.nwk'],
    ['../input.nwk', 'C:/input.nwk'],
    ['C:/', 'C:/'],
    ['D:/folder/file', 'D:/folder/file'],
) {
    is(windows_absolute_path($case->[0], 'C:/Work Dir'), $case->[1], "Windows lexical: $case->[0]");
}
for my $path ('', 'C:input.nwk', '/input.nwk', '\\input.nwk',
    '\\\\server\\share\\file', '\\\\?\\C:\\file', '\\\\.\\NUL',
    'C:/x:stream', 'C:/bad?', 'C:/bad*', 'C:/bad|', 'C:/bad"',
    'C:/NUL', 'C:/nul.txt', 'C:/COM1.log', 'C:/LPT9', 'C:/CON/child',
    'C:/bad.', 'C:/bad ', "C:/bad\tfile", 'C:/输入/file', 'C:/PROGRA~1/out', 'C:/' . ('x' x 240)) {
    my $ok = eval { windows_absolute_path($path, 'C:/Work Dir'); 1 };
    ok(!$ok && $@ =~ /\[Windows\]/, "reject ambiguous/unsupported Windows path: $path");
}
ok(!eval { windows_absolute_path('file', 'C:/用户'); 1 }, 'reject non-ASCII inherited working directory');
{
    no warnings 'redefine';
    local *SplitAligner::Platform::native_windows = sub { 1 };
    is(path_identity_key('C:\\Work\\TREE.NWK'), 'c:/work/tree.nwk', 'comparison folds Windows case and separator');
    is(SplitAligner::IOSafety::_input_destination_overlap(
        {canonical_path=>'C:/Work/tree.nwk'},
        {canonical_path=>'c:/work/TREE.NWK',type=>'file'}), 'exact-path', 'case aliases rejected by actual I/O guard');
    ok(SplitAligner::IOSafety::_is_same_or_descendant('C:\\Work\\out\\x', 'c:/work/OUT'), 'Windows descendant alias caught');
    ok(!SplitAligner::IOSafety::_is_same_or_descendant('C:/Work/output/x', 'C:/Work/out'), 'prefix sibling remains distinct');
}
{
    # Simulate only the OS response boundary to exercise fail-closed decisions.
    # This is not a native Win32 API or filesystem acceptance test.
    no warnings qw(redefine once);
    local *SplitAligner::Platform::native_windows = sub { 1 };
    local *SplitAligner::Platform::assert_native_runtime = sub { 1 };
    local $INC{'Win32API/File.pm'} = __FILE__;
    local $INC{'Win32.pm'} = __FILE__;
    my ($drive_type,$fs,$mapping,$attr,$error)=(3,'NTFS','\\Device\\HarddiskVolume3',16,0);
    local *Win32API::File::GetDriveType = sub { $drive_type };
    local *Win32API::File::GetVolumeInformation = sub { $_[6]=$fs; return 1 };
    local *Win32API::File::QueryDosDevice = sub { $_[1]=$mapping."\0\0"; return 1 };
    local *Win32API::File::GetFileAttributes = sub { $attr };
    local *Win32::GetLastError = sub { $error };
    local *Win32API::File::fileLastError = sub { die 'cached Win32API error must not be used' };
    ok(assert_windows_path('C:/work/file','C:/work'), 'simulated local NTFS long path accepted');
    $drive_type=4;
    ok(!eval { assert_windows_path('C:/work/file','C:/work'); 1 }, 'mapped remote drive rejected');
    $drive_type=3; $fs='ReFS';
    ok(!eval { assert_windows_path('C:/work/file','C:/work'); 1 }, 'unsupported filesystem rejected');
    $fs='NTFS'; $mapping='\\??\\D:\\work';
    ok(!eval { assert_windows_path('C:/work/file','C:/work'); 1 }, 'SUBST drive rejected');
    $mapping='\\Device\\HarddiskVolume3'; $attr=0x400;
    ok(!eval { assert_windows_path('C:/work/file','C:/work'); 1 }, 'reparse point rejected');
    $attr=4294967295; $error=5;
    ok(!eval { assert_windows_path('C:/work/file','C:/work'); 1 }, 'attribute access failure rejected');
    $error=2;
    ok(assert_windows_path('C:/work/new-file','C:/work'), 'not-yet-created normal path permitted');
    ok(!eval { assert_windows_stat([123,0],'zero fixture'); 1 }, 'zero file identity rejected');
    ok(assert_windows_stat([123,4294967301],'64-bit fixture'), 'complete nonzero 64-bit identity retained');
}
{
    # Directory entry names belong to their directory, not to process CWD.
    # A filesystem check on the bare basename could inspect an unrelated
    # CWD junction. Require all such checks to receive the full candidate.
    my $dir = tempdir(CLEANUP=>1);
    open my $fh, '>:raw', "$dir/child.txt" or die $!;
    print {$fh} 'entry'; close $fh or die $!;
    no warnings 'redefine';
    local *SplitAligner::TextIO::assert_windows_path = sub {
        die 'must not inspect a readdir basename relative to CWD'
            unless File::Spec->file_name_is_absolute($_[0]);
    };
    my @entries;
    my $ok = eval { @entries=SplitAligner::OutputOwnership::_directory_entries($dir); 1 };
    ok($ok, 'directory basename decoding does not inspect the CWD namespace');
    is(scalar(@entries), 1, 'directory inventory still records its actual child');
}
SKIP: {
    skip 'native Windows runtime/file identities require Windows host', 5 unless native_windows();
    ok(assert_native_runtime(), 'native 64-bit Perl runtime guard');
    my $dir = tempdir(CLEANUP=>1);
    my $path = "$dir/id.txt";
    open my $fh, '>:raw', $path or die $!;
    print {$fh} 'identity';
    close $fh or die $!;
    open $fh, '<:raw', $path or die $!;
    my @open = stat($fh);
    my @path = stat($path);
    ok(assert_windows_stat(\@open, 'test handle'), 'nonzero complete file identity');
    is_deeply([@open[0,1]], [@path[0,1]], 'stat and fstat identify same file');
    close $fh;
    my $alias = "$dir/hardlink.txt";
    ok(link($path,$alias), 'native filesystem can create a hard-link fixture');
    my @alias = stat($alias);
    is_deeply([@alias[0,1]], [@path[0,1]], 'hard-link identity is retained');
}
done_testing();
