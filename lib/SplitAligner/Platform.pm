package SplitAligner::Platform;

use strict;
use warnings;
use Config ();
use Exporter qw(import);

our @EXPORT_OK = qw(
    native_windows assert_native_runtime windows_absolute_path
    assert_windows_path assert_windows_stat path_identity_key
    protect_private_input
);

sub native_windows { return $^O eq 'MSWin32'; }

sub assert_native_runtime {
    return 1 unless native_windows();
    # Perl 5.40's Win32 stat/fstat use the volume serial and the complete
    # 64-bit file index. Older CRT stat implementations are not sufficient
    # for the existing alias and replacement guards.
    die "[ERROR][Windows] Native Windows requires 64-bit Perl 5.40 or newer; use WSL otherwise.\n"
        unless $] >= 5.040 && $Config::Config{ivsize} >= 8 && $Config::Config{ptrsize} >= 8;
    return 1;
}

# Deliberately bounded first native port: local drive paths and ASCII path
# characters. UTF-8 tree labels and gene IDs inside files remain supported.
# Do not silently reinterpret ANSI paths, UNC, device names, ADS, or extended
# path prefixes. This pure function is exercised on all development platforms.
sub windows_absolute_path {
    my ($path, $cwd) = @_;
    die "[ERROR][Windows] Empty filesystem path.\n"
        unless defined $path && length $path;
    die "[ERROR][Windows] Use an ASCII filesystem path (UTF-8 file contents are supported).\n"
        if $path =~ /[^\x20-\x7e]/;
    $path =~ tr{\\}{/};
    die "[ERROR][Windows] UNC, device, and root-relative paths are not supported; use a local drive path.\n"
        if $path =~ m{\A/};
    die "[ERROR][Windows] Drive-relative paths are ambiguous; use C:/path instead.\n"
        if $path =~ /\A[A-Za-z]:(?!\/)/;
    unless ($path =~ m{\A[A-Za-z]:/}) {
        die "[ERROR][Windows] An absolute local-drive working directory is required.\n"
            unless defined $cwd && $cwd =~ m{\A[A-Za-z]:[\\/]};
        $cwd =~ tr{\\}{/};
        $path = "$cwd/$path";
    }
    die "[ERROR][Windows] Use an ASCII working directory.\n"
        if $path =~ /[^\x20-\x7e]/;
    my ($drive, $tail) = $path =~ m{\A([A-Za-z]):/(.*)\z};
    die "[ERROR][Windows] Invalid drive path.\n" unless defined $drive;
    my @part;
    for my $part (split m{/+}, $tail, -1) {
        next if $part eq '' || $part eq '.';
        if ($part eq '..') { pop @part if @part; next; }
        die "[ERROR][Windows] Unsupported path component '$part' (reserved character or trailing dot/space).\n"
            if $part =~ /[<>:"|?*]/ || $part =~ /[. ]\z/;
        die "[ERROR][Windows] Reserved device name in path: $part.\n"
            if $part =~ /\A(?:CON|PRN|AUX|NUL|CLOCK\$|COM[1-9]|LPT[1-9])(?:\.|\z)/i;
        die "[ERROR][Windows] Short-name aliases are not supported; use long path components.\n"
            if $part =~ /~/;
        push @part, $part;
    }
    my $result = uc($drive) . ':/' . join('/', @part);
    die "[ERROR][Windows] Path is too long for the bounded native port; use a shorter work directory.\n"
        if length($result) >= 240;
    return $result;
}

sub assert_windows_path {
    my ($path, $cwd) = @_;
    return 1 unless native_windows();
    assert_native_runtime();
    require Cwd;
    $cwd //= Cwd::getcwd();
    my $absolute = windows_absolute_path($path, $cwd);
    # Reject junctions as well as symlinks. Cwd's lexical spelling alone is
    # insufficient to detect a not-yet-created output under a junction.
    eval { require Win32API::File; 1 }
        or die "[ERROR][Windows] Win32API::File is required for native filesystem guards (included in full Strawberry Perl): $@";
    require Win32;
    my ($drive, $tail) = $absolute =~ m{\A([A-Z]:)/(.*)\z};
    my $drive_root = $drive . '/';
    die "[ERROR][Windows] Native inputs and outputs require a local fixed NTFS volume.\n"
        unless Win32API::File::GetDriveType($drive_root) == 3;
    my ($volume, $serial, $max_component, $flags, $filesystem) = ('', 0, 0, 0, '');
    Win32API::File::GetVolumeInformation($drive_root, $volume, 256,
        $serial, $max_component, $flags, $filesystem, 256)
        or die "[ERROR][Windows] Cannot inspect volume: $drive_root.\n";
    die "[ERROR][Windows] Native inputs and outputs require NTFS (found $filesystem).\n"
        unless $filesystem eq 'NTFS';
    my $target = '';
    Win32API::File::QueryDosDevice($drive, $target, 4096)
        or die "[ERROR][Windows] Cannot inspect drive mapping: $drive.\n";
    ($target) = split /\x00/, $target;
    die "[ERROR][Windows] SUBST, mapped, and nonstandard drive aliases are not supported: $drive.\n"
        unless $target =~ /\A\\Device\\HarddiskVolume[0-9]+\z/i;
    my $prefix = $drive . '/';
    for my $part ('', split(m{/}, $tail)) {
        $prefix .= $part if length $part;
        my $attributes = Win32API::File::GetFileAttributes($prefix);
        # Win32API::File's GetFileAttributes wrapper does not reliably update
        # its cached fileLastError for INVALID_FILE_ATTRIBUTES (0xffffffff).
        # Capture the OS error before another Windows call can replace it.
        my $attribute_error = Win32::GetLastError();
        if ($attributes == -1 || $attributes == 4294967295) {
            my $error = 0 + $attribute_error;
            die "[ERROR][Windows] Cannot inspect filesystem attributes: $prefix (error $error).\n"
                unless $error == 2 || $error == 3;
        } else {
            die "[ERROR][Windows] Reparse points (including junctions/symlinks) are not supported: $prefix.\n"
                if $attributes & 0x400;
        }
        $prefix .= '/' unless $prefix =~ m{/\z};
    }
    return 1;
}

sub assert_windows_stat {
    my ($stat, $context) = @_;
    return 1 unless native_windows();
    assert_native_runtime();
    die "[ERROR][Windows] Filesystem has no usable volume/file identity for $context.\n"
        unless ref($stat) eq 'ARRAY' && @$stat && defined($stat->[0])
            && defined($stat->[1]) && $stat->[1] > 0;
    return 1;
}

sub path_identity_key {
    my ($path) = @_;
    return $path unless native_windows();
    $path =~ tr{\\}{/};
    # Conservatively reject aliases even in case-sensitive Windows folders.
    return lc($path);
}

my @protected;
sub protect_private_input {
    my ($path) = @_;
    chmod 0400, $path or die "[ERROR] Cannot protect private input snapshot: $!\n";
    if (native_windows()) {
        my @stat = stat($path);
        assert_windows_stat(\@stat, 'private input snapshot');
        push @protected, [$path, $stat[0], $stat[1]];
    }
}

# File::Temp's END cleanup runs later. Remove the Windows read-only attribute
# only from the same private snapshot, never from a replaced path.
END {
    for my $item (@protected) {
        my ($path, $device, $inode) = @$item;
        my @stat = lstat($path);
        next unless @stat;
        if (-l $path || $stat[0] != $device || $stat[1] != $inode) {
            warn "[WARNING][Windows] Private input changed; cannot clear read-only attribute: $path\n";
            next;
        }
        chmod 0600, $path
            or warn "[WARNING][Windows] Cannot make private input removable: $path: $!\n";
    }
}

1;
