package SplitAligner::TextIO;

use strict;
use warnings;
use Cwd ();
use Encode qw(decode encode FB_CROAK LEAVE_SRC);
use Exporter qw(import);

our @EXPORT_OK = qw(
    abs_path_utf8
    close_utf8_writer
    configure_utf8_stdio
    decode_argv_utf8
    decode_filesystem_path_utf8
    decode_utf8_strict
    encode_utf8_strict
    getcwd_utf8
    open_utf8_writer
    print_utf8
    read_raw_file
    read_utf8_file
    read_utf8_lines
    write_utf8_file
);

my $STRICT = FB_CROAK | LEAVE_SRC;

sub configure_utf8_stdio {
    binmode STDOUT, ':raw:encoding(UTF-8)'
        or die "[ERROR] Cannot configure STDOUT for UTF-8: $!\n";
    binmode STDERR, ':raw:encoding(UTF-8)'
        or die "[ERROR] Cannot configure STDERR for UTF-8: $!\n";
}

sub decode_argv_utf8 {
    for my $i (0 .. $#ARGV) {
        next if utf8::is_utf8($ARGV[$i]);
        $ARGV[$i] = decode_utf8_strict($ARGV[$i], "command-line argument " . ($i + 1));
    }
}

sub decode_filesystem_path_utf8 {
    my ($path, $context) = @_;
    return undef unless defined $path;
    $context //= 'filesystem path';
    return decode_utf8_strict($path, $context);
}

sub getcwd_utf8 {
    my $path = Cwd::getcwd();
    die "[ERROR] Cannot obtain current working directory: $!\n"
        unless defined $path;
    return decode_filesystem_path_utf8($path, 'current working directory');
}

sub abs_path_utf8 {
    my ($path, $context) = @_;
    $context //= 'absolute filesystem path';
    my $absolute = Cwd::abs_path($path);
    return undef unless defined $absolute;
    return decode_filesystem_path_utf8($absolute, $context);
}

sub decode_utf8_strict {
    my ($bytes, $context) = @_;
    $context //= 'UTF-8 text';
    return $bytes if utf8::is_utf8($bytes);

    my $text = eval { decode('UTF-8', $bytes, $STRICT) };
    die "[ERROR][$context] Invalid UTF-8: $@" if $@;
    return $text;
}

sub encode_utf8_strict {
    my ($text, $context) = @_;
    $context //= 'UTF-8 text';
    die "[ERROR][$context] Undefined text cannot be encoded.\n" unless defined $text;

    my $bytes = eval { encode('UTF-8', $text, $STRICT) };
    die "[ERROR][$context] Cannot encode UTF-8: $@" if $@;
    return $bytes;
}

sub read_raw_file {
    my ($path) = @_;
    open my $fh, '<:raw', $path or die "[ERROR] Cannot open $path: $!\n";
    local $/;
    my $bytes = <$fh>;
    close $fh or die "[ERROR] Cannot close $path: $!\n";
    return defined $bytes ? $bytes : '';
}

sub read_utf8_file {
    my ($path) = @_;
    my $bytes = read_raw_file($path);
    die "[ERROR][$path] UTF-8 BOM is not permitted.\n"
        if substr($bytes, 0, 3) eq "\xEF\xBB\xBF";
    return decode_utf8_strict($bytes, $path);
}

sub read_utf8_lines {
    my ($path) = @_;
    my $text = read_utf8_file($path);
    my @lines = split(/\n/, $text, -1);
    s/\r\z// for @lines;
    pop @lines if @lines && $lines[-1] eq '';
    return \@lines;
}

sub open_utf8_writer {
    my ($path) = @_;
    open my $fh, '>:raw', $path or die "[ERROR] Cannot write $path: $!\n";
    return $fh;
}

sub print_utf8 {
    my ($fh, @parts) = @_;
    my $text = join('', map { defined $_ ? $_ : '' } @parts);
    my $bytes = encode_utf8_strict($text, 'output text');
    print {$fh} $bytes or die "[ERROR] Cannot write UTF-8 output: $!\n";
}

sub close_utf8_writer {
    my ($fh, $path) = @_;
    $path //= 'UTF-8 output';
    close $fh or die "[ERROR] Cannot close $path: $!\n";
}

sub write_utf8_file {
    my ($path, $text) = @_;
    my $fh = open_utf8_writer($path);
    print_utf8($fh, $text);
    close_utf8_writer($fh, $path);
}

1;
