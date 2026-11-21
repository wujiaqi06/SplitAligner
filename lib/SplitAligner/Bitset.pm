package SplitAligner::Bitset;

use 5.010;
use strict;
use warnings;
use Exporter qw(import);
use utf8 ();

our @EXPORT_OK = qw(zero_bits ids_to_bits bits_to_ids bits_to_hex hex_to_bits
    validate_bits and_bits or_bits xor_bits min_bit_id canonical_bits);

sub _count {
    my ($n) = @_;
    die "[ERROR][Bitset] Invalid bit count.\n"
        unless defined($n) && !ref($n) && $n =~ /\A[0-9]+\z/ && $n > 0;
    return $n;
}

sub zero_bits { return "\0" x int((_count($_[0]) + 7) / 8); }

sub validate_bits {
    my ($bits, $n) = @_;
    _count($n);
    die "[ERROR][Bitset] Expected a fixed-length byte string.\n"
        if !defined($bits) || ref($bits) || utf8::is_utf8($bits)
            || length($bits) != int(($n + 7) / 8);
    if ($n % 8) {
        my $unused = 255 ^ ((1 << ($n % 8)) - 1);
        die "[ERROR][Bitset] Nonzero unused tail bits.\n"
            if ord(substr($bits, -1, 1)) & $unused;
    }
    return 1;
}

sub ids_to_bits {
    my ($n, $ids) = @_;
    my $bits = zero_bits($n);
    die "[ERROR][Bitset] Expected an ID array.\n" unless ref($ids) eq 'ARRAY';
    my %seen;
    for my $id (@{$ids}) {
        die "[ERROR][Bitset] Invalid or duplicate taxon ID.\n"
            unless defined($id) && !ref($id) && $id =~ /\A[0-9]+\z/
                && $id < $n && !$seen{$id}++;
        vec($bits, $id, 1) = 1;
    }
    return $bits;
}

sub bits_to_ids {
    my ($bits, $n) = @_;
    validate_bits($bits, $n);
    return [ grep { vec($bits, $_, 1) } 0 .. $n - 1 ];
}

sub bits_to_hex {
    my ($bits, $n) = @_;
    validate_bits($bits, $n);
    return unpack('H*', $bits);
}

sub hex_to_bits {
    my ($hex, $n) = @_;
    _count($n);
    die "[ERROR][Bitset] Invalid canonical hexadecimal bitset.\n"
        unless defined($hex) && !ref($hex)
            && length($hex) == 2 * int(($n + 7) / 8)
            && $hex =~ /\A[0-9a-f]+\z/;
    my $bits = pack('H*', $hex);
    validate_bits($bits, $n);
    return $bits;
}

sub _pair {
    my ($left, $right) = @_;
    die "[ERROR][Bitset] Operands must be equal-length byte strings.\n"
        if !defined($left) || !defined($right) || ref($left) || ref($right)
            || utf8::is_utf8($left) || utf8::is_utf8($right)
            || !length($left) || length($left) != length($right);
}

# This module deliberately does not enable the newer `bitwise` feature. In
# legacy Core Perl, interpolation forces string operations even if an operand
# acquired a numeric slot elsewhere. Keep these expressions inside this module.
sub and_bits { _pair(@_); my ($a, $b) = @_; return "$a" & "$b"; }
sub or_bits  { _pair(@_); my ($a, $b) = @_; return "$a" | "$b"; }
sub xor_bits { _pair(@_); my ($a, $b) = @_; return "$a" ^ "$b"; }

sub min_bit_id {
    my ($bits, $n) = @_;
    validate_bits($bits, $n);
    for my $byte (0 .. length($bits) - 1) {
        my $value = ord(substr($bits, $byte, 1));
        next unless $value;
        for my $bit (0 .. 7) { return 8 * $byte + $bit if $value & (1 << $bit); }
    }
    return undef;
}

sub canonical_bits {
    my ($side, $taxa, $minimum) = @_;
    _pair($side, $taxa);
    die "[ERROR][Bitset] Invalid canonical anchor.\n"
        unless defined($minimum) && !ref($minimum) && $minimum =~ /\A[0-9]+\z/
            && $minimum < 8 * length($taxa) && vec($taxa, $minimum, 1);
    die "[ERROR][Bitset] Split side is not a subset of its taxa.\n"
        unless and_bits($side, $taxa) eq $side;
    return vec($side, $minimum, 1) ? xor_bits($side, $taxa) : "$side";
}

1;
