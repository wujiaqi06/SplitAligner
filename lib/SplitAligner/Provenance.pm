package SplitAligner::Provenance;

use strict;
use warnings;
use Digest::SHA qw(sha256_hex);
use Exporter qw(import);
use JSON::PP;
use SplitAligner::Newick qw(canonicalize_split);
use SplitAligner::TextIO qw(
    decode_utf8_strict
    encode_utf8_strict
    read_raw_file
    read_utf8_lines
    write_utf8_file
);

our @EXPORT_OK = qw(
    axis_sha256
    axis_to_text
    file_sha256
    read_axis_ledger
    read_json_file
    verify_axis_equal
    write_axis_ledger
    write_json_file
);

sub axis_to_text {
    my ($axis) = @_;
    _validate_axis($axis, 'primitive axis');
    return join('', map { $_->{branch_id} . "\t" . $_->{split} . "\n" } @{$axis});
}

sub axis_sha256 {
    my ($axis) = @_;
    return sha256_hex(encode_utf8_strict(axis_to_text($axis), 'primitive axis serialization'));
}

sub write_axis_ledger {
    my ($path, $axis) = @_;
    write_utf8_file($path, axis_to_text($axis));
}

sub read_axis_ledger {
    my ($path) = @_;
    my @axis;
    my %seen_branch;
    my %seen_split;
    my $line_no = 0;
    for my $line (@{ read_utf8_lines($path) }) {
        $line_no++;
        next if $line eq '';
        my @field = split(/\t/, $line, -1);
        die "[ERROR][$path line $line_no] Expected exactly two tab-separated fields.\n"
            unless @field == 2;
        my ($branch_id, $split) = @field;
        die "[ERROR][$path line $line_no] Invalid branch ID '$branch_id'.\n"
            unless $branch_id =~ /^B\d+$/;
        die "[ERROR][$path line $line_no] Duplicate branch ID '$branch_id'.\n"
            if $seen_branch{$branch_id}++;
        die "[ERROR][$path line $line_no] Duplicate canonical split '$split'.\n"
            if $seen_split{$split}++;
        my $canonical = eval { canonicalize_split($split) };
        die "[ERROR][$path line $line_no] Invalid canonical split key '$split': $@"
            if $@;
        die "[ERROR][$path line $line_no] Non-canonical split key '$split'; expected '$canonical'.\n"
            unless $canonical eq $split;
        push @axis, { branch_id => $branch_id, split => $split };
    }
    die "[ERROR] Primitive axis ledger is empty: $path\n" unless @axis;
    return \@axis;
}

sub verify_axis_equal {
    my ($left, $right, $left_name, $right_name) = @_;
    $left_name  //= 'left axis';
    $right_name //= 'right axis';
    my $left_text  = axis_to_text($left);
    my $right_text = axis_to_text($right);
    return 1 if $left_text eq $right_text;

    my $max = @{$left} > @{$right} ? @{$left} : @{$right};
    for my $i (0 .. $max - 1) {
        my $l = $left->[$i];
        my $r = $right->[$i];
        my $l_text = $l ? "$l->{branch_id}\t$l->{split}" : '<missing>';
        my $r_text = $r ? "$r->{branch_id}\t$r->{split}" : '<missing>';
        next if $l_text eq $r_text;
        die "[ERROR] Primitive axis mismatch at coordinate " . ($i + 1)
            . ": $left_name='$l_text'; $right_name='$r_text'.\n";
    }
    die "[ERROR] Primitive axes differ: $left_name vs $right_name.\n";
}

sub file_sha256 {
    my ($path) = @_;
    open my $fh, '<:raw', $path or die "[ERROR] Cannot hash $path: $!\n";
    binmode $fh;
    my $sha = Digest::SHA->new(256);
    $sha->addfile($fh);
    close $fh;
    return $sha->hexdigest;
}

sub write_json_file {
    my ($path, $data) = @_;
    my $json = JSON::PP->new->canonical(1)->pretty(1)->utf8(1);
    open my $fh, '>:raw', $path or die "[ERROR] Cannot write $path: $!\n";
    print {$fh} $json->encode($data);
    close $fh or die "[ERROR] Cannot close $path: $!\n";
}

sub read_json_file {
    my ($path) = @_;
    my $raw = read_raw_file($path);
    die "[ERROR][$path] UTF-8 BOM is not permitted.\n"
        if substr($raw, 0, 3) eq "\xEF\xBB\xBF";
    my $text = decode_utf8_strict($raw, $path);
    my $data = eval { JSON::PP->new->utf8(0)->decode($text) };
    die "[ERROR] Invalid JSON manifest $path: $@\n" if $@;
    return $data;
}

sub _validate_axis {
    my ($axis, $name) = @_;
    die "[ERROR] $name must be a non-empty array.\n"
        unless ref($axis) eq 'ARRAY' && @{$axis};

    my (%seen_branch, %seen_split);
    for my $coordinate (@{$axis}) {
        die "[ERROR] Invalid coordinate in $name.\n" unless ref($coordinate) eq 'HASH';
        my $branch_id = $coordinate->{branch_id} // '';
        my $split = $coordinate->{split};
        die "[ERROR] Invalid branch ID '$branch_id' in $name.\n"
            unless $branch_id =~ /^B\d+$/;
        die "[ERROR] Duplicate branch ID '$branch_id' in $name.\n"
            if $seen_branch{$branch_id}++;
        die "[ERROR] Undefined canonical split for '$branch_id' in $name.\n"
            unless defined $split;
        die "[ERROR] Duplicate canonical split '$split' in $name.\n"
            if $seen_split{$split}++;
        my $canonical = eval { canonicalize_split($split) };
        die "[ERROR] Invalid canonical split key '$split' in $name: $@" if $@;
        die "[ERROR] Non-canonical split key '$split' in $name; expected '$canonical'.\n"
            unless $canonical eq $split;
    }
}

1;
