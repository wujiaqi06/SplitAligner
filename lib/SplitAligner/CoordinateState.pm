package SplitAligner::CoordinateState;

use strict;
use warnings;
use Exporter qw(import);
use SplitAligner::Newick qw(is_numeric_value);
use SplitAligner::TextIO qw(read_utf8_lines write_utf8_file);

our @EXPORT_OK = qw(
    read_matrix_table
    read_state_matrix
    state_code_legend
    state_schema_version
    validate_state_for_matrix
    write_state_matrix
);

my $STATE_SCHEMA = 'SplitAligner-primitive-state-v1';
my %VALID_STATE = map { $_ => 1 } qw(S D F U);

sub state_schema_version {
    return $STATE_SCHEMA;
}

sub state_code_legend {
    return {
        S => 'STRUCT_ABSENT: at least one projected split side is empty',
        D => 'DIRECT_MAPPED: projected primitive split is present directly',
        F => 'FUSED_MAPPED: primitive belongs to a mapped composite coordinate',
        U => 'PROJECTED_UNMAPPED: projected identity remains but is not mapped',
    };
}

sub write_state_matrix {
    my ($path, $axis, $rows) = @_;
    _validate_axis_ids($axis, 'state output axis');
    die "[ERROR] State rows must be a non-empty array.\n"
        unless ref($rows) eq 'ARRAY' && @{$rows};

    my @text = (join("\t", 'gene', @{$axis}));
    my %seen_gene;
    for my $row (@{$rows}) {
        die "[ERROR] Invalid state output row.\n" unless ref($row) eq 'HASH';
        my $gene = $row->{gene_id};
        my $states = $row->{states};
        die "[ERROR] Empty gene ID in state output.\n" unless defined $gene && $gene ne '';
        die "[ERROR] Duplicate gene ID '$gene' in state output.\n" if $seen_gene{$gene}++;
        die "[ERROR] State width mismatch for gene '$gene'.\n"
            unless ref($states) eq 'ARRAY' && @{$states} == @{$axis};
        for my $state (@{$states}) {
            die "[ERROR] Unknown coordinate state '$state' for gene '$gene'.\n"
                unless defined $state && $VALID_STATE{$state};
        }
        push @text, join("\t", $gene, @{$states});
    }
    write_utf8_file($path, join("\n", @text) . "\n");
}

sub read_state_matrix {
    my ($path) = @_;
    my $lines = read_utf8_lines($path);
    die "[ERROR] Coordinate-state sidecar is empty: $path\n" unless @{$lines};

    my @header = split(/\t/, $lines->[0], -1);
    die "[ERROR] Coordinate-state header must begin with 'gene': $path\n"
        unless @header >= 2 && shift(@header) eq 'gene';
    _validate_axis_ids(\@header, "coordinate-state header in $path");

    my (@genes, %rows, %seen_gene);
    for my $index (1 .. $#{$lines}) {
        my $line_no = $index + 1;
        my $line = $lines->[$index];
        next if $line eq '';
        my @field = split(/\t/, $line, -1);
        die "[ERROR][$path line $line_no] Expected " . (1 + @header)
            . " tab-separated fields, found " . scalar(@field) . ".\n"
            unless @field == 1 + @header;
        my $gene = shift @field;
        die "[ERROR][$path line $line_no] Empty gene ID.\n" if $gene eq '';
        die "[ERROR][$path line $line_no] Duplicate gene ID '$gene'.\n"
            if $seen_gene{$gene}++;
        for my $state (@field) {
            die "[ERROR][$path line $line_no] Unknown coordinate state '$state'.\n"
                unless $VALID_STATE{$state};
        }
        push @genes, $gene;
        $rows{$gene} = \@field;
    }
    die "[ERROR] Coordinate-state sidecar has no gene rows: $path\n" unless @genes;

    return {
        path       => $path,
        schema     => $STATE_SCHEMA,
        axis       => \@header,
        genes      => \@genes,
        rows       => \%rows,
        gene_count => scalar(@genes),
        cell_count => scalar(@genes) * scalar(@header),
    };
}

sub read_matrix_table {
    my ($path) = @_;
    my $lines = read_utf8_lines($path);
    die "[ERROR] Matrix is empty: $path\n" unless @{$lines};

    my @header = split(/\t/, $lines->[0], -1);
    die "[ERROR] Matrix header must begin with 'gene': $path\n"
        unless @header >= 2 && shift(@header) eq 'gene';
    my %seen_column;
    for my $column (@header) {
        die "[ERROR] Empty matrix column in $path.\n" if $column eq '';
        die "[ERROR] Duplicate matrix column '$column' in $path.\n" if $seen_column{$column}++;
        die "[ERROR] Invalid matrix branch column '$column' in $path.\n"
            unless $column =~ /^B\d+(?:\|B\d+)*$/;
    }

    my (@genes, %rows, %seen_gene);
    for my $index (1 .. $#{$lines}) {
        my $line_no = $index + 1;
        my $line = $lines->[$index];
        next if $line eq '';
        my @field = split(/\t/, $line, -1);
        die "[ERROR][$path line $line_no] Expected " . (1 + @header)
            . " tab-separated fields, found " . scalar(@field) . ".\n"
            unless @field == 1 + @header;
        my $gene = shift @field;
        die "[ERROR][$path line $line_no] Empty gene ID.\n" if $gene eq '';
        die "[ERROR][$path line $line_no] Duplicate gene ID '$gene'.\n"
            if $seen_gene{$gene}++;
        push @genes, $gene;
        $rows{$gene} = \@field;
    }
    die "[ERROR] Matrix has no gene rows: $path\n" unless @genes;

    my @primitive = grep { /^B\d+$/ } @header;
    return {
        path       => $path,
        branches   => \@header,
        primitive  => \@primitive,
        genes      => \@genes,
        rows       => \%rows,
        gene_count => scalar(@genes),
    };
}

sub validate_state_for_matrix {
    my (%arg) = @_;
    my $state = $arg{state};
    my $matrix = $arg{matrix};
    my $axis = $arg{axis};
    my $role = $arg{role} // 'matrix';
    my $allow_na_fuse = $arg{allow_na_fuse} ? 1 : 0;
    my $allow_classified = $arg{allow_classified} ? 1 : 0;

    die "[ERROR] Invalid state object for $role.\n" unless ref($state) eq 'HASH';
    die "[ERROR] Invalid matrix object for $role.\n" unless ref($matrix) eq 'HASH';
    _validate_axis_ids($axis, "$role expected axis");

    _assert_array_equal($state->{axis}, $axis, "$role state axis", "$role expected axis");
    _assert_array_equal($matrix->{primitive}, $axis, "$role matrix primitive axis", "$role expected axis");
    _assert_array_equal($state->{genes}, $matrix->{genes}, "$role state gene order", "$role matrix gene order");

    my %matrix_index = map { $matrix->{branches}[$_] => $_ } 0 .. $#{$matrix->{branches}};
    for my $gene (@{$matrix->{genes}}) {
        my $states = $state->{rows}{$gene};
        my $values = $matrix->{rows}{$gene};
        for my $i (0 .. $#{$axis}) {
            my $branch = $axis->[$i];
            my $value = $values->[ $matrix_index{$branch} ];
            my $code = $states->[$i];

            if (is_numeric_value($value)) {
                die "[ERROR][$role gene $gene branch $branch] Numeric primitive value is incompatible with state '$code'.\n"
                    unless $code eq 'D';
                next;
            }
            next if defined $value && $value eq 'NA';
            if ($allow_na_fuse && defined $value && $value eq 'NA_fuse') {
                die "[ERROR][$role gene $gene branch $branch] NA_fuse is incompatible with state '$code'.\n"
                    unless $code eq 'F';
                next;
            }
            if ($allow_classified && defined $value && $value eq 'NA_struct') {
                die "[ERROR][$role gene $gene branch $branch] NA_struct is incompatible with state '$code'.\n"
                    unless $code eq 'S';
                next;
            }
            if ($allow_classified && defined $value && $value eq 'NA_topo') {
                die "[ERROR][$role gene $gene branch $branch] NA_topo is incompatible with state '$code'.\n"
                    unless $code eq 'U';
                next;
            }
            # Other nonnumeric spellings (for example NaN, Inf, or an empty
            # field in a legacy hand-built table) carry no numeric evidence.
            # They are preserved as unavailable values; state remains the
            # authority for structural and mapping classification.
            next;
        }
    }
    return 1;
}

sub _validate_axis_ids {
    my ($axis, $name) = @_;
    die "[ERROR] $name must be a non-empty array.\n"
        unless ref($axis) eq 'ARRAY' && @{$axis};
    my %seen;
    for my $branch (@{$axis}) {
        die "[ERROR] Invalid branch ID '$branch' in $name.\n"
            unless defined $branch && $branch =~ /^B\d+$/;
        die "[ERROR] Duplicate branch ID '$branch' in $name.\n" if $seen{$branch}++;
    }
}

sub _assert_array_equal {
    my ($left, $right, $left_name, $right_name) = @_;
    die "[ERROR] Missing $left_name or $right_name.\n"
        unless ref($left) eq 'ARRAY' && ref($right) eq 'ARRAY';
    my $max = @{$left} > @{$right} ? @{$left} : @{$right};
    for my $i (0 .. $max - 1) {
        my $l = $i < @{$left} ? $left->[$i] : '<missing>';
        my $r = $i < @{$right} ? $right->[$i] : '<missing>';
        next if $l eq $r;
        die "[ERROR] $left_name differs from $right_name at position " . ($i + 1)
            . ": '$l' vs '$r'.\n";
    }
    return 1;
}

1;
