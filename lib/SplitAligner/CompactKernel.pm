package SplitAligner::CompactKernel;

use 5.010;
use strict;
use warnings;
use Exporter qw(import);
use Digest::SHA qw(sha256_hex);
use JSON::PP;
use POSIX qw(isfinite);
use Scalar::Util qw(refaddr);
use SplitAligner::Bitset qw(zero_bits ids_to_bits bits_to_ids bits_to_hex
    hex_to_bits and_bits or_bits xor_bits min_bit_id canonical_bits);
use SplitAligner::Newick qw(canonical_split split_taxon_sets branch_id_num
    is_numeric_value sum_component_lengths format_pooled_length);

our @EXPORT_OK = qw(prepare_reference compute_gene compute_gene_materialized materialize_gene
    construct_cell compare_patterns);
my $SCHEMA = 'splitaligner-compact-gene-v1';

sub compare_patterns {
    my ($a, $b) = @_;
    my @a = split(/\|/, $a);
    my @b = split(/\|/, $b);
    for my $i (0 .. (@a < @b ? $#a : $#b)) {
        my $order = branch_id_num($a[$i]) <=> branch_id_num($b[$i]);
        return $order if $order;
    }
    return @a <=> @b;
}

sub prepare_reference {
    my ($axis, $taxa) = @_;
    die "[ERROR][Compact] Reference requires a nonempty axis and at least two taxa.\n"
        unless ref($axis) eq 'ARRAY' && @{$axis}
            && ref($taxa) eq 'ARRAY' && @{$taxa} >= 2;
    my @taxa = sort @{$taxa};
    my (%taxon_ids, %axis_index, %seen_key, %terminal, @rows);
    for my $i (0 .. $#taxa) {
        die "[ERROR][Compact] Invalid or duplicate reference taxon.\n"
            unless defined($taxa[$i]) && !ref($taxa[$i]) && length($taxa[$i])
                && !exists $taxon_ids{$taxa[$i]};
        $taxon_ids{$taxa[$i]} = $i;
    }
    my $n = @taxa;
    my $full = ids_to_bits($n, [0 .. $n - 1]);
    for my $i (0 .. $#{$axis}) {
        my $entry = $axis->[$i];
        die "[ERROR][Compact] Invalid reference coordinate.\n"
            unless ref($entry) eq 'HASH' && defined($entry->{branch_id})
                && $entry->{branch_id} =~ /\AB[0-9]+\z/
                && !exists $axis_index{$entry->{branch_id}} && defined($entry->{split});
        my ($left, $right) = split_taxon_sets($entry->{split});
        die "[ERROR][Compact] Empty reference split side.\n" unless @{$left} && @{$right};
        my %members;
        for my $tip (@{$left}, @{$right}) {
            die "[ERROR][Compact] Reference split does not partition the taxon dictionary.\n"
                unless exists($taxon_ids{$tip}) && !$members{$tip}++;
        }
        die "[ERROR][Compact] Reference split omits taxa.\n" unless keys(%members) == $n;
        die "[ERROR][Compact] Reference split is not canonically encoded.\n"
            unless canonical_split($left, $right) eq $entry->{split};
        my $bits = ids_to_bits($n, [map { $taxon_ids{$_} } @{$left}]);
        my $key = canonical_bits($bits, $full, 0);
        die "[ERROR][Compact] Repeated reference split.\n" if $seen_key{$key}++;
        $terminal{$left->[0]} = 1 if @{$left} == 1;
        $terminal{$right->[0]} = 1 if @{$right} == 1;
        $axis_index{$entry->{branch_id}} = $i;
        push @rows, { branch_id => $entry->{branch_id}, bits => $bits,
            terminal => (@{$left} == 1 || @{$right} == 1) ? 1 : 0 };
    }
    die "[ERROR][Compact] Incomplete reference axis: missing terminal split.\n"
        if grep { !$terminal{$_} } @taxa;
    my $reference_id = sha256_hex(JSON::PP->new->canonical(1)->utf8(1)->encode({
        taxa => \@taxa, axis => [map { {branch_id => $_->{branch_id}, split => $_->{split}} } @{$axis}],
    }));
    return { axis => $axis, axis_ids => [map { $_->{branch_id} } @{$axis}],
        axis_index => \%axis_index, taxa => \@taxa, taxon_ids => \%taxon_ids,
        reference_id => $reference_id, n => $n, byte_length => length($full),
        full => $full, zero => zero_bits($n), rows => \@rows };
}

sub _nodes {
    my ($root) = @_;
    my (@nodes, @stack, %seen);
    push @stack, $root;
    while (@stack) {
        my $node = pop @stack;
        die "[ERROR][Compact] Invalid, cyclic, or shared tree node.\n"
            unless ref($node) eq 'HASH' && ref($node->{children}) eq 'ARRAY'
                && !$seen{refaddr($node)}++;
        push @nodes, $node;
        push @stack, reverse @{$node->{children}};
    }
    return \@nodes;
}

sub _reference_components {
    my ($ctx, $tg, $minimum, $lookup) = @_;
    my %group;
    for my $row (@{$ctx->{rows}}) {
        my $side = and_bits($row->{bits}, $tg);
        next if $side eq $ctx->{zero} || $side eq $tg;
        my $key = bits_to_hex(canonical_bits($side, $tg, $minimum), $ctx->{n});
        push @{$group{$key}}, $row;
    }
    my @components;
    for my $key (sort keys %group) {
        my @rows = sort { branch_id_num($a->{branch_id}) <=> branch_id_num($b->{branch_id}) }
            @{$group{$key}};
        if (@rows == 1 && !$rows[0]{terminal}) {
            my $bits = hex_to_bits($key, $ctx->{n});
            my $left_count = unpack('%32b*', $bits);
            my $right_count = unpack('%32b*', $tg) - $left_count;
            die "[ERROR][Compact] Invariant failure: isolated endpoint-collapsed internal branch.\n"
                if $left_count == 1 || $right_count == 1;
        }
        push @components, { members => join('|', map { $_->{branch_id} } @rows),
            key => $key, observation_id => $lookup->{$key} };
    }
    return [sort { compare_patterns($a->{members}, $b->{members}) } @components];
}

sub compute_gene {
    my ($record) = compute_gene_materialized(@_);
    return $record;
}

# Explicit paired API: callers can reuse the validated, freshly computed
# scientific rows without rerunning reference projection or decoding JSON.
# Keep compute_gene's historical single-record return in every context.
sub compute_gene_materialized {
    my ($root, $ctx) = @_;
    my $nodes = _nodes($root);
    my (%node_index, %tips, @bits);
    $node_index{refaddr($nodes->[$_])} = $_ for 0 .. $#{$nodes};
    for my $i (reverse 0 .. $#{$nodes}) {
        my $node = $nodes->[$i];
        my $side = $ctx->{zero};
        if (@{$node->{children}}) {
            $side = or_bits($side, $bits[$node_index{refaddr($_)}]) for @{$node->{children}};
        } else {
            my $label = $node->{label};
            die "[ERROR][Compact] Gene tree has an undefined terminal taxon.\n"
                unless defined($label);
            die "[ERROR][Compact] Unknown taxon '$label' is absent from the species tree.\n"
                unless exists($ctx->{taxon_ids}{$label});
            die "[ERROR][Compact] Duplicate terminal taxon '$label'.\n" if $tips{$label}++;
            vec($side, $ctx->{taxon_ids}{$label}, 1) = 1;
        }
        $bits[$i] = $side;
    }
    die "[ERROR][Compact] Gene tree must retain at least two taxa.\n" unless keys(%tips) >= 2;
    my $tg = $bits[0];
    my $minimum = min_bit_id($tg, $ctx->{n});
    my %pooled;
    for my $i (1 .. $#{$nodes}) {
        my $key = bits_to_hex(canonical_bits($bits[$i], $tg, $minimum), $ctx->{n});
        push @{$pooled{$key}{source_nodes}}, $i;
        push @{$pooled{$key}{values}}, $nodes->[$i]{length_raw};
    }
    my (@observations, %lookup);
    for my $key (sort keys %pooled) {
        my $pool = $pooled{$key};
        my @values = @{$pool->{values}};
        my $length = 'NA';
        if (!grep { !is_numeric_value($_) } @values) {
            $length = @values == 1 ? $values[0]
                : format_pooled_length(sum_component_lengths(@values));
        }
        push @observations, { key => $key, length => "$length",
            source_nodes => $pool->{source_nodes} };
        $lookup{$key} = scalar(@observations);
    }
    my $record = { schema => $SCHEMA, reference_id => $ctx->{reference_id},
        tg => bits_to_hex($tg, $ctx->{n}), source_node_count => $#{$nodes},
        observations => \@observations,
        components => _reference_components($ctx, $tg, $minimum, \%lookup) };
    my $materialized = _materialize($record, $ctx, 0);
    $record->{$_} = $materialized->{$_} for qw(truth state primitive);
    return ($record, $materialized);
}

sub _same_json {
    return JSON::PP->new->canonical(1)->encode($_[0]) eq JSON::PP->new->canonical(1)->encode($_[1]);
}

sub materialize_gene { return _materialize($_[0], $_[1], 1); }

sub _materialize {
    my ($record, $ctx, $check_saved) = @_;
    die "[ERROR][Compact] Unsupported record or reference identity.\n"
        unless ref($record) eq 'HASH' && ($record->{schema} // '') eq $SCHEMA
            && ($record->{reference_id} // '') eq $ctx->{reference_id}
            && ref($record->{observations}) eq 'ARRAY' && ref($record->{components}) eq 'ARRAY';
    my $tg = hex_to_bits($record->{tg}, $ctx->{n});
    die "[ERROR][Compact] Record has fewer than two retained taxa.\n" if unpack('%32b*', $tg) < 2;
    my $minimum = min_bit_id($tg, $ctx->{n});
    my $count = $record->{source_node_count};
    die "[ERROR][Compact] Invalid source-node count.\n"
        unless defined($count) && !ref($count) && $count =~ /\A[0-9]+\z/ && $count >= 2;
    my (%lookup, %sources, $previous);
    for my $i (0 .. $#{$record->{observations}}) {
        my $obs = $record->{observations}[$i];
        die "[ERROR][Compact] Invalid observation.\n"
            unless ref($obs) eq 'HASH' && ref($obs->{source_nodes}) eq 'ARRAY' && @{$obs->{source_nodes}};
        my $bits = hex_to_bits($obs->{key}, $ctx->{n});
        die "[ERROR][Compact] Noncanonical or unordered observation key.\n"
            if canonical_bits($bits, $tg, $minimum) ne $bits
                || (defined($previous) && $previous ge $obs->{key});
        $previous = $obs->{key};
        my $last = 0;
        for my $node (@{$obs->{source_nodes}}) {
            die "[ERROR][Compact] Invalid, repeated, or unordered source node.\n"
                unless defined($node) && !ref($node) && $node =~ /\A[0-9]+\z/
                    && $node > $last && $node <= $count && !$sources{$node}++;
            $last = $node;
        }
        die "[ERROR][Compact] Invalid observation length.\n"
            unless defined($obs->{length}) && !ref($obs->{length});
        if ($obs->{length} ne 'NA') {
            # Check a copy: numeric coercion must not change a saved raw token's
            # scalar flags and make JSON::PP serialize it as a JSON number.
            my $copy = $obs->{length};
            die "[ERROR][Compact] Invalid observation length.\n" unless is_numeric_value($copy);
            my $number = unpack('d', pack('d', $copy));
            my ($significand) = split(/[eE]/, $copy);
            die "[ERROR][Compact] Observation length is outside binary64 range.\n"
                unless isfinite($number) && !($number == 0 && $significand =~ /[1-9]/);
        }
        $lookup{$obs->{key}} = $i + 1;
    }
    die "[ERROR][Compact] Missing source nodes in observation ledger.\n" unless keys(%sources) == $count;
    my $expected = $check_saved ? _reference_components($ctx, $tg, $minimum, \%lookup)
        : $record->{components};
    die "[ERROR][Compact] Stored reference components differ from reference restriction.\n"
        if $check_saved && !_same_json($record->{components}, $expected);
    my @truth = map { 'NA_struct' } @{$ctx->{axis_ids}};
    my @state = map { 'S' } @{$ctx->{axis_ids}};
    my @primitive = map { 'NA' } @{$ctx->{axis_ids}};
    my (%owners, @groups, @patterns);
    for my $component (@{$expected}) {
        my @members = split(/\|/, $component->{members});
        my $recovered = defined($component->{observation_id}) ? 1 : 0;
        my $length = $recovered ? $record->{observations}[$component->{observation_id} - 1]{length} : 'NA';
        my $owner = { id => $component->{members}, size => scalar(@members),
            recovered => $recovered, length => $length };
        for my $member (@members) {
            my $i = $ctx->{axis_index}{$member};
            $truth[$i] = @members > 1 ? 'NA_fuse' : $recovered ? 'Mapped' : 'not_displayed';
            $state[$i] = $recovered ? (@members > 1 ? 'F' : 'D') : 'U';
            $primitive[$i] = $length if @members == 1 && $recovered;
            $owners{$member} = $owner;
        }
        if (@members > 1) {
            my $side = hex_to_bits($component->{key}, $ctx->{n});
            my $other = xor_bits($side, $tg);
            push @groups, { members => $component->{members},
                split => canonical_split([map { $ctx->{taxa}[$_] } @{bits_to_ids($side, $ctx->{n})}],
                    [map { $ctx->{taxa}[$_] } @{bits_to_ids($other, $ctx->{n})}]),
                status => $recovered ? 'recovered' : 'not_displayed', length => $length };
        }
        push @patterns, $component->{members} if $recovered;
    }
    my $result = { truth => \@truth, state => \@state, primitive => \@primitive,
        groups => [sort { $a->{split} cmp $b->{split} } @groups],
        owners => \%owners, patterns => [sort @patterns] };
    if ($check_saved) {
        for my $field (qw(truth state primitive)) {
            die "[ERROR][Compact] Stored $field differs from component evidence.\n"
                unless ref($record->{$field}) eq 'ARRAY' && _same_json($record->{$field}, $result->{$field});
        }
    }
    return $result;
}

sub construct_cell {
    my ($target, $owner, $view) = @_;
    die "[ERROR][Compact] Fusion view must be whole or direct.\n"
        unless defined($view) && ($view eq 'whole' || $view eq 'direct');
    die "[ERROR][Compact] Invalid fusion target or component owners.\n"
        unless ref($target) eq 'ARRAY' && @{$target} && ref($owner) eq 'HASH';
    my %target_seen;
    for my $branch (@{$target}) {
        die "[ERROR][Compact] Invalid or repeated fusion target member.\n"
            unless defined($branch) && !ref($branch) && $branch =~ /\AB[0-9]+\z/
                && !$target_seen{$branch}++;
    }
    my (@parts, @deleted, %intersection_size);
    for my $branch (@{$target}) {
        my $c = $owner->{$branch};
        if (!defined $c) { push @deleted, $branch; next; }
        push @parts, $c unless $intersection_size{$c->{id}}++;
    }
    my @cut = map { $_->{id} } grep { $intersection_size{$_->{id}} != $_->{size} } @parts;
    my $ids = join(',', map { $_->{id} } @parts);
    if (@deleted || @cut) {
        return ('unavailable', 'not_applicable', 'not_applicable', 'none',
            $ids, join(',', @deleted), join(',', @cut), 'NA');
    }
    my $active = @parts == 1;
    my $reference = $active ? 'active' : 'decomposable';
    my $mode = $active ? 'single_edge_lookup' : 'selected_component_edge_sum';
    return ($reference, 'not_recovered', 'not_recovered', $mode, $ids, '', '', 'NA')
        if grep { !$_->{recovered} } @parts;
    return ($reference, 'recovered', 'direct_only', $mode, $ids, '', '', 'NA')
        if !$active && $view eq 'direct';
    return ($reference, 'recovered', 'missing_length', $mode, $ids, '', '', 'NA')
        if grep { !is_numeric_value($_->{length}) || !isfinite(0 + $_->{length}) } @parts;
    my $value = $parts[0]{length};
    if (!$active) {
        $value = sum_component_lengths(map { $_->{length} } @parts);
        return ($reference, 'recovered', 'nonfinite_sum', $mode, $ids, '', '', 'NA') unless isfinite($value);
        $value = format_pooled_length($value);
    }
    return ($reference, 'recovered', 'available', $mode, $ids, '', '', $value);
}

1;
