package SplitAligner::Newick;

use strict;
use warnings;
use Exporter qw(import);
use Scalar::Util qw(refaddr);
use SplitAligner::TextIO qw(
    decode_utf8_strict
    encode_utf8_strict
    read_utf8_file
    read_utf8_lines
);

our @EXPORT_OK = qw(
    assign_branch_ids
    branch_id_num
    canonical_split
    canonical_taxon_set_key
    canonicalize_split
    edge_records
    gene_split_records
    is_numeric_value
    parse_gene_tree_file
    parse_labeled_species_tree_file
    parse_newick
    parse_species_tree_file
    retained_primitive_axis
    serialize_figtree
    serialize_branch_subtree
    serialize_for_split
    serialize_plain
    split_taxon_sets
    tip_labels
    validate_gene_taxa
);

my $NUMERIC_RE = qr/[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?/;
my $STRUCTURED_SPLIT_PREFIX = 'HX1:';
my $STRUCTURED_SET_PREFIX = 'HS1:';

sub is_numeric_value {
    my ($value) = @_;
    return 0 unless defined $value;
    return $value =~ /\A$NUMERIC_RE\z/ ? 1 : 0;
}

sub parse_species_tree_file {
    my ($path) = @_;
    my $raw = _slurp($path);
    die "[ERROR][$path] Empty file.\n" if $raw !~ /\S/;
    return parse_newick($raw, source => $path);
}

sub parse_labeled_species_tree_file {
    my ($path) = @_;
    my $raw = _slurp($path);
    die "[ERROR][$path] Empty file.\n" if $raw !~ /\S/;

    my $tree_pos = index($raw, '(');
    die "[ERROR][$path] No '(' found.\n" if $tree_pos < 0;
    my $record_id = substr($raw, 0, $tree_pos);
    $record_id =~ s/^\s+|\s+$//g;
    my $tree = substr($raw, $tree_pos);
    my $root = parse_newick($tree, source => $path, allow_branch_tags => 1);
    return ($record_id, $root);
}

sub parse_gene_tree_file {
    my ($path) = @_;
    my @records;
    my %seen_gene;
    my $line_no = 0;
    for my $line (@{ read_utf8_lines($path) }) {
        $line_no++;
        next if $line =~ /^\s*\z/;

        my $tree_pos = index($line, '(');
        die "[ERROR][$path line $line_no] No '(' found to start a gene tree.\n"
            if $tree_pos < 0;

        my $gene_id = substr($line, 0, $tree_pos);
        $gene_id =~ s/^\s+|\s+$//g;
        die "[ERROR][$path line $line_no] Empty gene ID.\n" if $gene_id eq '';
        die "[ERROR][$path line $line_no gene $gene_id] Duplicate gene ID; first seen on line $seen_gene{$gene_id}.\n"
            if exists $seen_gene{$gene_id};

        my $tree = substr($line, $tree_pos);
        my $root = parse_newick(
            $tree,
            source => "$path line $line_no gene $gene_id",
        );
        $seen_gene{$gene_id} = $line_no;
        push @records, {
            gene_id => $gene_id,
            line_no => $line_no,
            root    => $root,
        };
    }

    die "[ERROR][$path] No non-empty gene-tree records found.\n" unless @records;
    return \@records;
}

sub parse_newick {
    my ($text, %opt) = @_;
    my $source = $opt{source} // 'Newick';

    die "[ERROR][$source] Empty tree string.\n" unless defined $text && $text =~ /\S/;
    die "[ERROR][$source] Quoted labels are not supported.\n" if $text =~ /['\"]/;
    die "[ERROR][$source] Bracket comments are not supported.\n" if $text =~ /[\[\]]/;

    my $state = {
        text              => $text,
        pos               => 0,
        len               => length($text),
        source            => $source,
        allow_branch_tags => $opt{allow_branch_tags} ? 1 : 0,
    };

    _skip_ws($state);
    my $root = _parse_subtree($state);
    _skip_ws($state);
    _expect_char($state, ';');
    _skip_ws($state);
    _error($state, "Extra content after the terminating ';'.")
        if $state->{pos} != $state->{len};

    my @tips = tip_labels($root);
    die "[ERROR][$source] Tree must contain at least two terminal taxa.\n" if @tips < 2;
    my %seen;
    for my $tip (@tips) {
        die "[ERROR][$source] Duplicate terminal taxon '$tip'.\n" if $seen{$tip}++;
    }
    return $root;
}

sub validate_gene_taxa {
    my ($records, $species_root, $source) = @_;
    $source //= 'gene trees';
    my %species = map { $_ => 1 } tip_labels($species_root);

    for my $record (@{$records}) {
        for my $taxon (tip_labels($record->{root})) {
            die "[ERROR][$source line $record->{line_no} gene $record->{gene_id}] Unknown taxon '$taxon' is absent from the species tree.\n"
                unless $species{$taxon};
        }
    }
    return 1;
}

sub assign_branch_ids {
    my ($root) = @_;
    my $counter = 0;

    for my $leaf (_leaf_nodes($root)) {
        $leaf->{branch_id} = 'B' . ++$counter;
    }
    for my $node (_postorder_internal_nodes($root)) {
        next if $node == $root;
        $node->{branch_id} = 'B' . ++$counter;
    }
    delete $root->{branch_id};
    return $counter;
}

sub tip_labels {
    my ($root) = @_;
    return map { $_->{label} } _leaf_nodes($root);
}

sub branch_id_num {
    my ($branch_id) = @_;
    return undef unless defined $branch_id && $branch_id =~ /^B(\d+)$/;
    return 0 + $1;
}

sub canonical_split {
    my ($left_ref, $right_ref) = @_;
    my @left  = sort @{$left_ref};
    my @right = sort @{$right_ref};

    # Preserve the historical readable representation only when every label
    # has an unambiguous boundary under the legacy '..' set separator.
    if (!_requires_structured_split_key(@left, @right)) {
        my $left  = join('..', @left);
        my $right = join('..', @right);
        return $left le $right ? "$left||$right" : "$right||$left";
    }

    # Labels that are not provably legacy-safe use an injective byte encoding.
    # Colons cannot occur in an unquoted taxon label under the supported Newick
    # grammar, so the versioned prefix is disjoint from every legal legacy key.
    # Hex payloads cannot contain ',' or '||'.
    my $left_key  = canonical_taxon_set_key(@left);
    my $right_key = canonical_taxon_set_key(@right);
    return $left_key le $right_key
        ? "$STRUCTURED_SPLIT_PREFIX$left_key||$right_key"
        : "$STRUCTURED_SPLIT_PREFIX$right_key||$left_key";
}

sub canonical_taxon_set_key {
    my @hex = sort map { unpack('H*', _label_octets($_)) } @_;
    return $STRUCTURED_SET_PREFIX . join(',', @hex);
}

sub split_taxon_sets {
    my ($split) = @_;
    die "[ERROR] Undefined canonical split key.\n" unless defined $split;

    if (index($split, $STRUCTURED_SPLIT_PREFIX) == 0) {
        my $payload = substr($split, length($STRUCTURED_SPLIT_PREFIX));
        my @side = split(/\|\|/, $payload, -1);
        die "[ERROR] Malformed structured canonical split key '$split'.\n"
            unless @side == 2;
        return (_decode_taxon_set_key($side[0]), _decode_taxon_set_key($side[1]));
    }

    my @side = split(/\|\|/, $split, -1);
    die "[ERROR] Malformed legacy canonical split key '$split'.\n"
        unless @side == 2;
    my @left  = grep { $_ ne '' } split(/\.\./, $side[0], -1);
    my @right = grep { $_ ne '' } split(/\.\./, $side[1], -1);
    return (\@left, \@right);
}

sub canonicalize_split {
    my ($split) = @_;
    my ($left, $right) = split_taxon_sets($split);
    return canonical_split($left, $right);
}

sub edge_records {
    my ($root) = @_;
    my @all_tips = sort(tip_labels($root));
    my %descendant_tips;
    my $cache_descendants;
    $cache_descendants = sub {
        my ($node) = @_;
        my @tips;
        if (@{$node->{children}}) {
            push @tips, @{ $cache_descendants->($_) } for @{$node->{children}};
            @tips = sort @tips;
        } else {
            @tips = ($node->{label});
        }
        $descendant_tips{refaddr($node)} = \@tips;
        return \@tips;
    };
    $cache_descendants->($root);

    my @records;

    for my $node (_all_nodes_preorder($root)) {
        next if $node == $root;
        my @left = @{ $descendant_tips{refaddr($node)} };
        my %left = map { $_ => 1 } @left;
        my @right = grep { !$left{$_} } @all_tips;
        push @records, {
            node      => $node,
            branch_id => $node->{branch_id},
            split     => canonical_split(\@left, \@right),
            type      => @{$node->{children}} ? 'internal' : 'terminal',
            left      => \@left,
            right     => \@right,
            value     => $node->{length_raw},
        };
    }

    return \@records;
}

sub retained_primitive_axis {
    my ($root) = @_;
    my $edges = edge_records($root);
    my %winner;
    for my $edge (@{$edges}) {
        my $branch_id = $edge->{branch_id};
        die "[ERROR] A non-root species-tree edge lacks a B alias.\n"
            unless defined $branch_id && $branch_id =~ /^B\d+$/;
        if (!exists $winner{$edge->{split}}
            || branch_id_num($branch_id) < branch_id_num($winner{$edge->{split}})) {
            $winner{$edge->{split}} = $branch_id;
        }
    }
    my @axis = map {
        { branch_id => $winner{$_}, split => $_ }
    } sort {
        branch_id_num($winner{$a}) <=> branch_id_num($winner{$b})
            || $a cmp $b
    } keys %winner;
    return \@axis;
}

sub gene_split_records {
    my ($root) = @_;
    my $edges = edge_records($root);
    my %by_split;
    for my $edge (@{$edges}) {
        push @{ $by_split{$edge->{split}} }, $edge;
    }

    my @records;
    for my $split (sort keys %by_split) {
        my @members = @{ $by_split{$split} };
        my $value;

        if (@members == 1) {
            $value = defined $members[0]{value} && is_numeric_value($members[0]{value})
                ? $members[0]{value}
                : 'NA';
        } else {
            my $sum = 0;
            my $complete = 1;
            for my $member (@members) {
                if (!defined $member->{value} || !is_numeric_value($member->{value})) {
                    $complete = 0;
                    last;
                }
                $sum += 0 + $member->{value};
            }
            $value = $complete ? $sum : 'NA';
        }

        push @records, {
            split           => $split,
            value           => $value,
            component_count => scalar(@members),
        };
    }
    return \@records;
}

sub serialize_plain {
    my ($root) = @_;
    return _serialize_node($root, sub {
        my ($node, $is_root) = @_;
        my $label = @{$node->{children}} ? '' : $node->{label};
        my $length = (!$is_root && defined $node->{length_raw}) ? ':' . $node->{length_raw} : '';
        return ($label, $length);
    }) . ';';
}

sub serialize_for_split {
    my ($root, $record_id) = @_;
    $record_id //= 'species_tree';
    my $tree = _serialize_node($root, sub {
        my ($node, $is_root) = @_;
        my $label = @{$node->{children}} ? '' : $node->{label};
        my $branch = !$is_root ? ':' . ($node->{branch_id} // '') : '';
        return ($label, $branch);
    });
    return $record_id . $tree . ';';
}

sub serialize_figtree {
    my ($root) = @_;
    return _serialize_node($root, sub {
        my ($node, $is_root) = @_;
        return ('', '') if $is_root && @{$node->{children}};
        if (@{$node->{children}}) {
            return ('"' . ($node->{branch_id} // '') . '"', '');
        }
        return ($node->{label} . '_' . ($node->{branch_id} // ''), '');
    }) . ';';
}

sub serialize_branch_subtree {
    my ($node) = @_;
    return $node->{label} unless @{$node->{children}};

    my $child_text = sub {
        my ($child) = @_;
        my $body = serialize_branch_subtree($child);
        my $branch_id = $child->{branch_id};
        die "[ERROR] Missing B alias while serializing a species-tree subtree.\n"
            unless defined $branch_id && $branch_id =~ /^B\d+$/;
        return $body . ':' . $branch_id;
    };
    return '(' . join(',', map { $child_text->($_) } @{$node->{children}}) . ')';
}

sub _serialize_node {
    my ($root, $decorator) = @_;
    my $walk;
    $walk = sub {
        my ($node, $is_root) = @_;
        my $body = '';
        if (@{$node->{children}}) {
            $body = '(' . join(',', map { $walk->($_, 0) } @{$node->{children}}) . ')';
        }
        my ($label, $suffix) = $decorator->($node, $is_root);
        return $body . ($label // '') . ($suffix // '');
    };
    return $walk->($root, 1);
}

sub _parse_subtree {
    my ($state) = @_;
    _skip_ws($state);
    my $ch = _peek($state);
    _error($state, 'Unexpected end of input while parsing a subtree.') unless defined $ch;

    my $node = { children => [], label => '', length_raw => undef };
    if ($ch eq '(') {
        $state->{pos}++;
        push @{$node->{children}}, _parse_subtree($state);
        while (1) {
            _skip_ws($state);
            my $next = _peek($state);
            if (defined $next && $next eq ',') {
                $state->{pos}++;
                push @{$node->{children}}, _parse_subtree($state);
                next;
            }
            last;
        }
        _error($state, 'An internal node must contain at least two children.')
            if @{$node->{children}} < 2;
        _skip_ws($state);
        _expect_char($state, ')');
        _skip_ws($state);
        $node->{label} = _read_optional_token($state);
    } else {
        $node->{label} = _read_required_token($state, 'terminal taxon label');
    }

    _skip_ws($state);
    if (defined _peek($state) && _peek($state) eq ':') {
        $state->{pos}++;
        _skip_ws($state);
        my $value = _read_required_token($state, 'branch length or branch tag');
        if ($state->{allow_branch_tags} && $value =~ /^B\d+$/) {
            $node->{branch_id} = $value;
        } elsif (is_numeric_value($value)) {
            $node->{length_raw} = $value;
        } else {
            _error($state, "Invalid branch value '$value'; expected a finite numeric value"
                . ($state->{allow_branch_tags} ? ' or B<number>' : '') . '.');
        }
    }

    return $node;
}

sub _read_optional_token {
    my ($state) = @_;
    my $start = $state->{pos};
    while ($state->{pos} < $state->{len}) {
        my $ch = substr($state->{text}, $state->{pos}, 1);
        last if $ch =~ /[\s:,();]/;
        $state->{pos}++;
    }
    return substr($state->{text}, $start, $state->{pos} - $start);
}

sub _read_required_token {
    my ($state, $what) = @_;
    my $token = _read_optional_token($state);
    _error($state, "Missing $what.") if $token eq '';
    return $token;
}

sub _leaf_nodes {
    my ($root) = @_;
    my @nodes;
    my $walk;
    $walk = sub {
        my ($node) = @_;
        if (!@{$node->{children}}) {
            push @nodes, $node;
            return;
        }
        $walk->($_) for @{$node->{children}};
    };
    $walk->($root);
    return @nodes;
}

sub _postorder_internal_nodes {
    my ($root) = @_;
    my @nodes;
    my $walk;
    $walk = sub {
        my ($node) = @_;
        $walk->($_) for @{$node->{children}};
        push @nodes, $node if @{$node->{children}};
    };
    $walk->($root);
    return @nodes;
}

sub _all_nodes_preorder {
    my ($root) = @_;
    my @nodes;
    my $walk;
    $walk = sub {
        my ($node) = @_;
        push @nodes, $node;
        $walk->($_) for @{$node->{children}};
    };
    $walk->($root);
    return @nodes;
}

sub _skip_ws {
    my ($state) = @_;
    while ($state->{pos} < $state->{len}
        && substr($state->{text}, $state->{pos}, 1) =~ /\s/) {
        $state->{pos}++;
    }
}

sub _peek {
    my ($state) = @_;
    return undef if $state->{pos} >= $state->{len};
    return substr($state->{text}, $state->{pos}, 1);
}

sub _expect_char {
    my ($state, $expected) = @_;
    my $actual = _peek($state);
    _error($state, "Expected '$expected' but found " . (defined $actual ? "'$actual'" : 'end of input') . '.')
        unless defined $actual && $actual eq $expected;
    $state->{pos}++;
}

sub _error {
    my ($state, $message) = @_;
    die "[ERROR][$state->{source} at character $state->{pos}] $message\n";
}

sub _slurp {
    my ($path) = @_;
    return read_utf8_file($path);
}

sub _requires_structured_split_key {
    for my $label (@_) {
        return 1 unless _is_legacy_safe_label($label);
    }
    return 0;
}

sub _is_legacy_safe_label {
    my ($label) = @_;
    return 0 if $label =~ /\|/;
    return 0 if $label =~ /\.\./;
    return 0 if $label =~ /\A\./;
    return 0 if $label =~ /\.\z/;
    return 1;
}

sub _label_octets {
    my ($label) = @_;
    return encode_utf8_strict($label, 'canonical taxon label');
}

sub _decode_taxon_set_key {
    my ($key) = @_;
    die "[ERROR] Malformed structured taxon-set key '$key'.\n"
        unless index($key, $STRUCTURED_SET_PREFIX) == 0;
    my $payload = substr($key, length($STRUCTURED_SET_PREFIX));
    return [] if $payload eq '';

    my @taxa;
    for my $hex (split(/,/, $payload, -1)) {
        die "[ERROR] Malformed hexadecimal taxon payload in '$key'.\n"
            unless $hex ne '' && $hex =~ /\A(?:[0-9a-f]{2})+\z/;
        my $bytes = pack('H*', $hex);
        push @taxa, decode_utf8_strict($bytes, 'structured canonical taxon payload');
    }
    return \@taxa;
}

1;
