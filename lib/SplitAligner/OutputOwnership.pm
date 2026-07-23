package SplitAligner::OutputOwnership;

use strict;
use warnings;
use Digest::SHA qw(sha256_hex);
use Exporter qw(import);
use Fcntl qw(:mode);
use File::Basename qw(dirname);
use File::Spec;
use JSON::PP;
use SplitAligner::IOSafety qw(build_destination_record);
use SplitAligner::Provenance qw(file_sha256 read_json_file);
use SplitAligner::TextIO qw(
    decode_filesystem_path_utf8
);

our @EXPORT_OK = qw(
    attach_output_ownership
    build_publication_plan
    object_matches_snapshot
    path_exists
    preflight_output_ownership
    verify_publication_plan
    verify_staged_replacement
);

my $OWNERSHIP_SCHEMA = 'SplitAligner-output-ownership-v1';

sub attach_output_ownership {
    my ($manifest, %arg) = @_;
    die "[ERROR][Output ownership] Manifest must be a hash.\n"
        unless ref($manifest) eq 'HASH';
    die "[ERROR][Output ownership] Manifest already contains an ownership record.\n"
        if exists $manifest->{ownership};

    my $mode = _required_text($arg{mode}, 'owner mode');
    my $label = _required_text($arg{owner_label}, 'owner label');
    my $source_root = _required_text($arg{source_root}, 'ownership source root');
    my ($items, $owner) = _validated_items(
        $arg{items}, $arg{owner_manifest_path});

    my @outputs;
    for my $item (@{$items}) {
        next if $item->{owner_manifest};
        my $source = _source_path($source_root, $item->{source});
        push @outputs, _output_entry(
            path          => $source,
            type          => $item->{type},
            role          => $item->{role},
            shared        => $item->{shared},
            relative_path => _relative_publication_path(
                $item->{destination}, dirname($owner->{destination})),
        );
    }
    @outputs = sort { $a->{relative_path} cmp $b->{relative_path} } @outputs;

    $manifest->{ownership} = {
        schema_version          => $OWNERSHIP_SCHEMA,
        producer                => 'SplitAligner',
        owner_mode              => $mode,
        owner_label             => $label,
        manifest_relative_path  => _relative_publication_path(
            $owner->{destination}, dirname($owner->{destination})),
        manifest_role           => $owner->{role},
        manifest_payload_sha256 => _json_sha256($manifest),
        output_inventory_sha256 => _json_sha256(\@outputs),
        output_count            => scalar(@outputs),
        outputs                 => \@outputs,
    };
    return $manifest;
}

sub preflight_output_ownership {
    my (%arg) = @_;
    my $mode = _required_text($arg{mode}, 'owner mode');
    my $label = _required_text($arg{owner_label}, 'owner label');
    my ($items, $owner) = _validated_items(
        $arg{items}, $arg{owner_manifest_path});
    my $force = $arg{force} ? 1 : 0;

    my @existing_nonshared = grep {
        !$_->{shared} && path_exists($_->{destination})
    } @{$items};

    if (@existing_nonshared && !$force) {
        my $item = $existing_nonshared[0];
        die "[ERROR][Output ownership] Output already exists for role "
            . "'$item->{role}': $item->{destination}. Use --force only for "
            . "an intact prior SplitAligner-owned run.\n";
    }

    my $authorization = {
        mode                => $mode,
        owner_label         => $label,
        owner_manifest_path => $owner->{destination},
        has_prior_owner     => 0,
        items               => $items,
    };
    return $authorization unless @existing_nonshared;

    my $verified = _verify_existing_owner(
        mode                => $mode,
        owner_label         => $label,
        owner_manifest_path => $owner->{destination},
        items               => $items,
    );
    $authorization->{has_prior_owner} = 1;
    $authorization->{owner_manifest_sha256}
        = $verified->{owner_manifest_sha256};
    $authorization->{snapshots} = $verified->{snapshots};
    return $authorization;
}

sub build_publication_plan {
    my (%arg) = @_;
    my $source_root = _required_text($arg{source_root}, 'publication source root');
    my $authorization = $arg{authorization};
    die "[ERROR][Output ownership] Missing publication authorization.\n"
        unless ref($authorization) eq 'HASH';
    my $force = $arg{force} ? 1 : 0;
    my ($items) = _validated_items(
        $arg{items}, $authorization->{owner_manifest_path});

    _verify_authorization($authorization) if $authorization->{has_prior_owner};

    my @plan;
    for my $item (@{$items}) {
        my $source = _source_path($source_root, $item->{source});
        my $source_snapshot = _runtime_snapshot($source, $item->{type});
        my $destination = $item->{destination};
        my ($action, $existing_snapshot);

        if (!path_exists($destination)) {
            $action = 'create';
        } elsif ($item->{shared}) {
            my $current = _runtime_snapshot($destination, $item->{type});
            if (_fingerprints_equal(
                    $source_snapshot->{fingerprint}, $current->{fingerprint})) {
                $action = 'reuse';
                $existing_snapshot = $current;
            } else {
                die "[ERROR][Output ownership] Shared output for role "
                    . "'$item->{role}' differs from the generated object: "
                    . "$destination. Shared outputs are immutable and may only "
                    . "be reused byte-for-byte; --force cannot replace them.\n";
            }
        } else {
            die "[ERROR][Output ownership] Output already exists for role "
                . "'$item->{role}': $destination. Use --force only for an "
                . "intact prior SplitAligner-owned run.\n"
                unless $force;
            die "[ERROR][Output ownership] Existing output for role "
                . "'$item->{role}' is unowned or unverified: $destination. "
                . "Move or remove the conflicting path manually; --force "
                . "cannot authorize it.\n"
                unless $authorization->{has_prior_owner}
                    && exists $authorization->{snapshots}{$destination};
            $existing_snapshot = _runtime_snapshot($destination, $item->{type});
            _snapshot_changed($item, $destination)
                unless _snapshots_equal(
                    $existing_snapshot,
                    $authorization->{snapshots}{$destination});
            $action = _fingerprints_equal(
                $source_snapshot->{fingerprint},
                $existing_snapshot->{fingerprint}) ? 'reuse' : 'replace';
        }

        push @plan, {
            %{$item},
            action            => $action,
            source_abs        => $source,
            source_snapshot   => $source_snapshot,
            existing_snapshot => $existing_snapshot,
        };
    }

    return {
        authorization => $authorization,
        items         => \@plan,
    };
}

sub verify_publication_plan {
    my ($plan) = @_;
    die "[ERROR][Output ownership] Invalid publication plan.\n"
        unless ref($plan) eq 'HASH' && ref($plan->{items}) eq 'ARRAY';
    my $authorization = $plan->{authorization};
    _verify_authorization($authorization)
        if $authorization->{has_prior_owner};

    for my $item (@{$plan->{items}}) {
        my $source_now = _runtime_snapshot($item->{source_abs}, $item->{type});
        die "[ERROR][Output ownership] Generated transaction source changed "
            . "before publication for role '$item->{role}': "
            . "$item->{source_abs}.\n"
            unless _snapshots_equal($source_now, $item->{source_snapshot});

        my $destination = $item->{destination};
        if ($item->{action} eq 'create') {
            die "[ERROR][Output ownership] An unowned destination appeared "
                . "after preflight for role '$item->{role}': $destination. "
                . "Publication was aborted without replacing it.\n"
                if path_exists($destination);
        } elsif ($item->{action} eq 'reuse') {
            my $current = _runtime_snapshot($destination, $item->{type});
            _snapshot_changed($item, $destination)
                unless _snapshots_equal($current, $item->{existing_snapshot});
            die "[ERROR][Output ownership] Shared output is no longer identical "
                . "to the generated object for role '$item->{role}': "
                . "$destination.\n"
                unless _fingerprints_equal(
                    $source_now->{fingerprint}, $current->{fingerprint});
        } elsif ($item->{action} eq 'replace') {
            my $current = _runtime_snapshot($destination, $item->{type});
            _snapshot_changed($item, $destination)
                unless _snapshots_equal($current, $item->{existing_snapshot});
        } else {
            die "[ERROR][Output ownership] Unknown publication action "
                . "'$item->{action}'.\n";
        }
    }
    return 1;
}

sub verify_staged_replacement {
    my ($item, $staged_path) = @_;
    die "[ERROR][Output ownership] Cannot verify a non-replacement backup.\n"
        unless ref($item) eq 'HASH' && ($item->{action} // '') eq 'replace';
    my $staged = _runtime_snapshot($staged_path, $item->{type});
    die "[ERROR][Output ownership] The object moved to transactional backup "
        . "does not match the authorized prior output for role "
        . "'$item->{role}': $item->{destination}.\n"
        unless _snapshots_equal($staged, $item->{existing_snapshot});
    return 1;
}

sub object_matches_snapshot {
    my ($path, $type, $snapshot) = @_;
    return 0 unless path_exists($path) && ref($snapshot) eq 'HASH';
    my $current = eval { _runtime_snapshot($path, $type) };
    return 0 if $@ || !defined $current;
    return _snapshots_equal($current, $snapshot);
}

sub path_exists {
    my ($path) = @_;
    return scalar(lstat($path)) ? 1 : 0;
}

sub _verify_authorization {
    my ($authorization) = @_;
    my $verified = _verify_existing_owner(
        mode                => $authorization->{mode},
        owner_label         => $authorization->{owner_label},
        owner_manifest_path => $authorization->{owner_manifest_path},
        items               => $authorization->{items},
    );
    die "[ERROR][Output ownership] Ownership manifest changed after "
        . "preflight: $authorization->{owner_manifest_path}.\n"
        unless $verified->{owner_manifest_sha256}
            eq ($authorization->{owner_manifest_sha256} // '');
    for my $path (sort keys %{$authorization->{snapshots}}) {
        die "[ERROR][Output ownership] Authorized output set changed after "
            . "preflight: $path.\n"
            unless exists $verified->{snapshots}{$path}
                && _snapshots_equal(
                    $verified->{snapshots}{$path},
                    $authorization->{snapshots}{$path});
    }
    return 1;
}

sub _verify_existing_owner {
    my (%arg) = @_;
    my $mode = $arg{mode};
    my $label = $arg{owner_label};
    my ($items, $owner) = _validated_items(
        $arg{items}, $arg{owner_manifest_path});
    my $owner_path = $owner->{destination};

    die "[ERROR][Output ownership] Existing outputs have no ownership "
        . "manifest at $owner_path. Move or remove the conflicting paths "
        . "manually; --force cannot authorize them.\n"
        unless path_exists($owner_path);
    my $owner_snapshot = _runtime_snapshot($owner_path, 'file');
    my $manifest = eval { read_json_file($owner_path) };
    if ($@) {
        die "[ERROR][Output ownership] Invalid ownership manifest "
            . "$owner_path: $@";
    }
    my $ownership = $manifest->{ownership};
    die "[ERROR][Output ownership] Manifest $owner_path is legacy or lacks "
        . "a complete output ownership inventory. Regenerate in a clean "
        . "namespace or remove conflicts manually; --force is refused.\n"
        unless ref($ownership) eq 'HASH';

    _require_exact_keys($ownership, qw(
        schema_version producer owner_mode owner_label
        manifest_relative_path manifest_role manifest_payload_sha256
        output_inventory_sha256 output_count outputs
    ));
    die "[ERROR][Output ownership] Unsupported ownership schema in "
        . "$owner_path.\n"
        unless ($ownership->{schema_version} // '') eq $OWNERSHIP_SCHEMA;
    die "[ERROR][Output ownership] Invalid ownership producer in "
        . "$owner_path.\n"
        unless ($ownership->{producer} // '') eq 'SplitAligner';
    die "[ERROR][Output ownership] Ownership mode mismatch in "
        . "$owner_path.\n"
        unless ($ownership->{owner_mode} // '') eq $mode;
    die "[ERROR][Output ownership] Ownership label mismatch in "
        . "$owner_path.\n"
        unless ($ownership->{owner_label} // '') eq $label;
    die "[ERROR][Output ownership] Ownership manifest role mismatch in "
        . "$owner_path.\n"
        unless ($ownership->{manifest_role} // '') eq $owner->{role};
    my $manifest_relative = _relative_publication_path(
        $owner_path, dirname($owner_path));
    die "[ERROR][Output ownership] Ownership manifest path mismatch in "
        . "$owner_path.\n"
        unless ($ownership->{manifest_relative_path} // '') eq $manifest_relative;

    my %payload = %{$manifest};
    delete $payload{ownership};
    die "[ERROR][Output ownership] Manifest payload binding failed for "
        . "$owner_path.\n"
        unless ($ownership->{manifest_payload_sha256} // '')
            eq _json_sha256(\%payload);
    die "[ERROR][Output ownership] Invalid output inventory in "
        . "$owner_path.\n"
        unless ref($ownership->{outputs}) eq 'ARRAY';
    die "[ERROR][Output ownership] Output inventory count mismatch in "
        . "$owner_path.\n"
        unless 0 + ($ownership->{output_count} // -1)
            == scalar(@{$ownership->{outputs}});
    die "[ERROR][Output ownership] Output inventory checksum mismatch in "
        . "$owner_path.\n"
        unless ($ownership->{output_inventory_sha256} // '')
            eq _json_sha256($ownership->{outputs});

    my %stored;
    for my $entry (@{$ownership->{outputs}}) {
        die "[ERROR][Output ownership] Invalid output inventory entry in "
            . "$owner_path.\n"
            unless ref($entry) eq 'HASH';
        my $relative = $entry->{relative_path} // '';
        die "[ERROR][Output ownership] Duplicate or empty ownership path "
            . "'$relative' in $owner_path.\n"
            if $relative eq '' || exists $stored{$relative};
        $stored{$relative} = $entry;
    }

    my %snapshots = (
        $owner_path => $owner_snapshot,
    );
    my $expected_count = 0;
    for my $item (@{$items}) {
        next if $item->{owner_manifest};
        $expected_count++;
        my $relative = _relative_publication_path(
            $item->{destination}, dirname($owner_path));
        die "[ERROR][Output ownership] Ownership inventory in $owner_path "
            . "does not cover role '$item->{role}' at $item->{destination}.\n"
            unless exists $stored{$relative};
        my $current_entry = _output_entry(
            path          => $item->{destination},
            type          => $item->{type},
            role          => $item->{role},
            shared        => $item->{shared},
            relative_path => $relative,
        );
        die "[ERROR][Output ownership] Owned output is missing, altered, "
            . "type-swapped, hard-linked, or contains unknown descendants "
            . "for role '$item->{role}': $item->{destination}.\n"
            unless _json_equal($current_entry, $stored{$relative});
        $snapshots{$item->{destination}}
            = _runtime_snapshot($item->{destination}, $item->{type});
    }
    die "[ERROR][Output ownership] Ownership inventory in $owner_path "
        . "contains unexpected output paths.\n"
        unless $expected_count == scalar(keys %stored);

    return {
        owner_manifest_sha256 => file_sha256($owner_path),
        snapshots             => \%snapshots,
    };
}

sub _validated_items {
    my ($items, $owner_manifest_path) = @_;
    die "[ERROR][Output ownership] Publication items must be an array.\n"
        unless ref($items) eq 'ARRAY' && @{$items};
    my $owner_path = _canonical_path(
        _required_text($owner_manifest_path, 'owner manifest path'), 'file');
    my (@validated, %seen_destination);
    my $owner_count = 0;

    for my $raw (@{$items}) {
        die "[ERROR][Output ownership] Invalid publication item.\n"
            unless ref($raw) eq 'HASH';
        my $role = _required_text($raw->{role}, 'publication role');
        my $source = _required_text($raw->{source}, "source for $role");
        my $type = $raw->{type} // '';
        die "[ERROR][Output ownership] Invalid publication type for role "
            . "'$role'.\n"
            unless $type eq 'file' || $type eq 'directory';
        die "[ERROR][Output ownership] Publication source must be a safe "
            . "relative path for role '$role': $source.\n"
            if File::Spec->file_name_is_absolute($source)
                || $source =~ m{(?:\A|[\\/])\.\.(?:[\\/]|\z)}
                || $source =~ /[\x00-\x1f]/;
        my $destination = _canonical_path(
            _required_text($raw->{destination}, "destination for $role"), $type);
        die "[ERROR][Output ownership] Duplicate publication destination: "
            . "$destination.\n"
            if $seen_destination{$destination}++;
        my $is_owner = $raw->{owner_manifest} ? 1 : 0;
        $owner_count += $is_owner;
        die "[ERROR][Output ownership] Owner manifest must be a non-shared file.\n"
            if $is_owner && ($type ne 'file' || $raw->{shared});
        push @validated, {
            %{$raw},
            role           => $role,
            source         => $source,
            destination    => $destination,
            type           => $type,
            shared         => $raw->{shared} ? 1 : 0,
            owner_manifest => $is_owner,
        };
    }
    die "[ERROR][Output ownership] Exactly one owner manifest item is required.\n"
        unless $owner_count == 1;
    my ($owner) = grep { $_->{owner_manifest} } @validated;
    die "[ERROR][Output ownership] Owner manifest destination mismatch: "
        . "$owner->{destination} vs $owner_path.\n"
        unless $owner->{destination} eq $owner_path;
    return (\@validated, $owner);
}

sub _source_path {
    my ($root, $relative) = @_;
    my $path = File::Spec->catfile($root, $relative);
    return _canonical_path($path, 'file');
}

sub _canonical_path {
    my ($path, $type) = @_;
    return build_destination_record(
        role => 'output ownership path',
        path => $path,
        type => $type,
    )->{canonical_path};
}

sub _relative_publication_path {
    my ($path, $base) = @_;
    my $relative = File::Spec->canonpath(File::Spec->abs2rel($path, $base));
    die "[ERROR][Output ownership] Publication path is not relative: $path.\n"
        if File::Spec->file_name_is_absolute($relative);
    die "[ERROR][Output ownership] Invalid publication path: $relative.\n"
        if $relative eq '' || $relative =~ /[\x00-\x1f]/;
    return $relative;
}

sub _output_entry {
    my (%arg) = @_;
    return {
        role          => $arg{role},
        relative_path => $arg{relative_path},
        type          => $arg{type},
        shared        => $arg{shared} ? 1 : 0,
        fingerprint   => _fingerprint($arg{path}, $arg{type}),
    };
}

sub _runtime_snapshot {
    my ($path, $expected_type) = @_;
    my @stat = lstat($path);
    die "[ERROR][Output ownership] Expected object is missing: $path.\n"
        unless @stat;
    return {
        device      => 0 + $stat[0],
        inode       => 0 + $stat[1],
        fingerprint => _fingerprint_from_stat($path, $expected_type, \@stat),
    };
}

sub _fingerprint {
    my ($path, $expected_type) = @_;
    my @stat = lstat($path);
    die "[ERROR][Output ownership] Expected object is missing: $path.\n"
        unless @stat;
    return _fingerprint_from_stat($path, $expected_type, \@stat);
}

sub _fingerprint_from_stat {
    my ($path, $expected_type, $stat) = @_;
    my $actual_type = _mode_type($stat->[2]);
    die "[ERROR][Output ownership] Expected $expected_type but found "
        . "$actual_type at $path.\n"
        unless $actual_type eq $expected_type;
    my $mode = sprintf('%04o', $stat->[2] & 07777);

    if ($actual_type eq 'file') {
        die "[ERROR][Output ownership] Hard-linked output is not eligible for "
            . "ownership: $path.\n"
            unless $stat->[3] == 1;
        return {
            type   => 'file',
            mode   => $mode,
            size   => 0 + $stat->[7],
            nlink  => 0 + $stat->[3],
            sha256 => file_sha256($path),
        };
    }

    my @entry = _directory_entries($path);
    @entry = sort { $a->{relative_path} cmp $b->{relative_path} } @entry;
    return {
        type        => 'directory',
        mode        => $mode,
        entry_count => scalar(@entry),
        tree_sha256 => _json_sha256(\@entry),
        entries     => \@entry,
    };
}

sub _directory_entries {
    my ($root) = @_;
    my @entry;
    _walk_directory($root, '', \@entry);
    return @entry;
}

sub _walk_directory {
    my ($directory, $relative_parent, $entry) = @_;
    opendir my $dh, $directory
        or die "[ERROR][Output ownership] Cannot inspect output directory "
            . "$directory: $!\n";
    my @name;
    while (defined(my $raw = readdir $dh)) {
        next if $raw eq '.' || $raw eq '..';
        push @name, decode_filesystem_path_utf8(
            $raw, 'output inventory directory entry');
    }
    closedir $dh
        or die "[ERROR][Output ownership] Cannot close output directory "
            . "$directory: $!\n";

    for my $name (sort { $a cmp $b } @name) {
        my $candidate = File::Spec->catfile($directory, $name);
        my $relative = $relative_parent eq ''
            ? $name
            : File::Spec->catfile($relative_parent, $name);
        $relative = File::Spec->canonpath($relative);
        my @child_stat = lstat($candidate);
        die "[ERROR][Output ownership] Cannot inspect directory entry: "
            . "$candidate.\n"
            unless @child_stat;
        my $child_type = _mode_type($child_stat[2]);
        my $child = {
            relative_path => $relative,
            type          => $child_type,
            mode          => sprintf('%04o', $child_stat[2] & 07777),
        };
        if ($child_type eq 'file') {
            die "[ERROR][Output ownership] Hard-linked directory descendant "
                . "is not eligible for ownership: $candidate.\n"
                unless $child_stat[3] == 1;
            $child->{size} = 0 + $child_stat[7];
            $child->{nlink} = 0 + $child_stat[3];
            $child->{sha256} = file_sha256($candidate);
        } elsif ($child_type eq 'symlink') {
            my $target = readlink($candidate);
            die "[ERROR][Output ownership] Cannot read symlink target: "
                . "$candidate.\n"
                unless defined $target;
            $child->{target} = decode_filesystem_path_utf8(
                $target, 'output inventory symlink target');
        } elsif ($child_type eq 'directory') {
            _walk_directory($candidate, $relative, $entry);
        } else {
            die "[ERROR][Output ownership] Special output object is not "
                . "supported: $candidate.\n";
        }
        push @{$entry}, $child;
    }
}

sub _mode_type {
    my ($mode) = @_;
    return 'file'      if S_ISREG($mode);
    return 'directory' if S_ISDIR($mode);
    return 'symlink'   if S_ISLNK($mode);
    return 'special';
}

sub _snapshot_changed {
    my ($item, $destination) = @_;
    die "[ERROR][Output ownership] Previously authorized output changed "
        . "before publication for role '$item->{role}': $destination. "
        . "--force cannot replace it.\n";
}

sub _snapshots_equal {
    my ($left, $right) = @_;
    return _json_equal($left, $right);
}

sub _fingerprints_equal {
    my ($left, $right) = @_;
    return _json_equal($left, $right);
}

sub _json_equal {
    my ($left, $right) = @_;
    return _canonical_json($left) eq _canonical_json($right);
}

sub _json_sha256 {
    my ($value) = @_;
    return sha256_hex(_canonical_json($value));
}

sub _canonical_json {
    my ($value) = @_;
    return JSON::PP->new->canonical(1)->utf8(1)->encode($value);
}

sub _require_exact_keys {
    my ($hash, @expected) = @_;
    my %expected = map { $_ => 1 } @expected;
    for my $key (keys %{$hash}) {
        die "[ERROR][Output ownership] Unexpected ownership field '$key'.\n"
            unless $expected{$key};
    }
    for my $key (@expected) {
        die "[ERROR][Output ownership] Missing ownership field '$key'.\n"
            unless exists $hash->{$key};
    }
}

sub _required_text {
    my ($value, $name) = @_;
    die "[ERROR][Output ownership] Missing $name.\n"
        unless defined $value && $value ne '';
    return $value;
}

1;
