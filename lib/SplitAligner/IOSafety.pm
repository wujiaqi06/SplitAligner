package SplitAligner::IOSafety;

use strict;
use warnings;
use Digest::SHA ();
use File::Basename qw(dirname basename);
use File::Spec;
use Fcntl qw(O_NONBLOCK O_RDONLY S_ISREG);
use Exporter qw(import);
use SplitAligner::TextIO qw(abs_path_utf8 getcwd_utf8);

our @EXPORT_OK = qw(
    assert_io_namespaces_disjoint
    assert_publication_device_records
    build_destination_record
    build_input_record
    preflight_publication_devices
    verify_input_snapshots
    verify_publication_device_guard
    verify_publication_device_records
);

sub build_input_record {
    my (%arg) = @_;
    my $role = $arg{role} // '';
    my $path = $arg{path};
    die "[ERROR] Input role is required for I/O safety registration.\n"
        if $role eq '';
    die "[ERROR] Input path is required for role '$role'.\n"
        unless defined $path && $path ne '';
    die "[ERROR] Filesystem path contains a NUL byte.\n" if $path =~ /\x00/;

    my $identity = _strict_existing_input_lookup(
        role => $role,
        path => $path,
    );

    return {
        role             => $role,
        path_as_provided => $path,
        display_absolute_path => _absolute_path_as_provided_for_display($path),
        canonical_path   => $identity->{resolved_path},
        resolved_path    => $identity->{resolved_path},
        sha256           => $identity->{sha256},
        size             => $identity->{size},
        device           => $identity->{device},
        inode            => $identity->{inode},
    };
}

sub build_destination_record {
    my (%arg) = @_;
    my $role = $arg{role} // '';
    my $path = $arg{path};
    my $type = $arg{type} // '';
    die "[ERROR] Destination role is required for I/O safety registration.\n"
        if $role eq '';
    die "[ERROR] Destination path is required for role '$role'.\n"
        unless defined $path && $path ne '';
    die "[ERROR] Destination type for role '$role' must be 'file' or 'directory'.\n"
        unless $type eq 'file' || $type eq 'directory';

    my $canonical = _lexical_publication_absolute_path($path);
    my $resolved = _resolve_destination_path($canonical);
    my ($device, $inode) = _device_inode($canonical);

    return {
        role           => $role,
        path           => $path,
        type           => $type,
        canonical_path => $canonical,
        resolved_path  => $resolved,
        device         => $device,
        inode          => $inode,
    };
}

sub assert_io_namespaces_disjoint {
    my (%arg) = @_;
    my $inputs = $arg{inputs};
    my $destinations = $arg{destinations};
    die "[ERROR] I/O safety input registry must be an array.\n"
        unless ref($inputs) eq 'ARRAY';
    die "[ERROR] I/O safety destination registry must be an array.\n"
        unless ref($destinations) eq 'ARRAY';

    my @current_inputs = map { _refresh_input_identity($_) } @{$inputs};
    my @current_destinations = map { _refresh_destination_identity($_) } @{$destinations};

    for my $input (@current_inputs) {
        for my $destination (@current_destinations) {
            my $reason = _input_destination_overlap($input, $destination);
            next unless defined $reason;
            die _input_overlap_error($input, $destination, $reason);
        }
    }

    if (@current_destinations > 1) {
        for my $left_index (0 .. $#current_destinations - 1) {
            for my $right_index ($left_index + 1 .. $#current_destinations) {
                my $left = $current_destinations[$left_index];
                my $right = $current_destinations[$right_index];
                next unless _destinations_overlap($left, $right);
                die _destination_overlap_error($left, $right);
            }
        }
    }
    return 1;
}

sub verify_input_snapshots {
    my ($inputs) = @_;
    die "[ERROR] I/O safety input registry must be an array.\n"
        unless ref($inputs) eq 'ARRAY';

    for my $input (@{$inputs}) {
        my $role = $input->{role} // '<unknown>';
        my $path = $input->{path_as_provided};
        _input_changed($input, 'missing') unless defined $path && $path ne '';
        _input_changed($input, 'NUL-byte') if $path =~ /\x00/;

        my $identity = _strict_existing_input_lookup(
            role     => $role,
            path     => $path,
            snapshot => $input,
        );
        my $resolved = $identity->{resolved_path};
        _input_changed($input, 'resolved-path-changed')
            unless $resolved eq ($input->{resolved_path} // '');

        _input_changed($input, 'device')
            unless defined $identity->{device} && defined $input->{device}
                && $identity->{device} == $input->{device};
        _input_changed($input, 'inode')
            unless defined $identity->{inode} && defined $input->{inode}
                && $identity->{inode} == $input->{inode};

        _input_changed($input, 'content-sha256')
            unless $identity->{sha256} eq ($input->{sha256} // '');
        _input_changed($input, 'size')
            unless $identity->{size} == 0 + ($input->{size} // -1);
    }
    return 1;
}

sub preflight_publication_devices {
    my (%arg) = @_;
    my $snapshot = _publication_device_snapshot(%arg);
    assert_publication_device_records(
        mode         => $snapshot->{mode},
        transaction  => $snapshot->{transaction},
        destinations => $snapshot->{destinations},
    );
    return $snapshot;
}

sub verify_publication_device_guard {
    my (%arg) = @_;
    my $guard = $arg{guard};
    die "[ERROR][Publication device] Publication-device guard is required.\n"
        unless ref($guard) eq 'HASH';

    my $current = _publication_device_snapshot(
        mode             => $guard->{mode},
        transaction_root => $arg{transaction_root},
        destinations     => $arg{destinations},
    );
    verify_publication_device_records(
        initial => $guard,
        current => $current,
    );
    assert_publication_device_records(
        mode         => $current->{mode},
        transaction  => $current->{transaction},
        destinations => $current->{destinations},
    );
    return 1;
}

sub assert_publication_device_records {
    my (%arg) = @_;
    my $mode = $arg{mode} // '';
    my $transaction = $arg{transaction};
    my $destinations = $arg{destinations};
    die "[ERROR][Publication device] Mode is required for device preflight.\n"
        if $mode eq '';
    die "[ERROR][Publication device] Transaction device record is required.\n"
        unless ref($transaction) eq 'HASH';
    die "[ERROR][Publication device] Destination device registry must be an array.\n"
        unless ref($destinations) eq 'ARRAY';
    _require_device_record($transaction, 'transaction parent');

    for my $destination (@{$destinations}) {
        _require_device_record($destination, 'publication destination');
        next if $destination->{device} == $transaction->{device};
        die _publication_device_mismatch_error(
            $mode, $transaction, $destination);
    }
    return 1;
}

sub verify_publication_device_records {
    my (%arg) = @_;
    my $initial = $arg{initial};
    my $current = $arg{current};
    die "[ERROR][Publication device] Initial device snapshot is required.\n"
        unless ref($initial) eq 'HASH';
    die "[ERROR][Publication device] Current device snapshot is required.\n"
        unless ref($current) eq 'HASH';

    my $mode = $initial->{mode} // '';
    die "[ERROR][Publication device] Device snapshot mode changed before publication.\n"
        unless $mode ne '' && $mode eq ($current->{mode} // '');
    my $initial_transaction = $initial->{transaction};
    my $current_transaction = $current->{transaction};
    _require_device_record($initial_transaction, 'initial transaction parent');
    _require_device_record($current_transaction, 'current transaction parent');
    if ($initial_transaction->{canonical_path} ne $current_transaction->{canonical_path}
            || $initial_transaction->{publication_parent_path}
                ne $current_transaction->{publication_parent_path}
            || $initial_transaction->{resolved_publication_parent_path}
                ne $current_transaction->{resolved_publication_parent_path}) {
        die _publication_transaction_parent_change_error(
            $mode, $initial_transaction, $current_transaction,
            'transaction-parent-changed-after-preflight');
    }
    if ($initial_transaction->{device} != $current_transaction->{device}) {
        die _publication_transaction_device_change_error(
            $mode,
            $current_transaction,
            $initial_transaction->{device},
            'transaction-device-changed-after-preflight',
        );
    }

    my $initial_destinations = $initial->{destinations};
    my $current_destinations = $current->{destinations};
    die "[ERROR][Publication device] Initial destination registry must be an array.\n"
        unless ref($initial_destinations) eq 'ARRAY';
    die "[ERROR][Publication device] Current destination registry must be an array.\n"
        unless ref($current_destinations) eq 'ARRAY';

    my %initial_by_key;
    for my $destination (@{$initial_destinations}) {
        _require_device_record($destination, 'initial publication destination');
        my $key = _publication_device_record_key($destination);
        die "[ERROR][Publication device] Duplicate destination in initial device snapshot: $key\n"
            if exists $initial_by_key{$key};
        $initial_by_key{$key} = $destination;
    }
    die "[ERROR][Publication device] Destination registry changed before publication.\n"
        unless @{$initial_destinations} == @{$current_destinations};

    for my $destination (@{$current_destinations}) {
        _require_device_record($destination, 'current publication destination');
        my $key = _publication_device_record_key($destination);
        my $before = delete $initial_by_key{$key};
        die "[ERROR][Publication device] Destination registry changed before publication for role "
            . "'$destination->{role}' path '$destination->{canonical_path}'.\n"
            unless defined $before;
        if ($before->{publication_parent_path}
                ne $destination->{publication_parent_path}
                || $before->{resolved_publication_parent_path}
                    ne $destination->{resolved_publication_parent_path}) {
            die _publication_parent_change_error(
                $mode, $current_transaction, $before, $destination,
                'destination-publication-parent-changed-after-preflight');
        }
        next if $before->{device} == $destination->{device};
        die _publication_device_change_error(
            $mode,
            $current_transaction,
            $destination,
            $before->{device},
            'destination-publication-parent-device-changed-after-preflight',
        );
    }
    die "[ERROR][Publication device] Destination registry changed before publication.\n"
        if keys %initial_by_key;
    return 1;
}

sub _refresh_input_identity {
    my ($input) = @_;
    my %current = %{$input};
    my $path = $input->{path_as_provided};
    if (defined $path && $path ne '' && $path !~ /\x00/) {
        my $identity = _strict_existing_input_lookup(
            role => $input->{role},
            path => $path,
        );
        $current{canonical_path} = $identity->{resolved_path};
        $current{resolved_path} = $identity->{resolved_path};
        @current{qw(device inode sha256 size)} =
            @{$identity}{qw(device inode sha256 size)};
    }
    return \%current;
}

sub _strict_existing_input_lookup {
    my (%arg) = @_;
    my $role = $arg{role} // '<unknown>';
    my $path = $arg{path};
    my $snapshot = $arg{snapshot};

    _strict_input_lookup_failure($role, $path, 'missing', $snapshot)
        unless defined $path && $path ne '';
    _strict_input_lookup_failure($role, $path, 'NUL-byte', $snapshot)
        if $path =~ /\x00/;

    my $fh;
    unless (sysopen($fh, $path, O_RDONLY | O_NONBLOCK)) {
        my $reason = "$!";
        _strict_input_lookup_failure(
            $role, $path, "direct-open-failed: $reason", $snapshot);
    }
    binmode($fh, ':raw') or _strict_input_lookup_failure(
        $role, $path, "raw-mode-failed: $!", $snapshot, $fh);

    my @opened = stat($fh);
    _strict_input_lookup_failure(
        $role, $path, "fstat-failed: $!", $snapshot, $fh)
        unless @opened;
    _strict_input_lookup_failure(
        $role, $path, 'not-regular-file', $snapshot, $fh)
        unless S_ISREG($opened[2]);

    # Cwd::abs_path is used only after strict kernel lookup has accepted the
    # exact provided pathname. The resolved path must name the same open inode.
    my $resolved = abs_path_utf8($path, "$role resolved input path");
    _strict_input_lookup_failure(
        $role, $path, 'canonical-resolution-failed', $snapshot, $fh)
        unless defined $resolved;
    my @resolved = stat($resolved);
    _strict_input_lookup_failure(
        $role, $path, 'resolved-stat-failed', $snapshot, $fh)
        unless @resolved;
    _strict_input_lookup_failure(
        $role, $path, 'resolved-object-not-regular', $snapshot, $fh)
        unless S_ISREG($resolved[2]);
    _strict_input_lookup_failure(
        $role, $path, 'opened-resolved-identity-mismatch', $snapshot, $fh)
        unless $opened[0] == $resolved[0] && $opened[1] == $resolved[1];

    seek($fh, 0, 0) or _strict_input_lookup_failure(
        $role, $path, "seek-failed: $!", $snapshot, $fh);
    my $sha = Digest::SHA->new(256);
    my $hashed = eval { $sha->addfile($fh); 1 };
    _strict_input_lookup_failure(
        $role, $path, "read-for-hash-failed: $@", $snapshot, $fh)
        unless $hashed;
    my $sha256 = $sha->hexdigest;

    my @after = stat($fh);
    _strict_input_lookup_failure(
        $role, $path, "post-hash-fstat-failed: $!", $snapshot, $fh)
        unless @after;
    _strict_input_lookup_failure(
        $role, $path, 'opened-object-changed-during-lookup', $snapshot, $fh)
        unless $opened[0] == $after[0]
            && $opened[1] == $after[1]
            && $opened[7] == $after[7];
    my @resolved_after = stat($resolved);
    _strict_input_lookup_failure(
        $role, $path, 'resolved-object-changed-during-lookup', $snapshot, $fh)
        unless @resolved_after
            && $after[0] == $resolved_after[0]
            && $after[1] == $resolved_after[1]
            && $after[7] == $resolved_after[7];

    close($fh) or _strict_input_lookup_failure(
        $role, $path, "close-failed: $!", $snapshot);
    return {
        resolved_path => $resolved,
        sha256         => $sha256,
        size           => 0 + $after[7],
        device         => 0 + $after[0],
        inode          => 0 + $after[1],
    };
}

sub _strict_input_lookup_failure {
    my ($role, $path, $reason, $snapshot, $fh) = @_;
    close($fh) if defined $fh;
    if (ref($snapshot) eq 'HASH') {
        _input_changed($snapshot, "strict-lookup-$reason");
    }
    $path = '<undefined>' unless defined $path;
    die "[ERROR] Cannot open input path for role '$role' path '$path'. "
        . "Cannot resolve input path under strict POSIX lookup "
        . "(reason: $reason).\n";
}

sub _refresh_destination_identity {
    my ($destination) = @_;
    return build_destination_record(
        role => $destination->{role},
        path => $destination->{canonical_path},
        type => $destination->{type},
    );
}

sub _publication_device_snapshot {
    my (%arg) = @_;
    my $mode = $arg{mode} // '';
    my $transaction_root = $arg{transaction_root};
    my $destinations = $arg{destinations};
    die "[ERROR][Publication device] Mode is required for device preflight.\n"
        if $mode eq '';
    die "[ERROR][Publication device] Transaction parent is required for mode '$mode'.\n"
        unless defined $transaction_root && $transaction_root ne '';
    die "[ERROR][Publication device] Destination registry must be an array for mode '$mode'.\n"
        unless ref($destinations) eq 'ARRAY';

    my $transaction = _transaction_device_record(
        mode => $mode,
        role => 'transaction parent',
        path => $transaction_root,
    );
    my @destination_records = map {
        _destination_device_record(
            mode => $mode,
            role => $_->{role},
            path => $_->{canonical_path},
        )
    } @{$destinations};
    return {
        schema_version => 'SplitAligner-publication-device-guard-v2',
        mode            => $mode,
        transaction     => $transaction,
        destinations    => \@destination_records,
    };
}

sub _transaction_device_record {
    my (%arg) = @_;
    my $mode = $arg{mode};
    my $role = $arg{role} // '';
    my $path = $arg{path};
    die "[ERROR][Publication device] Device role is required for mode '$mode'.\n"
        if $role eq '';
    die "[ERROR][Publication device] Device path is required for role '$role' in mode '$mode'.\n"
        unless defined $path && $path ne '';

    my $canonical = _lexical_publication_absolute_path($path);
    die "[ERROR][Publication device] Transaction parent is not an existing "
        . "directory for mode '$mode': $canonical\n"
        unless -d $canonical;
    my $resolved = abs_path_utf8($canonical, "transaction parent for $role");
    die "[ERROR][Publication device] Cannot resolve transaction parent "
        . "'$canonical' for mode '$mode'.\n" unless defined $resolved;
    my ($device) = _device_inode($resolved);
    die "[ERROR][Publication device] Cannot determine filesystem device for "
        . "transaction parent '$canonical' (resolved via '$resolved') in mode "
        . "'$mode'.\n"
        unless defined $device;

    return {
        role                             => $role,
        path                             => $path,
        canonical_path                   => $canonical,
        publication_parent_path          => $canonical,
        resolved_publication_parent_path => $resolved,
        resolved_path                    => $resolved,
        existing_destination_path        => $canonical,
        device                           => $device,
    };
}

sub _destination_device_record {
    my (%arg) = @_;
    my $mode = $arg{mode};
    my $role = $arg{role} // '';
    my $path = $arg{path};
    die "[ERROR][Publication device] Device role is required for mode '$mode'.\n"
        if $role eq '';
    die "[ERROR][Publication device] Device path is required for role '$role' in mode '$mode'.\n"
        unless defined $path && $path ne '';

    my $canonical = _lexical_publication_absolute_path($path);
    my $publication_parent = dirname($canonical);
    my ($nearest_parent, $resolved_parent) =
        _nearest_existing_publication_parent(
            mode => $mode, role => $role, path => $publication_parent);
    my ($device) = _device_inode($resolved_parent);
    die "[ERROR][Publication device] Cannot determine publication-parent "
        . "filesystem device for role '$role' destination '$canonical' "
        . "(publication parent '$publication_parent', resolved via "
        . "'$resolved_parent') in mode '$mode'.\n" unless defined $device;

    return {
        role                             => $role,
        path                             => $path,
        canonical_path                   => $canonical,
        publication_parent_path          => $publication_parent,
        nearest_existing_parent_path     => $nearest_parent,
        resolved_publication_parent_path => $resolved_parent,
        resolved_path                    => $resolved_parent,
        existing_destination_path        =>
            ((-e $canonical || -l $canonical) ? $canonical : ''),
        device                           => $device,
    };
}

sub _nearest_existing_publication_parent {
    my (%arg) = @_;
    my $mode = $arg{mode};
    my $role = $arg{role};
    my $publication_parent = $arg{path};
    my $probe = $publication_parent;
    while (!-d $probe) {
        if (-e $probe || -l $probe) {
            die "[ERROR][Publication device] Publication parent component is "
                . "not a directory for role '$role' in mode '$mode': $probe\n";
        }
        my $parent = dirname($probe);
        die "[ERROR][Publication device] Cannot find an existing publication "
            . "parent for role '$role' path '$publication_parent' in mode "
            . "'$mode'.\n" if $parent eq $probe;
        $probe = $parent;
    }
    my $resolved = abs_path_utf8(
        $probe, "publication parent device path for $role");
    die "[ERROR][Publication device] Cannot resolve publication parent '$probe' "
        . "for role '$role' in mode '$mode'.\n" unless defined $resolved;
    die "[ERROR][Publication device] Resolved publication parent is not a "
        . "directory for role '$role' in mode '$mode': $resolved\n"
        unless -d $resolved;
    return ($probe, $resolved);
}

sub _publication_device_record_key {
    my ($record) = @_;
    return join("\x1e", $record->{role}, $record->{canonical_path});
}

sub _require_device_record {
    my ($record, $description) = @_;
    die "[ERROR][Publication device] Invalid $description device record.\n"
        unless ref($record) eq 'HASH'
            && defined $record->{role} && $record->{role} ne ''
            && defined $record->{canonical_path} && $record->{canonical_path} ne ''
            && defined $record->{publication_parent_path}
            && $record->{publication_parent_path} ne ''
            && defined $record->{resolved_publication_parent_path}
            && $record->{resolved_publication_parent_path} ne ''
            && exists $record->{existing_destination_path}
            && defined $record->{device} && $record->{device} =~ /\A\d+\z/;
}

sub _publication_device_mismatch_error {
    my ($mode, $transaction, $destination) = @_;
    my $existing = $destination->{existing_destination_path} ne ''
        ? "'$destination->{existing_destination_path}'" : '(absent)';
    return "[ERROR][Publication device] Mode '$mode' requires one publication "
        . "filesystem device. Transaction parent '$transaction->{canonical_path}' "
        . "resolves via '$transaction->{resolved_publication_parent_path}' on "
        . "device $transaction->{device}; destination role "
        . "'$destination->{role}' path '$destination->{canonical_path}' has "
        . "publication parent '$destination->{publication_parent_path}', resolved "
        . "via '$destination->{resolved_publication_parent_path}' on device "
        . "$destination->{device}; existing destination path: $existing. "
        . "Run the invocation on the destination filesystem; for finalization, "
        . "run on the matrix filesystem or copy the matrix run into the "
        . "transaction filesystem. --force cannot override publication-device "
        . "protection.\n";
}

sub _publication_device_change_error {
    my ($mode, $transaction, $destination, $previous_device, $reason) = @_;
    my $existing = $destination->{existing_destination_path} ne ''
        ? "'$destination->{existing_destination_path}'" : '(absent)';
    return "[ERROR][Publication device] Mode '$mode' publication-device identity "
        . "changed after preflight (reason: $reason). Transaction parent "
        . "'$transaction->{canonical_path}' resolves via "
        . "'$transaction->{resolved_publication_parent_path}' on device "
        . "$transaction->{device}; "
        . "destination role '$destination->{role}' path "
        . "'$destination->{canonical_path}' has publication parent "
        . "'$destination->{publication_parent_path}', resolved via "
        . "'$destination->{resolved_publication_parent_path}' on device "
        . "$destination->{device} (preflight device $previous_device); existing "
        . "destination path: $existing. Run the invocation on the "
        . "destination filesystem; for finalization, run on the matrix filesystem "
        . "or copy the matrix run into the transaction filesystem. --force cannot "
        . "override publication-device protection.\n";
}

sub _publication_transaction_device_change_error {
    my ($mode, $transaction, $previous_device, $reason) = @_;
    return "[ERROR][Publication device] Mode '$mode' publication-device identity "
        . "changed after preflight (reason: $reason). Transaction parent "
        . "'$transaction->{canonical_path}' resolves via "
        . "'$transaction->{resolved_publication_parent_path}' on device "
        . "$transaction->{device} "
        . "(preflight device $previous_device). Run the invocation on the "
        . "destination filesystem; for finalization, run on the matrix filesystem "
        . "or copy the matrix run into the transaction filesystem. --force cannot "
        . "override publication-device protection.\n";
}

sub _publication_parent_change_error {
    my ($mode, $transaction, $before, $current, $reason) = @_;
    return "[ERROR][Publication device] Mode '$mode' publication-parent identity "
        . "changed after preflight (reason: $reason). Transaction parent "
        . "'$transaction->{canonical_path}' resolves via "
        . "'$transaction->{resolved_publication_parent_path}' on device "
        . "$transaction->{device}; destination role '$current->{role}' path "
        . "'$current->{canonical_path}' had publication parent "
        . "'$before->{publication_parent_path}' resolved via "
        . "'$before->{resolved_publication_parent_path}' and now has publication "
        . "parent '$current->{publication_parent_path}' resolved via "
        . "'$current->{resolved_publication_parent_path}'. --force cannot override "
        . "publication-device protection.\n";
}

sub _publication_transaction_parent_change_error {
    my ($mode, $before, $current, $reason) = @_;
    return "[ERROR][Publication device] Mode '$mode' transaction-parent identity "
        . "changed after preflight (reason: $reason). Preflight parent "
        . "'$before->{canonical_path}' resolved via "
        . "'$before->{resolved_publication_parent_path}'; current parent "
        . "'$current->{canonical_path}' resolves via "
        . "'$current->{resolved_publication_parent_path}'. --force cannot override "
        . "publication-device protection.\n";
}

sub _input_destination_overlap {
    my ($input, $destination) = @_;
    return 'exact-path'
        if $input->{canonical_path} eq $destination->{canonical_path};
    return 'exact-path'
        if defined $input->{display_absolute_path}
            && $input->{display_absolute_path} eq $destination->{canonical_path};
    return 'resolved-path'
        if defined $input->{resolved_path}
            && defined $destination->{resolved_path}
            && $input->{resolved_path} eq $destination->{resolved_path};
    return 'same-inode' if _same_inode($input, $destination);

    if ($destination->{type} eq 'directory') {
        return 'inside-output-directory'
            if _is_same_or_descendant(
                $input->{canonical_path}, $destination->{canonical_path});
        return 'inside-output-directory'
            if defined $input->{resolved_path}
                && defined $destination->{resolved_path}
                && _is_same_or_descendant(
                    $input->{resolved_path}, $destination->{resolved_path});
    }
    return undef;
}

sub _destinations_overlap {
    my ($left, $right) = @_;
    return 1 if $left->{canonical_path} eq $right->{canonical_path};
    return 1 if defined $left->{resolved_path}
        && defined $right->{resolved_path}
        && $left->{resolved_path} eq $right->{resolved_path};
    return 1 if _same_inode($left, $right);

    if ($left->{type} eq 'directory') {
        return 1 if _is_same_or_descendant(
            $right->{canonical_path}, $left->{canonical_path});
        return 1 if defined $right->{resolved_path}
            && defined $left->{resolved_path}
            && _is_same_or_descendant(
                $right->{resolved_path}, $left->{resolved_path});
    }
    if ($right->{type} eq 'directory') {
        return 1 if _is_same_or_descendant(
            $left->{canonical_path}, $right->{canonical_path});
        return 1 if defined $left->{resolved_path}
            && defined $right->{resolved_path}
            && _is_same_or_descendant(
                $left->{resolved_path}, $right->{resolved_path});
    }
    return 0;
}

sub _lexical_publication_absolute_path {
    my ($path) = @_;
    die "[ERROR] Filesystem path contains a NUL byte.\n" if $path =~ /\x00/;
    my $absolute = File::Spec->file_name_is_absolute($path)
        ? $path
        : File::Spec->catfile(getcwd_utf8(), $path);
    $absolute = File::Spec->canonpath($absolute);

    # SplitAligner targets POSIX filesystems. Resolve lexical dot components
    # without resolving symlinks or applying Unicode normalization.
    die "[ERROR] Expected an absolute POSIX path, got: $absolute\n"
        unless $absolute =~ m{\A/};
    my @component;
    for my $part (split(m{/+}, $absolute, -1)) {
        next if $part eq '' || $part eq '.';
        if ($part eq '..') {
            pop @component if @component;
            next;
        }
        push @component, $part;
    }
    return '/' . join('/', @component);
}

sub _absolute_path_as_provided_for_display {
    my ($path) = @_;
    return $path if File::Spec->file_name_is_absolute($path);
    my $cwd = getcwd_utf8();
    return $cwd eq '/' ? "/$path" : "$cwd/$path";
}

sub _resolve_destination_path {
    my ($canonical) = @_;
    my $resolved = abs_path_utf8($canonical, 'existing publication destination');
    return $resolved if defined $resolved;

    my @tail;
    my $probe = $canonical;
    while (!-e $probe) {
        my $parent = dirname($probe);
        last if $parent eq $probe;
        unshift @tail, basename($probe);
        $probe = $parent;
    }
    my $parent_resolved = abs_path_utf8($probe, 'publication destination parent');
    return $canonical unless defined $parent_resolved;
    return _lexical_publication_absolute_path(
        File::Spec->catfile($parent_resolved, @tail));
}

sub _is_same_or_descendant {
    my ($child, $parent) = @_;
    my @child = _path_components($child);
    my @parent = _path_components($parent);
    return 0 if @child < @parent;
    for my $index (0 .. $#parent) {
        return 0 unless $child[$index] eq $parent[$index];
    }
    return 1;
}

sub _path_components {
    my ($path) = @_;
    return grep { $_ ne '' } split(m{/+}, $path, -1);
}

sub _device_inode {
    my ($path) = @_;
    my @stat = stat($path);
    return (undef, undef) unless @stat;
    return (0 + $stat[0], 0 + $stat[1]);
}

sub _same_inode {
    my ($left, $right) = @_;
    return 0 unless defined $left->{device} && defined $left->{inode};
    return 0 unless defined $right->{device} && defined $right->{inode};
    return $left->{device} == $right->{device}
        && $left->{inode} == $right->{inode};
}

sub _input_overlap_error {
    my ($input, $destination, $reason) = @_;
    return "[ERROR][I/O namespace] Input role '$input->{role}' path "
        . "'$input->{path_as_provided}' overlaps destination role "
        . "'$destination->{role}' path '$destination->{canonical_path}' "
        . "(reason: $reason). --force cannot override input protection.\n";
}

sub _destination_overlap_error {
    my ($left, $right) = @_;
    return "[ERROR][I/O namespace] Destination role '$left->{role}' path "
        . "'$left->{canonical_path}' overlaps destination role '$right->{role}' "
        . "path '$right->{canonical_path}' (reason: duplicate-destination). "
        . "--force cannot override I/O namespace protection.\n";
}

sub _input_changed {
    my ($input, $reason) = @_;
    die "[ERROR][Input immutability] Input role '$input->{role}' path "
        . "'$input->{path_as_provided}' changed before publication "
        . "(reason: $reason). --force cannot override input protection.\n";
}

1;
