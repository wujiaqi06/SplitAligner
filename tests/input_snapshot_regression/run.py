#!/usr/bin/env python3
"""The bytes read by helpers must be the bytes recorded in input provenance.

A test-only Perl hook changes an original file at a deterministic read boundary
and restores it afterwards. No product timing hooks or sleeps are required.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import sys

DRIVER = Path(__file__).resolve().parents[2] / "SplitAligner.pl"
REF = "((A,B),C,(D,E));\n"
GENES = "g((A:1,B:1):3,C:1,(D:1,E:1):2);\n"
REPLACEMENT = "g((A:7,B:7):9,C:7,(D:7,E:7):8);\n"
HOOK = r'''
package SnapshotRace;
use strict;
use warnings;
my ($changed, $old);
sub change {
    open my $fh, '<:raw', $ENV{SA_RACE_INPUT} or die $!;
    { local $/; $old = <$fh>; }
    close $fh or die $!;
    open $fh, '>:raw', $ENV{SA_RACE_INPUT} or die $!;
    print {$fh} $ENV{SA_RACE_NEW} or die $!;
    close $fh or die $!;
    open $fh, '>>', $ENV{SA_RACE_RECEIPT} or die $!;
    print {$fh} "changed\n";
    close $fh or die $!;
    $changed = 1;
}
sub restore {
    return unless $changed && !$ENV{SA_RACE_PERSIST};
    open my $fh, '>:raw', $ENV{SA_RACE_INPUT} or die $!;
    print {$fh} $old or die $!;
    close $fh or die $!;
    $changed = 0;
}
sub import {
    if ($ENV{SA_RACE_HOOK} eq 'helper' && $0 =~ /SplitAligner\.pl$/) {
        # v2 computes in-process: inject after snapshots exist and immediately
        # around the actual compact stage, not an obsolete helper subprocess.
        require SplitAligner::CompactPipeline;
        my $real = \&SplitAligner::CompactPipeline::run_compact_matrix;
        no warnings 'redefine';
        *SplitAligner::CompactPipeline::run_compact_matrix = sub {
            change();
            my $result = $real->(@_);
            restore();
            return $result;
        };
    }
}
CHECK {
    if ($ENV{SA_RACE_HOOK} eq 'copy' && $0 =~ /SplitAligner\.pl$/) {
        my $copy = \&main::copy;
        no warnings 'redefine';
        *main::copy = sub {
            if ($_[1] =~ /\Q$ENV{SA_RACE_DEST}\E$/) {
                change();
                my $result = $copy->(@_);
                restore();
                return $result;
            }
            return $copy->(@_);
        };
    }
}
END { restore(); }
1;
'''


def cli(root, *args, race=None, ok=True):
    env = {**os.environ, "LC_ALL": "C", "TMPDIR": str(root / "tmp")}
    if race:
        env.update(race)
        env["PERL5OPT"] = "-MSnapshotRace"
        env["PERL5LIB"] = os.pathsep.join((str(root), str(DRIVER.parent / "lib")))
    command = ["perl", str(DRIVER), *args]
    result = subprocess.run(command, cwd=root, env=env, capture_output=True, text=True)
    stem = root / f"command-{len(list(root.glob('command-*.json'))) + 1:03d}"
    stem.with_suffix(".json").write_text(json.dumps({"command": command, "exit": result.returncode}, indent=2) + "\n")
    stem.with_suffix(".stdout").write_text(result.stdout)
    stem.with_suffix(".stderr").write_text(result.stderr)
    assert (result.returncode == 0) == ok, result.stdout + result.stderr
    assert not list((root / "tmp").glob(".splitaligner-inputs-*")), "snapshot leaked"
    return result


def case(root, name):
    path = root / name
    path.mkdir()
    (path / "tmp").mkdir()
    (path / "SnapshotRace.pm").write_text(HOOK)
    (path / "reference").write_text(REF)
    (path / "genes").write_text(GENES)
    return path


def race_env(root, source, replacement, hook="helper", dest="", persist=False):
    return {"SA_RACE_INPUT": str(root / source), "SA_RACE_NEW": replacement,
            "SA_RACE_HOOK": hook, "SA_RACE_DEST": dest,
            "SA_RACE_PERSIST": "1" if persist else "0",
            "SA_RACE_RECEIPT": str(root / "race.receipt")}


def matrix(root, label, race=None, ok=True):
    return cli(root, "--mode", "matrix", "--species", "reference", "--gene", "genes",
               "--label", label, race=race, ok=ok)


def verify(root):
    transient = case(root, "transient")
    matrix(transient, "result", race_env(transient, "genes", REPLACEMENT))
    assert (transient / "race.receipt").read_text() == "changed\n"
    assert (transient / "genes").read_text() == GENES
    meta = json.loads((transient / "result.run_manifest.json").read_text())
    assert meta["inputs"]["gene_trees"]["sha256"] == hashlib.sha256(GENES.encode()).hexdigest()
    with (transient / "result.matrix_no_fuse.txt").open() as fh:
        row = next(csv.DictReader(fh, delimiter="\t"))
    assert [row[f"B{i}"] for i in range(1, 8)] == ["1"] * 5 + ["3", "2"], row

    persistent = case(root, "persistent")
    result = matrix(persistent, "result", race_env(persistent, "genes", REPLACEMENT, persist=True), ok=False)
    assert (persistent / "race.receipt").read_text() == "changed\n"
    assert "[Input immutability]" in result.stderr, result.stderr
    assert not (persistent / "result.run_manifest.json").exists()
    assert not (persistent / "result.matrix_no_fuse.txt").exists()

    # Finalize stages the already-validated snapshot, even if the original
    # matrix is changed and restored exactly around File::Copy's read.
    for mode in ("finalize_fix", "finalize"):
        path = case(root, mode)
        matrix(path, "fix")
        if mode == "finalize":
            matrix(path, "free")
        original = (path / "fix.matrix_with_fuse.txt").read_text()
        replacement = original.replace("g\t1\t", "g\t7\t")
        assert replacement != original
        args = ["--mode", mode, "--fix", "fix.matrix_with_fuse.txt", "--final_label", "final"]
        if mode == "finalize":
            args.extend(["--free", "free.matrix_with_fuse.txt"])
        cli(path, *args, race=race_env(path, "fix.matrix_with_fuse.txt", replacement,
                                     hook="copy", dest="fix.matrix_with_fuse.txt"))
        assert (path / "race.receipt").read_text() == "changed\n"
        assert (path / "fix.matrix_with_fuse.txt").read_text() == original
        with (path / "final.fix.na_classified.txt").open() as fh:
            row = next(csv.DictReader(fh, delimiter="\t"))
        assert row["B1"] == "1", row
        meta = json.loads((path / "final.finalize_manifest.json").read_text())
        assert meta["inputs"]["fix_matrix_sha256"] == hashlib.sha256(original.encode()).hexdigest()

    print("PASS: 4 transient/persistent matrix and finalize snapshot cases; snapshots cleaned", flush=True)
    # Also execute the current v2 before/after, species/gene and fresh/force
    # boundaries. Keep its independent command logs and eight-case receipt.
    subprocess.run([sys.executable, str(DRIVER.parent / "tests/compact_pipeline/input_snapshot.py"),
                    "--source", str(DRIVER.parent), "--work-dir", str(root / "compact-boundaries")],
                   check=True, env={**os.environ, "LC_ALL": "C", "LANG": "C"})
    receipt = json.loads((root / "compact-boundaries/RECEIPT.json").read_text())
    assert receipt["status"] == "PASS" and len(receipt["cases"]) == 8


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", type=Path, help="new directory retaining all fixture and injection evidence")
    args = parser.parse_args()
    if args.work_dir:
        args.work_dir.mkdir()
        verify(args.work_dir.resolve())
    else:
        with tempfile.TemporaryDirectory(prefix="perl-snapshots-") as tmp:
            verify(Path(tmp))
    print("PASS: 12 snapshot boundary cases, including real compact-stage injection and both finalizers")
