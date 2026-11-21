#!/usr/bin/env python3
"""An old intact run stays protected, with an actionable namespace hint."""
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import subprocess
import tarfile
import tempfile

HERE = Path(__file__).resolve().parent
DRIVER = HERE.parents[1] / "SplitAligner.pl"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def snapshot(root):
    return {str(p.relative_to(root)): digest(p) for p in root.rglob("*") if p.is_file()}


def cli(root, *args, ok=True):
    result = subprocess.run(["perl", str(DRIVER), *args], cwd=root,
                            env={**os.environ, "LC_ALL": "C"},
                            text=True, capture_output=True)
    assert (result.returncode == 0) == ok, result.stdout + result.stderr
    return result


def matrix(label, force=False):
    return ["--mode", "matrix", "--species", "reference.nwk", "--gene", "genes.nwk",
            "--label", label] + (["--force"] if force else [])


def finalize(label, force=False):
    return ["--mode", "finalize", "--free", "free.matrix_with_fuse.txt",
            "--fix", "fix.matrix_with_fuse.txt", "--final_label", label] + (["--force"] if force else [])


with tempfile.TemporaryDirectory(prefix="perl-legacy-layout-") as tmp:
    root = Path(tmp)
    provenance = json.loads((HERE / "fixture_provenance.json").read_text())
    archive = HERE / "baseline_25ab8e0.tar.gz"
    assert digest(archive) == provenance["archive_sha256"]
    with tarfile.open(archive) as tf:
        # This checked-in fixture contains only relative regular files/directories.
        for member in tf.getmembers():
            path = PurePosixPath(member.name)
            assert not path.is_absolute() and ".." not in path.parts
            assert member.isfile() or member.isdir()
        tf.extractall(root)
    assert snapshot(root) == provenance["files_sha256"]

    original = snapshot(root)
    result = cli(root, *matrix("free", force=True), ok=False)
    # The first missing role is now the compact ledger, before fusion cells.
    assert "does not cover role 'matrix compact source ledger'" in result.stderr, result.stderr
    assert "free.compact.jsonl" in result.stderr
    assert "--force cannot authorize replacement" in result.stderr
    assert "new --label value" in result.stderr
    assert "--final_label" not in result.stderr
    assert snapshot(root) == original
    print("PASS: real 25ab8e0 matrix rejected with --label guidance; every old file SHA unchanged")

    # The archived paired finalize included Support output; this invocation has
    # a smaller layout. The extra owned paths still prohibit replacement.
    result = cli(root, *finalize("old", force=True), ok=False)
    assert "contains unexpected output paths" in result.stderr
    assert "--force cannot authorize replacement" in result.stderr
    assert "new --final_label value" in result.stderr
    assert snapshot(root) == original
    print("PASS: incompatible finalize layout rejected with --final_label guidance; every old file SHA unchanged")

    # Following either hint succeeds. Identical force reruns at the current
    # version remain allowed and reuse all files byte-for-byte.
    cli(root, *matrix("current"))
    before = snapshot(root)
    cli(root, *matrix("current", force=True))
    assert snapshot(root) == before
    cli(root, *finalize("current_final"))
    before = snapshot(root)
    cli(root, *finalize("current_final", force=True))
    assert snapshot(root) == before
    fixed_args = ["--mode", "finalize_fix", "--fix", "fix.matrix_with_fuse.txt",
                  "--final_label", "current_fixed"]
    cli(root, *fixed_args)
    before = snapshot(root)
    cli(root, *fixed_args, "--force")
    assert snapshot(root) == before
    assert all(digest(root / p) == sha for p, sha in original.items())
    print("PASS: fresh labels work; current matrix/finalize/finalize_fix --force still succeeds and preserves old outputs")
