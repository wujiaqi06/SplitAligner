#!/usr/bin/env python3
"""RECERT-006 output-ownership and safe-force regression suite."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shutil
import stat
import subprocess
import sys
import tempfile
from typing import Any, Iterable


ROOT = Path(sys.argv[1]).resolve()
SPLITALIGNER = ROOT / "SplitAligner.pl"
PASS_COUNT = 0

COMMON_OUTPUTS = [
    "species_tree.forSplit.nwk",
    "species_tree.FigTree.tre",
    "species_tree.splits.txt",
    "species_tree.branch_map.txt",
    "species_tree.primitive_axis.tsv",
]


def fail(message: str) -> None:
    raise AssertionError(message)


def passed(message: str) -> None:
    global PASS_COUNT
    PASS_COUNT += 1
    print(f"[PASS] {message}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def mode_text(path: Path) -> str:
    return f"{stat.S_IMODE(path.lstat().st_mode):04o}"


def snapshot(path: Path, *, identity: bool = True) -> dict[str, Any]:
    info = path.lstat()
    base: dict[str, Any] = {
        "mode": mode_text(path),
    }
    if identity:
        base.update(device=info.st_dev, inode=info.st_ino)
    if path.is_symlink():
        base.update(type="symlink", target=os.readlink(path))
        return base
    if path.is_file():
        base.update(
            type="file",
            size=info.st_size,
            nlink=info.st_nlink,
            sha256=sha256(path),
        )
        return base
    if path.is_dir():
        entries = []
        for child in sorted(path.iterdir(), key=lambda item: item.name):
            entries.append({"name": child.name, "object": snapshot(child, identity=identity)})
        base.update(type="directory", entries=entries)
        return base
    base.update(type="special")
    return base


def snapshot_paths(paths: Iterable[Path], *, identity: bool = True) -> dict[str, Any]:
    result = {}
    for path in paths:
        result[str(path)] = snapshot(path, identity=identity) if path.exists() or path.is_symlink() else None
    return result


def assert_snapshot(paths: Iterable[Path], before: dict[str, Any], label: str) -> None:
    after = snapshot_paths(paths)
    if before != after:
        fail(f"{label}: protected path snapshot changed")


def write_seed(directory: Path, *, species: str | None = None, genes: str | None = None) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    (directory / "species.nwk").write_text(
        species or "((A:1,B:1):2,(C:1,D:1):3);\n", encoding="utf-8"
    )
    (directory / "genes.nwk").write_text(
        genes or "g((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);\n",
        encoding="utf-8",
    )


def run_sa(directory: Path, args: list[str], *, env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
    merged = os.environ.copy()
    if env:
        merged.update(env)
    return subprocess.run(
        ["perl", str(SPLITALIGNER), *args],
        cwd=directory,
        env=merged,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )


def assert_no_workspace(directory: Path) -> None:
    leftovers = [
        item.name
        for item in directory.iterdir()
        if item.name.startswith(".splitaligner-")
    ]
    if leftovers:
        fail(f"transaction workspaces remained: {leftovers}")


def expect_rejected(directory: Path, args: list[str], label: str) -> subprocess.CompletedProcess[str]:
    result = run_sa(directory, args)
    if result.returncode == 0:
        fail(f"{label}: command unexpectedly succeeded")
    if "Output ownership" not in result.stderr:
        fail(f"{label}: rejection did not identify output ownership\n{result.stderr}")
    assert_no_workspace(directory)
    return result


def matrix_args(label: str, *, force: bool = False, species: str = "species.nwk", genes: str = "genes.nwk") -> list[str]:
    args = ["--mode", "matrix", "--species", species, "--gene", genes, "--label", label]
    if force:
        args.append("--force")
    return args


def matrix(directory: Path, label: str, *, force: bool = False, species: str = "species.nwk", genes: str = "genes.nwk") -> None:
    result = run_sa(directory, matrix_args(label, force=force, species=species, genes=genes))
    if result.returncode != 0:
        fail(f"matrix {label} failed:\n{result.stderr}")
    assert_no_workspace(directory)


def finalize_args(final: str, *, force: bool = False, support: bool = False) -> list[str]:
    args = [
        "--mode", "finalize",
        "--free", "free.matrix_with_fuse.txt",
        "--fix", "fix.matrix_with_fuse.txt",
        "--final_label", final,
    ]
    if support:
        args.extend(["--species_tree", "species_tree.forSplit.nwk"])
    if force:
        args.append("--force")
    return args


def finalize(directory: Path, final: str, *, force: bool = False, support: bool = False, env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
    result = run_sa(directory, finalize_args(final, force=force, support=support), env=env)
    if result.returncode != 0:
        fail(f"finalize {final} failed:\n{result.stderr}")
    assert_no_workspace(directory)
    return result


def finalize_fix_args(final: str, *, force: bool = False) -> list[str]:
    args = [
        "--mode", "finalize_fix",
        "--fix", "fix.matrix_with_fuse.txt",
        "--final_label", final,
    ]
    if force:
        args.append("--force")
    return args


def finalize_fix(directory: Path, final: str, *, force: bool = False) -> None:
    result = run_sa(directory, finalize_fix_args(final, force=force))
    if result.returncode != 0:
        fail(f"finalize_fix {final} failed:\n{result.stderr}")
    assert_no_workspace(directory)


def matrix_outputs(directory: Path, label: str, *, include_common: bool = True) -> list[Path]:
    names = [
        f"{label}_splits",
        f"{label}_split_branch_label",
        f"{label}.gene_id_map.tsv",
        f"{label}.primitive_axis.tsv",
        f"{label}.primitive_state.tsv",
        f"{label}.matrix_no_fuse.txt",
        f"{label}.matrix_with_fuse.txt",
        f"{label}.run_manifest.json",
    ]
    if include_common:
        names = COMMON_OUTPUTS + names
    return [directory / name for name in names]


def make_pair(directory: Path) -> None:
    write_seed(directory)
    matrix(directory, "free")
    matrix(directory, "fix")


def load_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def write_json(path: Path, value: dict[str, Any]) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":")) + "\n",
        encoding="utf-8",
    )


def inventory_paths(manifest: Path) -> set[str]:
    return {entry["relative_path"] for entry in load_json(manifest)["ownership"]["outputs"]}


def first_regular_file(directory: Path) -> Path:
    for path in sorted(directory.rglob("*")):
        if path.is_file() and not path.is_symlink():
            return path
    fail(f"no regular descendant in {directory}")
    raise AssertionError


def case_dir(base: Path, name: str) -> Path:
    path = base / name
    path.mkdir(parents=True)
    return path


def test_unowned_matrix_files(base: Path) -> None:
    names = [
        "x.gene_id_map.tsv",
        "x.primitive_axis.tsv",
        "x.primitive_state.tsv",
        "x.matrix_no_fuse.txt",
        "x.matrix_with_fuse.txt",
        "x.run_manifest.json",
    ]
    for index, name in enumerate(names):
        directory = case_dir(base, f"own01_{index}")
        write_seed(directory)
        target = directory / name
        target.write_text(f"UNOWNED {name}\n", encoding="utf-8")
        paths = list(directory.iterdir())
        before = snapshot_paths(paths)
        expect_rejected(directory, matrix_args("x", force=True), f"OWN-01 {name}")
        assert_snapshot(paths, before, f"OWN-01 {name}")
    passed("OWN-01 arbitrary unowned matrix file destinations are preserved")


def test_unowned_matrix_directories(base: Path) -> None:
    for suffix in ("splits", "split_branch_label"):
        directory = case_dir(base, f"own02_{suffix}")
        write_seed(directory)
        target = directory / f"x_{suffix}"
        (target / "nested").mkdir(parents=True)
        (target / "important.txt").write_text("DO NOT DELETE\n", encoding="utf-8")
        (target / "nested" / "user.txt").write_text("USER DATA\n", encoding="utf-8")
        paths = list(directory.iterdir())
        before = snapshot_paths(paths)
        expect_rejected(directory, matrix_args("x", force=True), f"OWN-02 {suffix}")
        assert_snapshot(paths, before, f"OWN-02 {suffix}")
        swapped = case_dir(base, f"own02_{suffix}_as_file")
        write_seed(swapped)
        target = swapped / f"x_{suffix}"
        target.write_text("UNOWNED TYPE-SWAPPED DIRECTORY DESTINATION\n", encoding="utf-8")
        paths = list(swapped.iterdir())
        before = snapshot_paths(paths)
        expect_rejected(swapped, matrix_args("x", force=True), f"OWN-02 {suffix} as file")
        assert_snapshot(paths, before, f"OWN-02 {suffix} as file")
    passed("OWN-02 arbitrary unowned matrix directory trees are preserved")


def test_unowned_common_files(base: Path) -> None:
    for index, name in enumerate(COMMON_OUTPUTS):
        directory = case_dir(base, f"own03_{index}")
        write_seed(directory)
        target = directory / name
        target.write_text(f"UNOWNED COMMON {name}\n", encoding="utf-8")
        paths = list(directory.iterdir())
        before = snapshot_paths(paths)
        expect_rejected(directory, matrix_args("x", force=True), f"OWN-03 {name}")
        assert_snapshot(paths, before, f"OWN-03 {name}")
    passed("OWN-03 differing unowned common species artifacts are immutable")


def test_unowned_finalize_outputs(base: Path) -> None:
    names = [
        "free.matrix_with_fuse.na_fuse.txt",
        "fix.matrix_with_fuse.na_fuse.txt",
        "out.fix.na_classified.txt",
        "out.free.na_classified.txt",
        "out.support_b.txt",
        "species_tree.support_b.nwk",
        "out.finalize_manifest.json",
    ]
    for index, name in enumerate(names):
        directory = case_dir(base, f"own04_{index}")
        make_pair(directory)
        target = directory / name
        target.write_text(f"UNOWNED FINALIZE {name}\n", encoding="utf-8")
        paths = list(directory.iterdir())
        before = snapshot_paths(paths)
        expect_rejected(
            directory,
            finalize_args("out", force=True, support=True),
            f"OWN-04 {name}",
        )
        assert_snapshot(paths, before, f"OWN-04 {name}")
    passed("OWN-04 arbitrary paired-finalize destinations are preserved")


def test_unowned_finalize_fix_outputs(base: Path) -> None:
    names = [
        "fix.matrix_with_fuse.na_fuse.txt",
        "out.fix.na_classified.txt",
        "out.finalize_manifest.json",
    ]
    for index, name in enumerate(names):
        directory = case_dir(base, f"own05_{index}")
        write_seed(directory)
        matrix(directory, "fix")
        target = directory / name
        target.write_text(f"UNOWNED FINALIZE_FIX {name}\n", encoding="utf-8")
        paths = list(directory.iterdir())
        before = snapshot_paths(paths)
        expect_rejected(directory, finalize_fix_args("out", force=True), f"OWN-05 {name}")
        assert_snapshot(paths, before, f"OWN-05 {name}")
    passed("OWN-05 arbitrary fix-only destinations are preserved")


def test_owner_record_defects(base: Path) -> None:
    actions = (
        "missing",
        "invalid",
        "legacy",
        "payload_tamper",
        "ownership_tamper",
        "owned_file_tamper",
        "owned_file_mode_tamper",
    )
    for action in actions:
        directory = case_dir(base, f"own06_07_15_{action}")
        write_seed(directory)
        matrix(directory, "x")
        manifest = directory / "x.run_manifest.json"
        if action == "missing":
            manifest.unlink()
        elif action == "invalid":
            manifest.write_text("{not-json\n", encoding="utf-8")
        elif action == "legacy":
            data = load_json(manifest)
            data.pop("ownership")
            write_json(manifest, data)
        elif action == "payload_tamper":
            data = load_json(manifest)
            data["label"] = "tampered"
            write_json(manifest, data)
        elif action == "ownership_tamper":
            data = load_json(manifest)
            data["ownership"]["output_inventory_sha256"] = "0" * 64
            write_json(manifest, data)
        elif action == "owned_file_tamper":
            with (directory / "x.matrix_no_fuse.txt").open("a", encoding="utf-8") as handle:
                handle.write("TAMPER\n")
        elif action == "owned_file_mode_tamper":
            (directory / "x.matrix_no_fuse.txt").chmod(0o600)
        paths = list(directory.iterdir())
        before = snapshot_paths(paths)
        expect_rejected(directory, matrix_args("x", force=True), f"owner defect {action}")
        assert_snapshot(paths, before, f"owner defect {action}")
    passed("OWN-06/07/15 missing, invalid, legacy, and tampered owners fail closed")


def test_owned_directory_defects(base: Path) -> None:
    mutations = ("add", "modify", "delete", "type_swap", "mode")
    for dirname in ("x_splits", "x_split_branch_label"):
        for mutation in mutations:
            directory = case_dir(base, f"own08_09_{dirname}_{mutation}")
            write_seed(directory)
            matrix(directory, "x")
            owned = directory / dirname
            if mutation == "add":
                extra = owned / "nested" / "user.txt"
                extra.parent.mkdir()
                extra.write_text("USER DATA\n", encoding="utf-8")
            else:
                target = first_regular_file(owned)
                if mutation == "modify":
                    with target.open("a", encoding="utf-8") as handle:
                        handle.write("TAMPER\n")
                elif mutation == "delete":
                    target.unlink()
                elif mutation == "type_swap":
                    target.unlink()
                    target.mkdir()
                elif mutation == "mode":
                    owned.chmod(0o700)
            paths = list(directory.iterdir())
            before = snapshot_paths(paths)
            expect_rejected(directory, matrix_args("x", force=True), f"directory {dirname} {mutation}")
            assert_snapshot(paths, before, f"directory {dirname} {mutation}")
    passed("OWN-08/09 recursive directory additions, edits, deletions, and type swaps fail closed")


def test_link_substitutions(base: Path) -> None:
    for kind in ("top_symlink", "top_hardlink", "nested_symlink", "nested_hardlink"):
        directory = case_dir(base, f"links_{kind}")
        write_seed(directory)
        matrix(directory, "x")
        target = directory / "unrelated.txt"
        if kind.startswith("top"):
            owned = directory / "x.matrix_no_fuse.txt"
        else:
            owned = first_regular_file(directory / "x_splits")
        original = owned.read_bytes()
        target.write_bytes(original)
        owned.unlink()
        if kind.endswith("symlink"):
            owned.symlink_to(os.path.relpath(target, owned.parent))
        else:
            os.link(target, owned)
        paths = list(directory.iterdir())
        before = snapshot_paths(paths)
        expect_rejected(directory, matrix_args("x", force=True), f"link substitution {kind}")
        assert_snapshot(paths, before, f"link substitution {kind}")
    passed("symlink and hardlink substitutions cannot acquire destructive ownership")


def test_valid_matrix_force(base: Path) -> None:
    directory = case_dir(base, "own10_matrix_force")
    write_seed(directory)
    matrix(directory, "x")
    outputs = matrix_outputs(directory, "x")

    identical_before = snapshot_paths(outputs)
    matrix(directory, "x", force=True)
    assert_snapshot(outputs, identical_before, "OWN-10 identical matrix reuse")

    old_matrix = sha256(directory / "x.matrix_with_fuse.txt")
    (directory / "genes.nwk").write_text(
        "g((A:1.1,B:1.2):1.3,(C:1.4,D:1.5):1.6);\n",
        encoding="utf-8",
    )
    common_before = snapshot_paths(directory / name for name in COMMON_OUTPUTS)
    matrix(directory, "x", force=True)
    if sha256(directory / "x.matrix_with_fuse.txt") == old_matrix:
        fail("OWN-10 changed input did not replace the owned matrix")
    assert_snapshot(
        [directory / name for name in COMMON_OUTPUTS],
        common_before,
        "OWN-10 shared common reuse",
    )
    second_before = snapshot_paths(outputs)
    matrix(directory, "x", force=True)
    assert_snapshot(outputs, second_before, "OWN-10 second identical matrix reuse")
    passed("OWN-10 intact matrix runs replace transactionally and identical reruns do not mutate")


def test_valid_finalize_force(base: Path) -> None:
    paired = case_dir(base, "own11_paired")
    make_pair(paired)
    finalize(paired, "out", support=True)
    paired_outputs = [
        paired / "free.matrix_with_fuse.na_fuse.txt",
        paired / "fix.matrix_with_fuse.na_fuse.txt",
        paired / "out.fix.na_classified.txt",
        paired / "out.free.na_classified.txt",
        paired / "out.support_b.txt",
        paired / "species_tree.support_b.nwk",
        paired / "out.finalize_manifest.json",
    ]
    before = snapshot_paths(paired_outputs)
    finalize(paired, "out", force=True, support=True)
    assert_snapshot(paired_outputs, before, "OWN-11 paired identical reuse")
    old_manifest = sha256(paired / "out.finalize_manifest.json")
    with (paired / "free.run_manifest.json").open("a", encoding="utf-8") as handle:
        handle.write(" \n")
    finalize(paired, "out", force=True, support=True)
    if sha256(paired / "out.finalize_manifest.json") == old_manifest:
        fail("OWN-11 paired finalization manifest was not transactionally replaced")

    fixed = case_dir(base, "own11_fix_only")
    write_seed(fixed)
    matrix(fixed, "fix")
    finalize_fix(fixed, "out")
    fix_outputs = [
        fixed / "fix.matrix_with_fuse.na_fuse.txt",
        fixed / "out.fix.na_classified.txt",
        fixed / "out.finalize_manifest.json",
    ]
    before = snapshot_paths(fix_outputs)
    finalize_fix(fixed, "out", force=True)
    assert_snapshot(fix_outputs, before, "OWN-11 fix-only identical reuse")
    old_manifest = sha256(fixed / "out.finalize_manifest.json")
    with (fixed / "fix.run_manifest.json").open("a", encoding="utf-8") as handle:
        handle.write(" \n")
    finalize_fix(fixed, "out", force=True)
    if sha256(fixed / "out.finalize_manifest.json") == old_manifest:
        fail("OWN-11 fix-only finalization manifest was not transactionally replaced")
    passed("OWN-11 intact paired and fix-only finalizations support safe force replacement")


def test_shared_common_axis(base: Path) -> None:
    directory = case_dir(base, "own12_shared_axis")
    write_seed(directory)
    matrix(directory, "free")
    common = [directory / name for name in COMMON_OUTPUTS]
    common_before = snapshot_paths(common)
    matrix(directory, "fix")
    assert_snapshot(common, common_before, "OWN-12 common axis")
    matrix(directory, "free", force=True)
    matrix(directory, "fix", force=True)
    assert_snapshot(common, common_before, "OWN-12 manifests remain valid")
    passed("OWN-12 free/fix runs reuse one common axis without mutation")


def test_different_common_axis(base: Path) -> None:
    directory = case_dir(base, "own13_different_axis")
    write_seed(directory)
    matrix(directory, "x")
    (directory / "species_other.nwk").write_text(
        "((A:1,C:1):2,(B:1,D:1):3);\n", encoding="utf-8"
    )
    paths = list(directory.iterdir())
    before = snapshot_paths(paths)
    expect_rejected(
        directory,
        matrix_args("y", force=True, species="species_other.nwk"),
        "OWN-13 different axis",
    )
    assert_snapshot(paths, before, "OWN-13 different axis")
    if any(path.name.startswith("y.") or path.name.startswith("y_") for path in directory.iterdir()):
        fail("OWN-13 published label-owned outputs before rejecting the common-axis conflict")
    passed("OWN-13 a different species axis fails closed without disturbing prior runs")


def test_identical_unowned_content(base: Path) -> None:
    source = case_dir(base, "own14_source")
    write_seed(source)
    matrix(source, "x")

    shared = case_dir(base, "own14_shared")
    write_seed(shared)
    for name in COMMON_OUTPUTS:
        shutil.copy2(source / name, shared / name)
    common = [shared / name for name in COMMON_OUTPUTS]
    before = snapshot_paths(common)
    matrix(shared, "x", force=True)
    assert_snapshot(common, before, "OWN-14 identical unowned shared content")

    nonshared = case_dir(base, "own14_nonshared")
    write_seed(nonshared)
    shutil.copy2(source / "x.matrix_with_fuse.txt", nonshared / "x.matrix_with_fuse.txt")
    paths = list(nonshared.iterdir())
    before = snapshot_paths(paths)
    expect_rejected(nonshared, matrix_args("x", force=True), "OWN-14 identical unowned nonshared")
    assert_snapshot(paths, before, "OWN-14 identical unowned nonshared")
    passed("OWN-14 identical shared content is reused without mutation; nonshared content is not adopted")


def test_unicode_ownership(base: Path) -> None:
    directory = case_dir(base, "own16_\u8def\u5f84_\u6d77\u9a6c")
    species_name = "\u7269\u79cd\u6811.nwk"
    genes_name = "\u57fa\u56e0\u6811.nwk"
    label = "\u03b1_\u6d77\u9a6c"
    write_seed(directory)
    (directory / "species.nwk").rename(directory / species_name)
    (directory / "genes.nwk").rename(directory / genes_name)
    matrix(directory, label, species=species_name, genes=genes_name)
    outputs = matrix_outputs(directory, label)
    before = snapshot_paths(outputs)
    matrix(directory, label, force=True, species=species_name, genes=genes_name)
    assert_snapshot(outputs, before, "OWN-16 Unicode matrix reuse")
    ownership = load_json(directory / f"{label}.run_manifest.json")["ownership"]
    if ownership["owner_label"] != label:
        fail("OWN-16 owner label was not preserved as UTF-8")
    passed("OWN-16 Unicode CWD, input paths, labels, and ownership records remain exact")


def test_late_rollback(base: Path) -> None:
    matrix_case = case_dir(base, "own17_matrix")
    write_seed(matrix_case)
    matrix(matrix_case, "x")
    outputs = matrix_outputs(matrix_case, "x")
    before = snapshot_paths(outputs)
    (matrix_case / "genes.nwk").write_text(
        "g((A:2.1,B:2.2):2.3,(C:2.4,D:2.5):2.6);\n",
        encoding="utf-8",
    )
    result = run_sa(
        matrix_case,
        matrix_args("x", force=True),
        env={
            "SPLITALIGNER_INTERNAL_TESTING": "1",
            "SPLITALIGNER_TEST_FAIL_AFTER_PUBLISH": "2",
        },
    )
    if result.returncode == 0 or "Injected late publication failure" not in result.stderr:
        fail(f"OWN-17 matrix failpoint did not fire:\n{result.stderr}")
    assert_snapshot(outputs, before, "OWN-17 matrix rollback")
    assert_no_workspace(matrix_case)

    final_case = case_dir(base, "own17_finalize")
    make_pair(final_case)
    finalize(final_case, "out")
    final_outputs = [
        final_case / "free.matrix_with_fuse.na_fuse.txt",
        final_case / "fix.matrix_with_fuse.na_fuse.txt",
        final_case / "out.fix.na_classified.txt",
        final_case / "out.free.na_classified.txt",
        final_case / "out.finalize_manifest.json",
    ]
    before = snapshot_paths(final_outputs)
    with (final_case / "free.run_manifest.json").open("a", encoding="utf-8") as handle:
        handle.write(" \n")
    result = run_sa(
        final_case,
        finalize_args("out", force=True),
        env={
            "SPLITALIGNER_INTERNAL_TESTING": "1",
            "SPLITALIGNER_TEST_FAIL_AFTER_PUBLISH": "1",
        },
    )
    if result.returncode == 0 or "Injected late publication failure" not in result.stderr:
        fail(f"OWN-17 finalize failpoint did not fire:\n{result.stderr}")
    assert_snapshot(final_outputs, before, "OWN-17 finalize rollback")
    assert_no_workspace(final_case)
    passed("OWN-17 late publication failures restore prior owned outputs exactly")


def test_inventory_and_determinism(base: Path) -> None:
    first = case_dir(base, "inventory_first")
    second = case_dir(base, "inventory_second")
    for directory in (first, second):
        make_pair(directory)
        finalize(directory, "out", support=True)
        finalize_fix(directory, "fixed_only")

    matrix_expected = set(COMMON_OUTPUTS) | {
        "free_splits",
        "free_split_branch_label",
        "free.gene_id_map.tsv",
        "free.primitive_axis.tsv",
        "free.primitive_state.tsv",
        "free.matrix_no_fuse.txt",
        "free.matrix_with_fuse.txt",
    }
    if inventory_paths(first / "free.run_manifest.json") != matrix_expected:
        fail("matrix ownership inventory does not cover the exact publication set")
    paired_expected = {
        "free.matrix_with_fuse.na_fuse.txt",
        "fix.matrix_with_fuse.na_fuse.txt",
        "out.fix.na_classified.txt",
        "out.free.na_classified.txt",
        "out.support_b.txt",
        "species_tree.support_b.nwk",
    }
    if inventory_paths(first / "out.finalize_manifest.json") != paired_expected:
        fail("paired-finalize ownership inventory does not cover the exact publication set")
    fix_expected = {
        "fix.matrix_with_fuse.na_fuse.txt",
        "fixed_only.fix.na_classified.txt",
    }
    if inventory_paths(first / "fixed_only.finalize_manifest.json") != fix_expected:
        fail("fix-only ownership inventory does not cover the exact publication set")

    matrix_manifest = load_json(first / "free.run_manifest.json")
    ownership = matrix_manifest["ownership"]
    if ownership["schema_version"] != "SplitAligner-output-ownership-v1":
        fail("unexpected ownership schema")
    directory_entries = [
        entry for entry in ownership["outputs"] if entry["type"] == "directory"
    ]
    if len(directory_entries) != 2:
        fail("matrix ownership inventory does not include both output directories")
    for entry in directory_entries:
        fingerprint = entry["fingerprint"]
        paths = [item["relative_path"] for item in fingerprint["entries"]]
        if paths != sorted(paths) or fingerprint["entry_count"] != len(paths):
            fail("directory ownership inventory is not recursively deterministic")

    names = [
        "free.run_manifest.json",
        "fix.run_manifest.json",
        "out.finalize_manifest.json",
        "fixed_only.finalize_manifest.json",
    ]
    for name in names:
        if (first / name).read_bytes() != (second / name).read_bytes():
            fail(f"deterministic ownership manifest mismatch: {name}")
    passed("ownership inventories are complete, recursive, and byte-deterministic")


def main() -> None:
    if not SPLITALIGNER.is_file():
        fail(f"invalid source root: {ROOT}")
    with tempfile.TemporaryDirectory(prefix="splitaligner-recert006-") as temp:
        base = Path(temp)
        test_unowned_matrix_files(base)
        test_unowned_matrix_directories(base)
        test_unowned_common_files(base)
        test_unowned_finalize_outputs(base)
        test_unowned_finalize_fix_outputs(base)
        test_owner_record_defects(base)
        test_owned_directory_defects(base)
        test_link_substitutions(base)
        test_valid_matrix_force(base)
        test_valid_finalize_force(base)
        test_shared_common_axis(base)
        test_different_common_axis(base)
        test_identical_unowned_content(base)
        test_unicode_ownership(base)
        test_late_rollback(base)
        test_inventory_and_determinism(base)
    if PASS_COUNT != 16:
        fail(f"expected 16 regression groups, observed {PASS_COUNT}")
    print(f"[PASS] RECERT-006 output-ownership regression groups: {PASS_COUNT}/16")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"[FAIL] {error}", file=sys.stderr)
        raise
