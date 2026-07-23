#!/usr/bin/env python3
"""Regression coverage for bundled example workspace identity."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import shutil
import stat
import subprocess
import tempfile


SOURCE_ROOT = Path(__file__).resolve().parents[2]
ERROR_PREFIX = "[ERROR][Example workspace]"
CONTROLLER_MARKER = "[INFO] Preflighting species"


def fail(message: str) -> None:
    raise SystemExit(f"[FAIL] {message}")


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def snapshot(root: Path) -> str:
    items: list[list[str]] = []
    for path in sorted(root.rglob("*")):
        relative = path.relative_to(root).as_posix()
        metadata = path.lstat()
        if path.is_symlink():
            items.append([relative, "symlink", os.readlink(path)])
        elif path.is_dir():
            items.append([relative, "directory", oct(stat.S_IMODE(metadata.st_mode))])
        elif path.is_file():
            items.append([relative, "file", digest(path)])
        else:
            items.append([relative, "special", str(metadata.st_mode)])
    return json.dumps(items, ensure_ascii=False, separators=(",", ":"))


def copy_checkout(destination: Path) -> Path:
    shutil.copytree(SOURCE_ROOT, destination, symlinks=True)
    return destination


def run_example(root: Path, mode: str, *, via_shell_cd: bool = False) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    env.update({"LC_ALL": "C", "LANG": "C", "PYTHONDONTWRITEBYTECODE": "1"})
    if via_shell_cd:
        command = ["bash", "-c", 'cd "$1" && bash examples/run.sh "$2"', "bash", str(root), mode]
        cwd = root.parent
    else:
        command = ["bash", "examples/run.sh", mode]
        cwd = root
    return subprocess.run(
        command,
        cwd=cwd,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def assert_early_rejection(result: subprocess.CompletedProcess[str], label: str) -> None:
    combined = result.stdout + result.stderr
    if result.returncode == 0:
        fail(f"{label} unexpectedly succeeded")
    if ERROR_PREFIX not in result.stderr:
        fail(f"{label} lacks stable workspace diagnostic:\n{combined}")
    if CONTROLLER_MARKER in combined or "run_manifest.json" in combined:
        fail(f"{label} reached SplitAligner before rejection:\n{combined}")


def assert_expected_outputs(checkout: Path) -> None:
    expected_dir = checkout / "examples/302mammal/expected"
    run_dir = checkout / "examples/302mammal/run"
    references = sorted(path for path in expected_dir.iterdir() if path.is_file())
    if len(references) != 12:
        fail(f"expected fixture count is {len(references)}, not 12")
    for reference in references:
        actual = run_dir / reference.name
        if not actual.is_file() or digest(actual) != digest(reference):
            fail(f"toy output differs from expected: {reference.name}")


def reject_symlink_case(base: Path, name: str, mode: str, target_kind: str) -> None:
    checkout = copy_checkout(base / name)
    example_name = "302mammal" if mode == "toy" else "preprint_302mammal"
    example_root = checkout / "examples" / example_name
    run_entry = example_root / "run"
    external = base / f"{name}-external"
    if target_kind == "external":
        external.mkdir()
        target = external
        target_snapshot = snapshot(external)
        run_entry.symlink_to(target)
    else:
        target = example_root / target_kind
        target_snapshot = snapshot(target)
        run_entry.symlink_to(target_kind)

    result = run_example(checkout, mode)
    assert_early_rejection(result, f"{mode} run -> {target_kind} symlink")
    if snapshot(target) != target_snapshot:
        fail(f"{mode} run -> {target_kind} changed target inventory")


def reject_non_directory_case(base: Path, kind: str) -> None:
    checkout = copy_checkout(base / f"non-directory-{kind}")
    run_entry = checkout / "examples/302mammal/run"
    if kind == "file":
        run_entry.write_text("not a workspace\n", encoding="utf-8")
    elif kind == "fifo":
        os.mkfifo(run_entry)
    else:
        fail(f"unknown non-directory kind: {kind}")
    before = run_entry.lstat()
    result = run_example(checkout, "toy")
    assert_early_rejection(result, f"run is {kind}")
    after = run_entry.lstat()
    if stat.S_IFMT(before.st_mode) != stat.S_IFMT(after.st_mode):
        fail(f"run {kind} changed object type")


def valid_workflow_cases(base: Path) -> None:
    checkout = copy_checkout(base / "fresh-real-run")
    input_dir = checkout / "examples/302mammal/input"
    expected_dir = checkout / "examples/302mammal/expected"
    immutable_before = (snapshot(input_dir), snapshot(expected_dir))

    for ordinal in ("first", "second"):
        result = run_example(checkout, "toy")
        if result.returncode:
            fail(f"{ordinal} fresh toy run failed:\n{result.stdout}\n{result.stderr}")
        assert_expected_outputs(checkout)
    if (snapshot(input_dir), snapshot(expected_dir)) != immutable_before:
        fail("valid toy run changed immutable input/expected inventory")

    unicode_checkout = copy_checkout(base / "親ディレクトリ" / "源コード Δ")
    unicode_result = run_example(unicode_checkout, "toy")
    if unicode_result.returncode:
        fail(f"Unicode checkout failed:\n{unicode_result.stdout}\n{unicode_result.stderr}")
    assert_expected_outputs(unicode_checkout)

    symlink_checkout = base / "checkout-root-link"
    symlink_checkout.symlink_to(unicode_checkout, target_is_directory=True)
    shutil.rmtree(unicode_checkout / "examples/302mammal/run")
    linked_result = run_example(symlink_checkout, "toy", via_shell_cd=True)
    if linked_result.returncode:
        fail(f"checkout-root symlink failed:\n{linked_result.stdout}\n{linked_result.stderr}")
    assert_expected_outputs(unicode_checkout)


def main() -> int:
    print("[INFO] RECERT-012 regression: verify pristine source", flush=True)
    for relative in ("examples/302mammal/run", "examples/preprint_302mammal/run"):
        if (SOURCE_ROOT / relative).exists() or (SOURCE_ROOT / relative).is_symlink():
            fail(f"source tree contains generated workspace: {relative}")

    with tempfile.TemporaryDirectory(prefix="splitaligner-recert012-") as tmp:
        base = Path(tmp)
        reject_symlink_case(base, "toy-run-to-input", "toy", "input")
        reject_symlink_case(base, "toy-run-to-expected", "toy", "expected")
        reject_symlink_case(base, "toy-run-to-external", "toy", "external")
        reject_symlink_case(base, "preprint-run-to-input", "preprint", "input")
        reject_non_directory_case(base, "file")
        if hasattr(os, "mkfifo"):
            reject_non_directory_case(base, "fifo")
        valid_workflow_cases(base)

    print("[PASS] RECERT-012 example workspace identity regression")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
