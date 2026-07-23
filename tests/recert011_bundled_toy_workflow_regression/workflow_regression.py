#!/usr/bin/env python3
"""Regression coverage for ownership-compatible bundled example workflows."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import tempfile


SOURCE_ROOT = Path(__file__).resolve().parents[2]


def fail(message: str) -> None:
    raise SystemExit(f"[FAIL] {message}")


def run_example(root: Path, mode: str, extra_env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    env.update({"LC_ALL": "C", "LANG": "C", "PYTHONDONTWRITEBYTECODE": "1"})
    if extra_env:
        env.update(extra_env)
    return subprocess.run(
        ["bash", "examples/run.sh", mode],
        cwd=root,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def snapshot_files(root: Path, excluded: tuple[Path, ...] = ()) -> dict[str, str]:
    excluded_resolved = tuple(path.resolve() for path in excluded)
    result: dict[str, str] = {}
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        resolved = path.resolve()
        if any(resolved == item or item in resolved.parents for item in excluded_resolved):
            continue
        result[path.relative_to(root).as_posix()] = digest(path)
    return result


def compare_expected(root: Path) -> None:
    expected_dir = root / "examples/302mammal/expected"
    run_dir = root / "examples/302mammal/run"
    expected = sorted(path for path in expected_dir.iterdir() if path.is_file())
    if len(expected) != 12:
        fail(f"expected fixture count is {len(expected)}, not 12")
    for reference in expected:
        actual = run_dir / reference.name
        if not actual.is_file() or digest(actual) != digest(reference):
            fail(f"toy output differs from expected: {reference.name}")


def require_manifests(run_dir: Path) -> None:
    for name in ("free.run_manifest.json", "fix.run_manifest.json", "final.finalize_manifest.json"):
        if not (run_dir / name).is_file():
            fail(f"owned run manifest missing: {name}")


def assert_success(result: subprocess.CompletedProcess[str], label: str) -> None:
    if result.returncode:
        fail(f"{label} exited {result.returncode}:\n{result.stdout}\n{result.stderr}")


def run_example_logged(root: Path, mode: str, stdout_path: Path, stderr_path: Path) -> int:
    env = os.environ.copy()
    env.update({"LC_ALL": "C", "LANG": "C", "PYTHONDONTWRITEBYTECODE": "1"})
    with stdout_path.open("w", encoding="utf-8") as stdout_handle, stderr_path.open(
        "w", encoding="utf-8"
    ) as stderr_handle:
        result = subprocess.run(
            ["bash", "examples/run.sh", mode],
            cwd=root,
            env=env,
            text=True,
            stdout=stdout_handle,
            stderr=stderr_handle,
        )
    return result.returncode


def main() -> int:
    print("[INFO] RECERT-011 regression: verify pristine source", flush=True)
    if (SOURCE_ROOT / "examples/302mammal/run").exists() or (SOURCE_ROOT / "examples/preprint_302mammal/run").exists():
        fail("source tree contains a generated example run directory")

    with tempfile.TemporaryDirectory(prefix="splitaligner-recert011-") as tmp:
        base = Path(tmp)
        checkout = base / "fresh-checkout"
        print("[INFO] RECERT-011 regression: copy fresh checkout", flush=True)
        shutil.copytree(SOURCE_ROOT, checkout)

        toy_run = checkout / "examples/302mammal/run"
        preprint_run = checkout / "examples/preprint_302mammal/run"
        if toy_run.exists():
            fail("fresh checkout unexpectedly contains toy run/")

        static_before = snapshot_files(
            checkout,
            excluded=(toy_run, preprint_run),
        )
        fixture_before = snapshot_files(checkout / "examples/302mammal/input") | snapshot_files(
            checkout / "examples/302mammal/expected"
        )

        print("[INFO] RECERT-011 regression: first literal toy run", flush=True)
        first = run_example(checkout, "toy")
        assert_success(first, "first literal toy run")
        require_manifests(toy_run)
        compare_expected(checkout)

        print("[INFO] RECERT-011 regression: second literal toy run", flush=True)
        second = run_example(checkout, "toy")
        assert_success(second, "second literal toy run")
        require_manifests(toy_run)
        compare_expected(checkout)

        fixture_after = snapshot_files(checkout / "examples/302mammal/input") | snapshot_files(
            checkout / "examples/302mammal/expected"
        )
        if fixture_after != fixture_before:
            fail("toy input/expected fixtures changed")
        if snapshot_files(checkout, excluded=(toy_run, preprint_run)) != static_before:
            fail("toy workflow wrote outside its dedicated run directory")

        owned_before_failure = snapshot_files(toy_run)
        fixtures_before_failure = fixture_after.copy()
        alternate_gene = base / "alternate_free.nwk"
        original_gene = (checkout / "examples/302mammal/input/free_tree.examples.nwk").read_text(encoding="utf-8")
        alternate_gene.write_text(
            original_gene.replace("0.0041071867", "0.5041071867", 1),
            encoding="utf-8",
        )
        failure_env = os.environ.copy()
        failure_env.update(
            {
                "LC_ALL": "C",
                "LANG": "C",
                "PYTHONDONTWRITEBYTECODE": "1",
                "SPLITALIGNER_INTERNAL_TESTING": "1",
                "SPLITALIGNER_TEST_FAIL_AFTER_PUBLISH": "2",
            }
        )
        print("[INFO] RECERT-011 regression: injected late-failure rollback", flush=True)
        failed = subprocess.run(
            [
                "perl",
                str(checkout / "SplitAligner.pl"),
                "--mode",
                "matrix",
                "--species",
                str(checkout / "examples/302mammal/input/speciesTree302.nwk"),
                "--gene",
                str(alternate_gene),
                "--label",
                "free",
                "--force",
            ],
            cwd=toy_run,
            env=failure_env,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        if failed.returncode == 0:
            fail("late-failure toy invocation unexpectedly succeeded")
        if snapshot_files(toy_run) != owned_before_failure:
            fail("late failure changed prior owned toy outputs")
        if snapshot_files(checkout / "examples/302mammal/input") | snapshot_files(
            checkout / "examples/302mammal/expected"
        ) != fixtures_before_failure:
            fail("late failure changed immutable toy fixtures")

        preprint_input_before = snapshot_files(checkout / "examples/preprint_302mammal/input")
        print("[INFO] RECERT-011 regression: literal preprint run", flush=True)
        preprint_stdout = base / "preprint.stdout.txt"
        preprint_stderr = base / "preprint.stderr.txt"
        preprint_rc = run_example_logged(checkout, "preprint", preprint_stdout, preprint_stderr)
        if preprint_rc:
            fail(
                f"literal preprint run exited {preprint_rc}:\n"
                + preprint_stdout.read_text(encoding="utf-8")[-4000:]
                + preprint_stderr.read_text(encoding="utf-8")[-4000:]
            )
        require_manifests(preprint_run)
        if snapshot_files(checkout / "examples/preprint_302mammal/input") != preprint_input_before:
            fail("preprint input fixtures changed")
        if snapshot_files(checkout, excluded=(toy_run, preprint_run)) != static_before:
            fail("preprint workflow wrote outside its dedicated run directory")

        print("[INFO] RECERT-011 regression: Unicode checkout toy run", flush=True)
        unicode_root = base / "親ディレクトリ" / "源コード Δ"
        unicode_root.parent.mkdir(parents=True)
        shutil.copytree(SOURCE_ROOT, unicode_root)
        unicode_run = run_example(unicode_root, "toy")
        assert_success(unicode_run, "Unicode checkout literal toy run")
        require_manifests(unicode_root / "examples/302mammal/run")
        compare_expected(unicode_root)

    print("[PASS] RECERT-011 bundled toy and preprint workflow regression")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
