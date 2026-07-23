#!/usr/bin/env python3
"""Regression coverage for non-destructive bundled-example reset."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import stat
import subprocess
import tempfile
import time


SOURCE_ROOT = Path(__file__).resolve().parents[2]
ERROR_PREFIX = "[ERROR][Example archive]"
FIXED_TIMESTAMP = "20260721T120000Z"


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
        mode = oct(stat.S_IMODE(metadata.st_mode))
        if path.is_symlink():
            items.append([relative, "symlink", os.readlink(path), mode])
        elif path.is_dir():
            items.append([relative, "directory", mode])
        elif path.is_file():
            items.append([relative, "file", digest(path), mode])
        else:
            items.append([relative, "special", str(stat.S_IFMT(metadata.st_mode)), mode])
    return json.dumps(items, ensure_ascii=False, separators=(",", ":"))


def copy_checkout(destination: Path) -> Path:
    shutil.copytree(SOURCE_ROOT, destination, symlinks=True)
    return destination


def minimal_checkout(destination: Path) -> Path:
    (destination / "examples/302mammal").mkdir(parents=True)
    (destination / "examples/preprint_302mammal").mkdir(parents=True)
    shutil.copy2(SOURCE_ROOT / "examples/archive_run.sh", destination / "examples/archive_run.sh")
    return destination


def example_root(checkout: Path, mode: str) -> Path:
    name = "302mammal" if mode == "toy" else "preprint_302mammal"
    return checkout / "examples" / name


def populate_run(checkout: Path, mode: str) -> Path:
    run = example_root(checkout, mode) / "run"
    (run / "nested").mkdir(parents=True)
    (run / "important.txt").write_text("IMPORTANT USER CONTENT\n", encoding="utf-8")
    (run / "final.support_b.txt").write_text("branch_id\tsupport\nB1\t100\n", encoding="utf-8")
    (run / "nested/payload.bin").write_bytes(bytes(range(256)))
    os.chmod(run / "nested/payload.bin", 0o640)
    return run


def archive_command(checkout: Path, mode: str, *, env: dict[str, str] | None = None) -> subprocess.CompletedProcess[str]:
    merged = os.environ.copy()
    merged.update({"LC_ALL": "C", "LANG": "C", "PYTHONDONTWRITEBYTECODE": "1"})
    if env:
        merged.update(env)
    return subprocess.run(
        ["bash", "examples/archive_run.sh", mode],
        cwd=checkout,
        env=merged,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def run_toy(checkout: Path) -> subprocess.CompletedProcess[str]:
    env = os.environ.copy()
    env.update({"LC_ALL": "C", "LANG": "C", "PYTHONDONTWRITEBYTECODE": "1"})
    return subprocess.run(
        ["bash", "examples/run.sh", "toy"],
        cwd=checkout,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def archives(root: Path) -> list[Path]:
    return sorted(path for path in root.iterdir() if path.name.startswith("run.saved.") and path.is_dir())


def require_success(result: subprocess.CompletedProcess[str], label: str) -> None:
    if result.returncode:
        fail(f"{label} failed:\n{result.stdout}\n{result.stderr}")


def require_rejection(result: subprocess.CompletedProcess[str], label: str) -> None:
    if result.returncode == 0:
        fail(f"{label} unexpectedly succeeded")
    if ERROR_PREFIX not in result.stderr:
        fail(f"{label} lacks stable archive diagnostic:\n{result.stdout}\n{result.stderr}")


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


def test_basic_archives(base: Path) -> None:
    for mode in ("toy", "preprint"):
        checkout = minimal_checkout(base / f"basic-{mode}")
        run = populate_run(checkout, mode)
        before = snapshot(run)
        result = archive_command(checkout, mode)
        require_success(result, f"basic {mode} archive")
        if run.exists() or run.is_symlink():
            fail(f"{mode} source workspace still exists")
        saved = archives(example_root(checkout, mode))
        if len(saved) != 1 or snapshot(saved[0]) != before:
            fail(f"{mode} archive did not preserve every workspace byte and mode")
        if "inspect this archive" not in result.stdout:
            fail(f"{mode} archive output lacks inspection warning")


def test_absent_noop(base: Path) -> None:
    checkout = minimal_checkout(base / "absent")
    root = example_root(checkout, "toy")
    before = snapshot(root)
    result = archive_command(checkout, "toy")
    require_success(result, "absent run no-op")
    if snapshot(root) != before or (root / "run").exists() or archives(root):
        fail("absent run no-op created or changed content")


def test_invalid_entries(base: Path) -> None:
    for kind in ("symlink", "file", "fifo"):
        if kind == "fifo" and not hasattr(os, "mkfifo"):
            continue
        checkout = minimal_checkout(base / f"invalid-{kind}")
        root = example_root(checkout, "toy")
        run = root / "run"
        if kind == "symlink":
            external = base / "invalid-symlink-external"
            external.mkdir()
            (external / "important.txt").write_text("external\n", encoding="utf-8")
            before = snapshot(external)
            run.symlink_to(external, target_is_directory=True)
        elif kind == "file":
            run.write_text("not a directory\n", encoding="utf-8")
            before = digest(run)
        else:
            os.mkfifo(run)
            before = stat.S_IFMT(run.lstat().st_mode)
        result = archive_command(checkout, "toy")
        require_rejection(result, f"run {kind}")
        if kind == "symlink" and (not run.is_symlink() or snapshot(external) != before):
            fail("symlink rejection changed the link or its target")
        if kind == "file" and (not run.is_file() or digest(run) != before):
            fail("file rejection changed the run entry")
        if kind == "fifo" and stat.S_IFMT(run.lstat().st_mode) != before:
            fail("FIFO rejection changed the run entry")
        if archives(root):
            fail(f"run {kind} rejection created an archive")


def test_unicode_and_symlink_checkout(base: Path) -> None:
    physical = minimal_checkout(base / "親ディレクトリ" / "源コード Δ")
    run = populate_run(physical, "toy")
    before = snapshot(run)
    linked = base / "checkout-root-link"
    linked.symlink_to(physical, target_is_directory=True)
    result = archive_command(linked, "toy")
    require_success(result, "Unicode checkout through symlink ancestor")
    saved = archives(example_root(physical, "toy"))
    if len(saved) != 1 or snapshot(saved[0]) != before:
        fail("Unicode/symlink checkout archive did not preserve content")


def write_executable(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")
    path.chmod(0o755)


def wait_for(path: Path, timeout: float = 5.0) -> None:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if path.exists():
            return
        time.sleep(0.01)
    fail(f"timed out waiting for test controller marker: {path}")


def test_collision(base: Path) -> None:
    checkout = minimal_checkout(base / "collision")
    run = populate_run(checkout, "toy")
    before = snapshot(run)
    mock = base / "collision-mock"
    mock.mkdir()
    ready = mock / "date.ready"
    release = mock / "date.release"
    write_executable(
        mock / "date",
        "#!/bin/sh\n"
        'touch "$SPLITALIGNER_DATE_READY"\n'
        'while [ ! -e "$SPLITALIGNER_DATE_RELEASE" ]; do /bin/sleep 0.01; done\n'
        f"printf '%s\\n' '{FIXED_TIMESTAMP}'\n",
    )
    env = os.environ.copy()
    env.update(
        {
            "LC_ALL": "C",
            "LANG": "C",
            "PATH": f"{mock}:{env['PATH']}",
            "SPLITALIGNER_DATE_READY": str(ready),
            "SPLITALIGNER_DATE_RELEASE": str(release),
        }
    )
    process = subprocess.Popen(
        ["bash", "examples/archive_run.sh", "toy"],
        cwd=checkout,
        env=env,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    wait_for(ready)
    root = example_root(checkout, "toy")
    base_name = f"run.saved.{FIXED_TIMESTAMP}.{process.pid}"
    for suffix in ("", ".1"):
        occupied = root / f"{base_name}{suffix}"
        occupied.mkdir()
        (occupied / "preexisting.txt").write_text("do not replace\n", encoding="utf-8")
    release.touch()
    stdout, stderr = process.communicate(timeout=15)
    if process.returncode:
        fail(f"collision archive failed:\n{stdout}\n{stderr}")
    expected = root / f"{base_name}.2"
    if not expected.is_dir() or snapshot(expected) != before:
        fail("collision handling did not choose the next absent archive name")
    for suffix in ("", ".1"):
        occupied = root / f"{base_name}{suffix}/preexisting.txt"
        if occupied.read_text(encoding="utf-8") != "do not replace\n":
            fail("collision handling replaced a pre-existing archive")


def test_failed_rename(base: Path) -> None:
    checkout = minimal_checkout(base / "rename-failure")
    run = populate_run(checkout, "toy")
    before = snapshot(run)
    mock = base / "rename-failure-mock"
    mock.mkdir()
    real_perl = shutil.which("perl")
    if not real_perl:
        fail("Perl is unavailable for rename-failure test")
    write_executable(
        mock / "perl",
        "#!/bin/sh\n"
        "case \"$2\" in\n"
        "  *'atomic archive rename failed for'*) exit 71 ;;\n"
        "esac\n"
        f'exec "{real_perl}" "$@"\n',
    )
    result = archive_command(checkout, "toy", env={"PATH": f"{mock}:{os.environ['PATH']}"})
    require_rejection(result, "simulated rename failure")
    if not run.is_dir() or snapshot(run) != before:
        fail("failed rename changed the original workspace")
    if archives(example_root(checkout, "toy")):
        fail("failed rename created an archive")


def test_cross_device_preflight(base: Path) -> None:
    checkout = minimal_checkout(base / "cross-device")
    run = populate_run(checkout, "toy")
    before = snapshot(run)
    mock = base / "cross-device-mock"
    mock.mkdir()
    real_perl = shutil.which("perl")
    if not real_perl:
        fail("Perl is unavailable for cross-device test")
    wrapper = f'''#!/usr/bin/env python3
import os
import subprocess
import sys
real = {real_perl!r}
if len(sys.argv) >= 4 and 'cannot lstat' in sys.argv[2] and sys.argv[-1].endswith('/run'):
    result = subprocess.run([real, *sys.argv[1:]], text=True, stdout=subprocess.PIPE)
    if result.returncode:
        raise SystemExit(result.returncode)
    print('999999999:' + result.stdout.strip().split(':', 1)[1])
    raise SystemExit(0)
os.execv(real, [real, *sys.argv[1:]])
'''
    write_executable(mock / "perl", wrapper)
    result = archive_command(checkout, "toy", env={"PATH": f"{mock}:{os.environ['PATH']}"})
    require_rejection(result, "simulated cross-device/mounted workspace")
    if not run.is_dir() or snapshot(run) != before:
        fail("cross-device preflight changed the original workspace")
    if archives(example_root(checkout, "toy")):
        fail("cross-device preflight created an archive")


def test_literal_toy_archive_and_rerun(base: Path) -> None:
    checkout = copy_checkout(base / "literal-toy")
    first = run_toy(checkout)
    require_success(first, "first literal toy run")
    assert_expected_outputs(checkout)
    run = checkout / "examples/302mammal/run"
    (run / "important.txt").write_text("IMPORTANT USER CONTENT\n", encoding="utf-8")
    before = snapshot(run)
    archived = archive_command(checkout, "toy")
    require_success(archived, "literal toy archive")
    saved = archives(checkout / "examples/302mammal")
    if len(saved) != 1 or snapshot(saved[0]) != before:
        fail("literal toy archive did not preserve generated and unrelated content")
    second = run_toy(checkout)
    require_success(second, "post-archive literal toy rerun")
    assert_expected_outputs(checkout)
    if (saved[0] / "important.txt").read_text(encoding="utf-8") != "IMPORTANT USER CONTENT\n":
        fail("post-archive rerun changed preserved unrelated content")


def test_documentation() -> None:
    documents = [SOURCE_ROOT / "README.md", SOURCE_ROOT / "examples/README.md"]
    combined = "\n".join(path.read_text(encoding="utf-8") for path in documents)
    unsafe = re.compile(r"rm\s+-r[fF]\s+examples/(?:302mammal|preprint_302mammal)/run")
    if unsafe.search(combined):
        fail("documentation still contains recursive example-workspace deletion")
    for command in ("bash examples/archive_run.sh toy", "bash examples/archive_run.sh preprint"):
        if command not in combined:
            fail(f"documentation omits safe archive command: {command}")


def main() -> int:
    print("[INFO] RECERT-013 regression: safe example archive/reset", flush=True)
    test_documentation()
    with tempfile.TemporaryDirectory(prefix="splitaligner-recert013-") as tmp:
        base = Path(tmp)
        test_basic_archives(base)
        test_absent_noop(base)
        test_invalid_entries(base)
        test_unicode_and_symlink_checkout(base)
        test_collision(base)
        test_failed_rename(base)
        test_cross_device_preflight(base)
        test_literal_toy_archive_and_rerun(base)
    print("[PASS] RECERT-013 safe example archive/reset regression")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
