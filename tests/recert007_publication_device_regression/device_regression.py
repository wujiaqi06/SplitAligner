#!/usr/bin/env python3
"""RECERT-007 controller-level publication-device regression suite."""

from __future__ import annotations

import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


ROOT = Path(sys.argv[1]).resolve()
SA = ROOT / "SplitAligner.pl"


def fail(message: str) -> None:
    raise AssertionError(message)


def passed(message: str) -> None:
    print(f"[PASS] {message}")


def run_sa(cwd: Path, args: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["perl", str(SA), *args],
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )


def assert_success(result: subprocess.CompletedProcess[str], label: str) -> None:
    if result.returncode != 0:
        fail(f"{label} failed:\n{result.stderr}")


def assert_no_workdir(path: Path) -> None:
    leftovers = [entry.name for entry in path.iterdir() if entry.name.startswith(".splitaligner-")]
    if leftovers:
        fail(f"transaction workspaces remained in {path}: {leftovers}")


def write_inputs(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)
    (path / "species.nwk").write_text("((A:1,B:1):2,(C:1,D:1):3);\n", encoding="utf-8")
    (path / "free.nwk").write_text(
        "g((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);\n", encoding="utf-8"
    )
    (path / "fix.nwk").write_text(
        "g((A:1.1,B:1.2):1.3,(C:1.4,D:1.5):1.6);\n", encoding="utf-8"
    )


def generate_matrices(path: Path) -> None:
    write_inputs(path)
    assert_success(
        run_sa(path, ["--mode", "matrix", "--species", "species.nwk", "--gene", "free.nwk", "--label", "free"]),
        "FREE matrix generation",
    )
    assert_success(
        run_sa(path, ["--mode", "matrix", "--species", "species.nwk", "--gene", "fix.nwk", "--label", "fix"]),
        "FIX matrix generation",
    )


def finalize_args(matrix_root: Path, final: str = "final", force: bool = False) -> list[str]:
    args = [
        "--mode", "finalize",
        "--free", str(matrix_root / "free.matrix_with_fuse.txt"),
        "--fix", str(matrix_root / "fix.matrix_with_fuse.txt"),
        "--final_label", final,
    ]
    if force:
        args.append("--force")
    return args


def assert_early_device_rejection(
    cwd: Path,
    args: list[str],
    label: str,
    expected_role: str,
) -> None:
    before = {str(path): path.read_bytes() for path in cwd.rglob("*") if path.is_file()}
    result = run_sa(cwd, args)
    if result.returncode == 0:
        fail(f"{label} unexpectedly succeeded")
    required = [
        "[ERROR][Publication device]",
        expected_role,
        "Transaction parent",
        "device",
        "--force cannot override",
    ]
    for text in required:
        if text not in result.stderr:
            fail(f"{label} diagnostic omitted {text!r}:\n{result.stderr}")
    if "Mark NA_fuse" in result.stderr or "Invalid cross-device link" in result.stderr:
        fail(f"{label} was not rejected during early device preflight:\n{result.stderr}")
    after = {str(path): path.read_bytes() for path in cwd.rglob("*") if path.is_file()}
    if before != after:
        fail(f"{label} changed files in the transaction directory")
    assert_no_workdir(cwd)
    passed(label)


def candidate_roots() -> list[Path]:
    roots: list[Path] = []
    explicit = os.environ.get("SPLITALIGNER_SECOND_DEVICE_ROOT")
    if explicit:
        roots.append(Path(explicit))
    roots.extend(Path(path) for path in ("/dev/shm", "/run/shm", "/var/run/shm"))
    volumes = Path("/Volumes")
    if volumes.is_dir():
        roots.extend(sorted(volumes.iterdir(), key=lambda path: path.name))
    return roots


def find_second_device(local_device: int) -> Path | None:
    seen: set[Path] = set()
    for root in candidate_roots():
        try:
            resolved = root.resolve(strict=True)
            if resolved in seen or resolved.stat().st_dev == local_device:
                continue
            seen.add(resolved)
            probe = Path(tempfile.mkdtemp(prefix="splitaligner-recert007-probe-", dir=resolved))
            probe.rmdir()
            return resolved
        except (OSError, RuntimeError):
            continue
    return None


def same_device_checks(local: Path) -> None:
    matrix_root = local / "same_device_matrices"
    output_root = local / "same_device_outputs"
    output_root.mkdir(parents=True)
    generate_matrices(matrix_root)
    assert_success(run_sa(output_root, finalize_args(matrix_root)), "same-device absolute finalize")
    if not (output_root / "final.finalize_manifest.json").is_file():
        fail("same-device absolute finalize did not publish its manifest")
    assert_no_workdir(output_root)
    passed("DEV-07 same-device absolute matrix paths")

    external_inputs = local / "external_inputs"
    external_output = local / "external_matrix_output"
    write_inputs(external_inputs)
    external_output.mkdir()
    assert_success(
        run_sa(
            external_output,
            [
                "--mode", "matrix",
                "--species", str(external_inputs / "species.nwk"),
                "--gene", str(external_inputs / "free.nwk"),
                "--label", "external",
            ],
        ),
        "matrix with external same-device inputs",
    )
    passed("matrix inputs are not publication destinations")


def cross_device_checks(local: Path, remote_base: Path) -> None:
    remote = Path(tempfile.mkdtemp(prefix="splitaligner-recert007-remote-", dir=remote_base))
    try:
        matrix_root = remote / "矩阵_输入"
        generate_matrices(matrix_root)

        paired = local / "paired"
        paired.mkdir()
        assert_early_device_rejection(
            paired, finalize_args(matrix_root), "DEV-01 paired cross-device", "FREE matrix NA_fuse output"
        )

        fix_only = local / "fix_only"
        fix_only.mkdir()
        assert_early_device_rejection(
            fix_only,
            [
                "--mode", "finalize_fix",
                "--fix", str(matrix_root / "fix.matrix_with_fuse.txt"),
                "--final_label", "final",
            ],
            "DEV-02 fix-only cross-device",
            "FIX matrix NA_fuse output",
        )

        local_fix_root = local / "mixed_local_fix"
        generate_matrices(local_fix_root)
        mixed = local / "mixed"
        mixed.mkdir()
        mixed_args = [
            "--mode", "finalize",
            "--free", str(matrix_root / "free.matrix_with_fuse.txt"),
            "--fix", str(local_fix_root / "fix.matrix_with_fuse.txt"),
            "--final_label", "final",
        ]
        assert_early_device_rejection(
            mixed, mixed_args, "DEV-03 mixed publication devices", "FREE matrix NA_fuse output"
        )

        unicode_cwd = local / "事务_目录"
        unicode_cwd.mkdir()
        assert_early_device_rejection(
            unicode_cwd,
            finalize_args(matrix_root),
            "DEV-04 Unicode cross-device diagnostic",
            "矩阵_输入",
        )

        symlink_root = local / "remote_matrix_link"
        symlink_root.symlink_to(matrix_root, target_is_directory=True)
        symlink_case = local / "symlink_case"
        symlink_case.mkdir()
        assert_early_device_rejection(
            symlink_case,
            finalize_args(symlink_root),
            "DEV-05 symlink-resolved matrix device",
            "FREE matrix NA_fuse output",
        )

        forced = local / "forced"
        forced.mkdir()
        assert_early_device_rejection(
            forced,
            finalize_args(matrix_root, force=True),
            "DEV-06 --force cannot override",
            "FREE matrix NA_fuse output",
        )

        external = remote / "external_matrix_inputs"
        write_inputs(external)
        external_output = local / "external_cross_device_inputs"
        external_output.mkdir()
        assert_success(
            run_sa(
                external_output,
                [
                    "--mode", "matrix",
                    "--species", str(external / "species.nwk"),
                    "--gene", str(external / "free.nwk"),
                    "--label", "external",
                ],
            ),
            "DEV-08 matrix mode with cross-device read-only inputs",
        )
        passed("DEV-08 external species/gene inputs remain supported")

        for forbidden in (
            matrix_root / "free.matrix_with_fuse.na_fuse.txt",
            matrix_root / "fix.matrix_with_fuse.na_fuse.txt",
        ):
            if forbidden.exists():
                fail(f"cross-device rejection published {forbidden}")
    finally:
        shutil.rmtree(remote)


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="splitaligner-recert007-") as temp:
        local = Path(temp)
        same_device_checks(local)
        second = find_second_device(local.stat().st_dev)
        if second is None:
            print("[SKIP] RECERT-007 real cross-device integration: no writable second filesystem device.")
            print("[INFO] Synthetic device comparisons and same-device controller checks completed.")
            return 0
        print(f"[INFO] RECERT-007 second filesystem root: {second}")
        cross_device_checks(local, second)
        print("[PASS] RECERT-007 real cross-device integration suite")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
