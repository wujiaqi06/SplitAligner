#!/usr/bin/env python3
"""Controller regressions for RECERT-008 publication-parent device authority."""

from __future__ import annotations

from pathlib import Path
import subprocess
import sys
import tempfile


ROOT = Path(sys.argv[1]).resolve()
SA = ROOT / "SplitAligner.pl"


def run_sa(cwd: Path, args: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["perl", str(SA), *args], cwd=cwd, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
    )


def require_success(result: subprocess.CompletedProcess[str], label: str) -> None:
    if result.returncode != 0:
        raise AssertionError(f"{label} failed:\n{result.stderr}")
    print(f"[PASS] {label}")


def write_inputs(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)
    (path / "species.nwk").write_text(
        "((A:1,B:1):2,(C:1,D:1):3);\n", encoding="utf-8")
    (path / "free.nwk").write_text(
        "g((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);\n", encoding="utf-8")
    (path / "fix.nwk").write_text(
        "g((A:1.1,B:1.2):1.3,(C:1.4,D:1.5):1.6);\n", encoding="utf-8")


def matrix(path: Path, label: str, gene: str, force: bool = False) -> None:
    args = [
        "--mode", "matrix", "--species", "species.nwk",
        "--gene", gene, "--label", label,
    ]
    if force:
        args.append("--force")
    require_success(run_sa(path, args), f"matrix {label}{' --force' if force else ''}")


def paired_args(matrix_root: Path, label: str, force: bool = False) -> list[str]:
    args = [
        "--mode", "finalize",
        "--free", str(matrix_root / "free.matrix_with_fuse.txt"),
        "--fix", str(matrix_root / "fix.matrix_with_fuse.txt"),
        "--final_label", label,
    ]
    if force:
        args.append("--force")
    return args


def fix_only_args(matrix_root: Path, label: str, force: bool = False) -> list[str]:
    args = [
        "--mode", "finalize_fix",
        "--fix", str(matrix_root / "fix.matrix_with_fuse.txt"),
        "--final_label", label,
    ]
    if force:
        args.append("--force")
    return args


def assert_no_workdir(path: Path) -> None:
    leftovers = [p.name for p in path.iterdir() if p.name.startswith(".splitaligner-")]
    if leftovers:
        raise AssertionError(f"temporary workdirs remained: {leftovers}")


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="splitaligner-recert008-controller-") as tmp:
        base = Path(tmp)

        two_label = base / "two_label_same_cwd"
        write_inputs(two_label)
        matrix(two_label, "free", "free.nwk")
        matrix(two_label, "fix", "fix.nwk")
        for name in ("free.run_manifest.json", "fix.run_manifest.json"):
            if not (two_label / name).is_file():
                raise AssertionError(f"PARENT-01 missing {name}")
        assert_no_workdir(two_label)
        print("[PASS] PARENT-01 two-label same-CWD matrix workflow")

        matrix_force = base / "matrix_force"
        write_inputs(matrix_force)
        matrix(matrix_force, "free", "free.nwk")
        matrix(matrix_force, "free", "free.nwk", force=True)
        assert_no_workdir(matrix_force)
        print("[PASS] PARENT-04 ownership-verified matrix --force rerun")

        paired_root = base / "paired_outputs"
        paired_root.mkdir()
        require_success(run_sa(paired_root, paired_args(two_label, "paired")),
                        "paired finalize initial")
        require_success(run_sa(paired_root, paired_args(two_label, "paired", force=True)),
                        "paired finalize --force rerun")
        assert_no_workdir(paired_root)

        fix_root = base / "fix_only_outputs"
        fix_root.mkdir()
        require_success(run_sa(fix_root, fix_only_args(two_label, "fix_only")),
                        "fix-only finalize initial")
        require_success(run_sa(fix_root, fix_only_args(two_label, "fix_only", force=True)),
                        "fix-only finalize --force rerun")
        assert_no_workdir(fix_root)
        print("[PASS] PARENT-05 paired and fix-only ownership-verified force reruns")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
