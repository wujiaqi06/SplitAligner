#!/usr/bin/env python3
"""Controller-level regressions for RECERT-010 strict POSIX input lookup."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path
import shutil
import stat
import subprocess
import sys
import tempfile


ROOT = Path(sys.argv[1]).resolve()
SA = ROOT / "SplitAligner.pl"


def run_sa(cwd: Path, args: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["perl", str(SA), *args],
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )


def require_success(result: subprocess.CompletedProcess[str], label: str) -> None:
    if result.returncode != 0:
        raise AssertionError(f"{label} failed:\n{result.stderr}")
    print(f"[PASS] {label}")


def tree_snapshot(root: Path) -> tuple[tuple[str, str, str], ...]:
    records: list[tuple[str, str, str]] = []
    for path in sorted(root.rglob("*"), key=lambda item: str(item)):
        relative = str(path.relative_to(root))
        mode = path.lstat().st_mode
        if stat.S_ISLNK(mode):
            records.append((relative, "symlink", os.readlink(path)))
        elif stat.S_ISREG(mode):
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            records.append((relative, "file", digest))
        elif stat.S_ISDIR(mode):
            records.append((relative, "directory", ""))
        else:
            records.append((relative, "special", str(mode)))
    return tuple(records)


def require_failure_unchanged(
    cwd: Path,
    args: list[str],
    label: str,
    expected_role: str | None = None,
    expected_path: str | None = None,
) -> None:
    before = tree_snapshot(cwd)
    result = run_sa(cwd, args)
    if result.returncode == 0:
        raise AssertionError(f"{label} unexpectedly succeeded")
    if expected_role and expected_role not in result.stderr:
        raise AssertionError(
            f"{label} did not name input role {expected_role!r}:\n{result.stderr}")
    if expected_path and expected_path not in result.stderr:
        raise AssertionError(
            f"{label} did not preserve provided path {expected_path!r}:\n{result.stderr}")
    if "Cannot open input path" not in result.stderr:
        raise AssertionError(f"{label} did not report strict lookup failure:\n{result.stderr}")
    after = tree_snapshot(cwd)
    if after != before:
        raise AssertionError(f"{label} changed files before strict preflight completed")
    print(f"[PASS] {label} fails before workdir/output and preserves prior files")


def matrix(cwd: Path, species: str, gene: str, label: str) -> None:
    require_success(
        run_sa(cwd, [
            "--mode", "matrix", "--species", species,
            "--gene", gene, "--label", label,
        ]),
        f"matrix seed {label}",
    )


def copy_seed(seed: Path, destination: Path) -> None:
    shutil.copytree(seed, destination)


def invalid_forms(root: Path, filename: str, content: bytes) -> dict[str, str]:
    real = root / "real"
    real.mkdir(parents=True)
    (real / "notdir").write_bytes(b"regular intermediate\n")
    (real / filename).write_bytes(content)
    (root / "linkfile").symlink_to(Path("real") / "notdir")
    return {
        "symlink_file_parent": str(root / "linkfile" / ".." / filename),
        "regular_file_parent": str(real / "notdir" / ".." / filename),
        "file_dot": f"{real / filename}/.",
        "trailing_slash": f"{real / filename}/",
    }


def replace_with_invalid_symlink(path: Path, fixture_name: str) -> str:
    content = path.read_bytes()
    target = path.parent / fixture_name / "real" / path.name
    target.parent.mkdir(parents=True)
    target.write_bytes(content)
    path.unlink()
    invalid_target = f"{os.path.relpath(target, path.parent)}/"
    path.symlink_to(invalid_target)
    return str(path)


def finalize_args(
    free: str = "free.matrix_with_fuse.txt",
    fix: str = "fix.matrix_with_fuse.txt",
    final: str = "final",
    free_manifest: str | None = None,
    fix_manifest: str | None = None,
    species_tree: str | None = None,
) -> list[str]:
    args = [
        "--mode", "finalize", "--free", free, "--fix", fix,
        "--final_label", final,
    ]
    if free_manifest is not None:
        args.extend(["--free_manifest", free_manifest])
    if fix_manifest is not None:
        args.extend(["--fix_manifest", fix_manifest])
    if species_tree is not None:
        args.extend(["--species_tree", species_tree])
    return args


def finalize_fix_args(
    fix: str = "fix.matrix_with_fuse.txt",
    final: str = "fix_final",
    manifest: str | None = None,
) -> list[str]:
    args = ["--mode", "finalize_fix", "--fix", fix, "--final_label", final]
    if manifest is not None:
        args.extend(["--fix_manifest", manifest])
    return args


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="sa-recert010-controller-", dir="/tmp") as tmp:
        base = Path(tmp)
        species = b"((A:1,B:1):2,(C:1,D:1):3);\n"
        free_gene = b"g((A:11,B:12):13,(C:14,D:15):16);\n"
        fix_gene = b"g((A:21,B:22):23,(C:24,D:25):26);\n"

        seed = base / "seed"
        seed.mkdir()
        (seed / "species.nwk").write_bytes(species)
        (seed / "free.nwk").write_bytes(free_gene)
        (seed / "fix.nwk").write_bytes(fix_gene)
        matrix(seed, "species.nwk", "free.nwk", "free")
        matrix(seed, "species.nwk", "fix.nwk", "fix")

        # STRICT-01..05: species and gene roles reject all invalid forms.
        for role, option, filename, content, valid_name in (
            ("matrix species tree", "--species", "species.nwk", species, "genes.nwk"),
            ("matrix gene trees", "--gene", "genes.nwk", free_gene, "species.nwk"),
        ):
            for form_name in ("symlink_file_parent", "regular_file_parent",
                              "file_dot", "trailing_slash"):
                case = base / f"matrix_{option[2:]}_{form_name}"
                case.mkdir()
                if option == "--species":
                    (case / valid_name).write_bytes(free_gene)
                else:
                    (case / valid_name).write_bytes(species)
                forms = invalid_forms(case / "invalid", filename, content)
                provided = forms[form_name]
                args = [
                    "--mode", "matrix", "--species", "species.nwk",
                    "--gene", "genes.nwk", "--label", "bad",
                ]
                args[args.index(option) + 1] = provided
                require_failure_unchanged(case, args,
                    f"STRICT matrix {option[2:]} {form_name}", role, provided)

        # STRICT-06: FREE and FIX matrix inputs reject every invalid form.
        for role_name, matrix_role in (("free", "FREE matrix"), ("fix", "FIX matrix")):
            source_name = f"{role_name}.matrix_with_fuse.txt"
            for form_name in ("symlink_file_parent", "regular_file_parent",
                              "file_dot", "trailing_slash"):
                case = base / f"finalize_{role_name}_matrix_{form_name}"
                copy_seed(seed, case)
                forms = invalid_forms(case / "invalid", source_name,
                                      (case / source_name).read_bytes())
                provided = forms[form_name]
                args = finalize_args(**{role_name: provided}, final="bad")
                require_failure_unchanged(case, args,
                    f"STRICT finalize {role_name.upper()} matrix {form_name}",
                    matrix_role, provided)

        # Explicit FREE/FIX manifests reject every invalid form.
        for role_name, manifest_role in (("free", "FREE run manifest"),
                                         ("fix", "FIX run manifest")):
            source_name = f"{role_name}.run_manifest.json"
            for form_name in ("symlink_file_parent", "regular_file_parent",
                              "file_dot", "trailing_slash"):
                case = base / f"explicit_{role_name}_manifest_{form_name}"
                copy_seed(seed, case)
                forms = invalid_forms(case / "invalid", source_name,
                                      (case / source_name).read_bytes())
                provided = forms[form_name]
                kwargs = {f"{role_name}_manifest": provided}
                args = finalize_args(final="bad", **kwargs)
                require_failure_unchanged(case, args,
                    f"STRICT explicit {role_name.upper()} manifest {form_name}",
                    manifest_role, provided)

        # Inferred manifests and manifest-declared state sidecars use strict lookup.
        for role_name, expected_role in (("free", "FREE run manifest"),
                                         ("fix", "FIX run manifest")):
            case = base / f"inferred_{role_name}_manifest"
            copy_seed(seed, case)
            provided = replace_with_invalid_symlink(
                case / f"{role_name}.run_manifest.json", f"bad-{role_name}-manifest")
            require_failure_unchanged(case, finalize_args(final="bad"),
                f"STRICT inferred {role_name.upper()} manifest", expected_role, provided)

        for role_name, expected_role in (("free", "FREE coordinate-state sidecar"),
                                         ("fix", "FIX coordinate-state sidecar")):
            case = base / f"state_{role_name}"
            copy_seed(seed, case)
            provided = replace_with_invalid_symlink(
                case / f"{role_name}.primitive_state.tsv", f"bad-{role_name}-state")
            require_failure_unchanged(case, finalize_args(final="bad"),
                f"STRICT {role_name.upper()} state sidecar", expected_role, provided)

        # STRICT-07: Support tree and branch-map paths use strict lookup.
        for form_name in ("symlink_file_parent", "regular_file_parent",
                          "file_dot", "trailing_slash"):
            case = base / f"support_tree_{form_name}"
            copy_seed(seed, case)
            forms = invalid_forms(case / "invalid", "species_tree.forSplit.nwk",
                                  (case / "species_tree.forSplit.nwk").read_bytes())
            provided = forms[form_name]
            require_failure_unchanged(case,
                finalize_args(final="bad", species_tree=provided),
                f"STRICT Support tree {form_name}", "Support species tree", provided)

        support_map = base / "support_branch_map"
        copy_seed(seed, support_map)
        provided_map = replace_with_invalid_symlink(
            support_map / "species_tree.branch_map.txt", "bad-support-map")
        require_failure_unchanged(support_map,
            finalize_args(final="bad", species_tree="species_tree.forSplit.nwk"),
            "STRICT Support branch map", "Support branch map", provided_map)

        # finalize_fix must apply the same rules to matrix, manifest, and state.
        for form_name in ("symlink_file_parent", "regular_file_parent",
                          "file_dot", "trailing_slash"):
            case = base / f"fix_only_matrix_{form_name}"
            copy_seed(seed, case)
            forms = invalid_forms(case / "invalid", "fix.matrix_with_fuse.txt",
                                  (case / "fix.matrix_with_fuse.txt").read_bytes())
            provided = forms[form_name]
            require_failure_unchanged(case, finalize_fix_args(fix=provided),
                f"STRICT finalize_fix matrix {form_name}", "FIX matrix", provided)

        for form_name in ("symlink_file_parent", "regular_file_parent",
                          "file_dot", "trailing_slash"):
            case = base / f"fix_only_manifest_{form_name}"
            copy_seed(seed, case)
            forms = invalid_forms(case / "invalid", "fix.run_manifest.json",
                                  (case / "fix.run_manifest.json").read_bytes())
            provided = forms[form_name]
            require_failure_unchanged(case,
                finalize_fix_args(manifest=provided),
                f"STRICT finalize_fix explicit manifest {form_name}",
                "FIX run manifest", provided)

        fix_inferred = base / "fix_only_inferred_manifest"
        copy_seed(seed, fix_inferred)
        provided = replace_with_invalid_symlink(
            fix_inferred / "fix.run_manifest.json", "bad-fix-only-manifest")
        require_failure_unchanged(fix_inferred, finalize_fix_args(),
            "STRICT finalize_fix inferred manifest", "FIX run manifest", provided)

        fix_state = base / "fix_only_state"
        copy_seed(seed, fix_state)
        provided = replace_with_invalid_symlink(
            fix_state / "fix.primitive_state.tsv", "bad-fix-only-state")
        require_failure_unchanged(fix_state, finalize_fix_args(),
            "STRICT finalize_fix state", "FIX coordinate-state sidecar", provided)

        # STRICT-08: a directory symlink followed by '..' remains valid.
        valid_link = base / "valid_directory_symlink"
        (valid_link / "real" / "sub").mkdir(parents=True)
        (valid_link / "real" / "species.nwk").write_bytes(species)
        (valid_link / "genes.nwk").write_bytes(free_gene)
        (valid_link / "link").symlink_to(Path("real") / "sub", target_is_directory=True)
        valid_species = str(valid_link / "link" / ".." / "species.nwk")
        matrix(valid_link, valid_species, "genes.nwk", "valid")
        if (valid_link / "species_tree.splits.txt").read_bytes() != \
                (seed / "species_tree.splits.txt").read_bytes():
            raise AssertionError("valid directory symlink selected the wrong species tree")
        print("[PASS] STRICT-08 valid directory symlink/.. selects intended input")

        # STRICT-09: Unicode and multiple directory-symlink hops remain valid.
        unicode = base / "路径_猫"
        (unicode / "真实" / "内层").mkdir(parents=True)
        (unicode / "真实" / "species.nwk").write_bytes(species)
        (unicode / "genes.nwk").write_bytes(free_gene)
        (unicode / "第二跳").symlink_to(Path("真实") / "内层", target_is_directory=True)
        (unicode / "第一跳").symlink_to("第二跳", target_is_directory=True)
        unicode_species = str(
            unicode / "第一跳" / ".." / ".." / "真实" / "species.nwk")
        matrix(unicode, unicode_species, "genes.nwk", "unicode")
        if (unicode / "species_tree.splits.txt").read_bytes() != \
                (seed / "species_tree.splits.txt").read_bytes():
            raise AssertionError("Unicode multiple-hop path selected the wrong input")
        print("[PASS] STRICT-09 Unicode multiple-hop path selects intended input")

        # Nonregular, dangling, looping, and permission-denied inputs fail early.
        special = base / "special_inputs"
        special.mkdir()
        (special / "genes.nwk").write_bytes(free_gene)
        (special / "directory").mkdir()
        os.mkfifo(special / "fifo", 0o600)
        (special / "dangling").symlink_to("missing")
        (special / "loop-a").symlink_to("loop-b")
        (special / "loop-b").symlink_to("loop-a")
        for name, provided in (
            ("directory", str(special / "directory")),
            ("FIFO", str(special / "fifo")),
            ("device", "/dev/null"),
            ("dangling", str(special / "dangling")),
            ("loop", str(special / "loop-a")),
        ):
            require_failure_unchanged(special, [
                "--mode", "matrix", "--species", provided,
                "--gene", "genes.nwk", "--label", "bad",
            ], f"STRICT nonregular {name}", "matrix species tree", provided)

        denied = special / "denied"
        denied.mkdir()
        (denied / "species.nwk").write_bytes(species)
        denied.chmod(0)
        try:
            require_failure_unchanged(special, [
                "--mode", "matrix", "--species", str(denied / "species.nwk"),
                "--gene", "genes.nwk", "--label", "bad",
            ], "STRICT permission-denied input", "matrix species tree",
                str(denied / "species.nwk"))
        finally:
            denied.chmod(0o700)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
