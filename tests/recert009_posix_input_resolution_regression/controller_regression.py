#!/usr/bin/env python3
"""Controller-level regressions for RECERT-009 POSIX input resolution."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile


ROOT = Path(sys.argv[1]).resolve()
SA = ROOT / "SplitAligner.pl"


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_sa(cwd: Path, args: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["perl", str(SA), *args], cwd=cwd, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
    )


def require_success(result: subprocess.CompletedProcess[str], label: str) -> None:
    if result.returncode != 0:
        raise AssertionError(f"{label} failed:\n{result.stderr}")
    print(f"[PASS] {label}")


def require_failure(result: subprocess.CompletedProcess[str], label: str) -> None:
    if result.returncode == 0:
        raise AssertionError(f"{label} unexpectedly succeeded")
    print(f"[PASS] {label} fails closed")


def write_matrix_inputs(path: Path, species: str, gene: str) -> None:
    path.mkdir(parents=True, exist_ok=True)
    (path / "species.nwk").write_text(species + "\n", encoding="utf-8")
    (path / "genes.nwk").write_text(gene + "\n", encoding="utf-8")


def matrix(path: Path, species: str, gene: str, label: str) -> None:
    require_success(run_sa(path, [
        "--mode", "matrix", "--species", species,
        "--gene", gene, "--label", label,
    ]), f"matrix {label}")


def make_parent_symlink(root: Path, filename: str, content: str,
                        decoy: str | None = None) -> str:
    (root / "real" / "sub").mkdir(parents=True)
    (root / "real" / filename).write_text(content, encoding="utf-8")
    if decoy is not None:
        (root / filename).write_text(decoy, encoding="utf-8")
    (root / "link").symlink_to(Path("real") / "sub", target_is_directory=True)
    return str(root / "link" / ".." / filename)


def manifest_input(path: Path, role: str) -> dict[str, str]:
    data = json.loads(path.read_text(encoding="utf-8"))
    return data["inputs"][role]


def copy_seed(seed: Path, destination: Path) -> None:
    shutil.copytree(seed, destination)


def main() -> int:
    with tempfile.TemporaryDirectory(prefix="splitaligner-recert009-controller-") as tmp:
        base = Path(tmp)
        species_intended = "((A:1,B:1):2,(C:1,D:1):3);"
        species_decoy = "((A:1,C:1):2,(B:1,D:1):3);"
        gene_intended = "g((A:11,B:12):13,(C:14,D:15):16);"
        gene_decoy = "g((A:101,B:102):103,(C:104,D:105):106);"

        reference = base / "reference"
        write_matrix_inputs(reference, species_intended, gene_intended)
        matrix(reference, "species.nwk", "genes.nwk", "ref")

        # PATH-009-01/03: species paths with a symlink before '..'.
        for name, decoy in (("species_decoy_present", species_decoy + "\n"),
                            ("species_decoy_absent", None)):
            case = base / name
            case.mkdir()
            provided = make_parent_symlink(
                case / "species_path", "species.nwk", species_intended + "\n", decoy)
            (case / "genes.nwk").write_text(gene_intended + "\n", encoding="utf-8")
            matrix(case, provided, "genes.nwk", "x")
            assert (case / "species_tree.splits.txt").read_bytes() == \
                (reference / "species_tree.splits.txt").read_bytes()
            record = manifest_input(case / "x.run_manifest.json", "species_tree")
            assert record["sha256"] == sha((case / "species_path" / "real" / "species.nwk"))
            print(f"[PASS] PATH-009 species controller {name}")

        # PATH-009-02/04: gene paths select intended values, never lexical decoys.
        for name, decoy in (("gene_decoy_present", gene_decoy + "\n"),
                            ("gene_decoy_absent", None)):
            case = base / name
            case.mkdir()
            (case / "species.nwk").write_text(species_intended + "\n", encoding="utf-8")
            provided = make_parent_symlink(
                case / "gene_path", "genes.nwk", gene_intended + "\n", decoy)
            matrix(case, "species.nwk", provided, "x")
            assert (case / "x.matrix_no_fuse.txt").read_bytes() == \
                (reference / "ref.matrix_no_fuse.txt").read_bytes()
            record = manifest_input(case / "x.run_manifest.json", "gene_trees")
            assert record["sha256"] == sha(case / "gene_path" / "real" / "genes.nwk")
            print(f"[PASS] PATH-009 gene controller {name}")

        # Create authoritative FREE/FIX artifacts for finalization path roles.
        seed = base / "seed"
        write_matrix_inputs(seed, species_intended, gene_intended)
        (seed / "free.nwk").write_text(gene_intended + "\n", encoding="utf-8")
        (seed / "fix.nwk").write_text(
            "g((A:21,B:22):23,(C:24,D:25):26);\n", encoding="utf-8")
        matrix(seed, "species.nwk", "free.nwk", "free")
        matrix(seed, "species.nwk", "fix.nwk", "fix")

        # PATH-009-05/06/07: FREE/FIX matrix and explicit manifests.
        explicit = base / "finalize_explicit"
        copy_seed(seed, explicit / "real")
        (explicit / "real" / "sub").mkdir()
        (explicit / "link").symlink_to(Path("real") / "sub", target_is_directory=True)
        provided_root = explicit / "link" / ".."
        result = run_sa(explicit, [
            "--mode", "finalize",
            "--free", str(provided_root / "free.matrix_with_fuse.txt"),
            "--fix", str(provided_root / "fix.matrix_with_fuse.txt"),
            "--free_manifest", str(provided_root / "free.run_manifest.json"),
            "--fix_manifest", str(provided_root / "fix.run_manifest.json"),
            "--final_label", "explicit",
        ])
        require_success(result, "PATH-009-05/06/07 matrices and explicit manifests")
        assert (explicit / "explicit.free.na_classified.txt").is_file()
        assert (explicit / "explicit.fix.na_classified.txt").is_file()

        # PATH-009-08: manifests are inferred from each resolved matrix directory.
        inferred = base / "finalize_inferred"
        copy_seed(seed, inferred / "real")
        (inferred / "real" / "sub").mkdir()
        (inferred / "link").symlink_to(Path("real") / "sub", target_is_directory=True)
        inferred_root = inferred / "link" / ".."
        require_success(run_sa(inferred, [
            "--mode", "finalize",
            "--free", str(inferred_root / "free.matrix_with_fuse.txt"),
            "--fix", str(inferred_root / "fix.matrix_with_fuse.txt"),
            "--final_label", "inferred",
        ]), "PATH-009-08 inferred manifests use resolved matrix directory")

        # PATH-009-09: a manifest-declared state filename may itself be a symlink.
        state_case = base / "state_sidecar"
        copy_seed(seed, state_case / "real")
        (state_case / "real" / "sub").mkdir()
        (state_case / "link").symlink_to(Path("real") / "sub", target_is_directory=True)
        state_target = state_case / "state_targets"
        state_target.mkdir()
        for label in ("free", "fix"):
            source = state_case / "real" / f"{label}.primitive_state.tsv"
            target = state_target / f"{label}.primitive_state.tsv"
            source.rename(target)
            source.symlink_to(target)
        state_root = state_case / "link" / ".."
        require_success(run_sa(state_case, [
            "--mode", "finalize",
            "--free", str(state_root / "free.matrix_with_fuse.txt"),
            "--fix", str(state_root / "fix.matrix_with_fuse.txt"),
            "--final_label", "state",
        ]), "PATH-009-09 state sidecars resolve to hashed targets")

        # PATH-009-10: Support tree uses POSIX resolution; branch map is derived
        # from the resolved tree directory and may itself be a symlink.
        support = base / "support"
        copy_seed(seed, support / "real")
        (support / "real" / "sub").mkdir()
        (support / "link").symlink_to(Path("real") / "sub", target_is_directory=True)
        support_map_target = support / "map-target.tsv"
        map_path = support / "real" / "species_tree.branch_map.txt"
        map_path.rename(support_map_target)
        map_path.symlink_to(support_map_target)
        support_root = support / "link" / ".."
        require_success(run_sa(support, [
            "--mode", "finalize",
            "--free", str(support_root / "free.matrix_with_fuse.txt"),
            "--fix", str(support_root / "fix.matrix_with_fuse.txt"),
            "--species_tree", str(support_root / "species_tree.forSplit.nwk"),
            "--final_label", "support",
        ]), "PATH-009-10 Support tree and derived branch map")
        assert (support / "support.support_b.txt").is_file()

        # PATH-009-06 finalize_fix path and explicit manifest.
        fix_only = base / "fix_only"
        copy_seed(seed, fix_only / "real")
        (fix_only / "real" / "sub").mkdir()
        (fix_only / "link").symlink_to(Path("real") / "sub", target_is_directory=True)
        fix_root = fix_only / "link" / ".."
        require_success(run_sa(fix_only, [
            "--mode", "finalize_fix",
            "--fix", str(fix_root / "fix.matrix_with_fuse.txt"),
            "--fix_manifest", str(fix_root / "fix.run_manifest.json"),
            "--final_label", "fix_only",
        ]), "PATH-009-06 finalize_fix input chain")

        # PATH-009-11/12: Unicode and multiple symlink hops reach the same tree.
        unicode_case = base / "路径_猫"
        (unicode_case / "真实" / "内层").mkdir(parents=True)
        shutil.copy2(seed / "species.nwk", unicode_case / "真实" / "species.nwk")
        shutil.copy2(seed / "free.nwk", unicode_case / "genes.nwk")
        (unicode_case / "第二跳").symlink_to(Path("真实") / "内层", target_is_directory=True)
        (unicode_case / "第一跳").symlink_to("第二跳", target_is_directory=True)
        unicode_provided = unicode_case / "第一跳" / ".." / ".." / "真实" / "species.nwk"
        matrix(unicode_case, str(unicode_provided), "genes.nwk", "unicode")
        assert (unicode_case / "species_tree.splits.txt").read_bytes() == \
            (seed / "species_tree.splits.txt").read_bytes()
        print("[PASS] PATH-009-11/12 Unicode and multiple-hop controller path")

        # PATH-009-14: controller rejects dangling and looping inputs before output.
        failures = base / "failures"
        failures.mkdir()
        (failures / "genes.nwk").write_text(gene_intended + "\n", encoding="utf-8")
        (failures / "dangling").symlink_to("missing")
        require_failure(run_sa(failures, [
            "--mode", "matrix", "--species", "dangling",
            "--gene", "genes.nwk", "--label", "dangling",
        ]), "PATH-009-14 dangling controller input")
        (failures / "loop-a").symlink_to("loop-b")
        (failures / "loop-b").symlink_to("loop-a")
        require_failure(run_sa(failures, [
            "--mode", "matrix", "--species", "loop-a",
            "--gene", "genes.nwk", "--label", "loop",
        ]), "PATH-009-14 looping controller input")
        assert not (failures / "dangling.matrix_with_fuse.txt").exists()
        assert not (failures / "loop.matrix_with_fuse.txt").exists()

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
