#!/usr/bin/env python3
"""Reject invalid binary64 lengths before publishing a successful matrix."""
import csv
import os
from pathlib import Path
import subprocess
import tempfile

DRIVER = Path(__file__).resolve().parents[2] / "SplitAligner.pl"
REF = "((A,B),(C,D));"


def run(root, tree):
    (root / "reference").write_text(REF + "\n", encoding="utf-8")
    (root / "genes").write_text("gene" + tree + "\n", encoding="utf-8")
    return subprocess.run(
        ["perl", str(DRIVER), "--mode", "matrix", "--species", "reference",
         "--gene", "genes", "--label", "result"],
        cwd=root, capture_output=True, text=True,
        env={**os.environ, "LC_ALL": "C"},
    )


with tempfile.TemporaryDirectory(prefix="perl-numeric-") as tmp:
    root = Path(tmp)
    for i, (tree, diagnostic) in enumerate([
        ("((A:1e999,B:1):2,(C:1,D:1):3);", "not finite"),
        ("((A:1e-999,B:1):2,(C:1,D:1):3);", "underflows to zero"),
        ("((A:1e-400,B:1):2,(C:1,D:1):3);", "underflows to zero"),
        ("((A:2e-324,B:1):2,(C:1,D:1):3);", "underflows to zero"),
        ("((A:\u0661,B:1):2,(C:1,D:1):3);", "Invalid branch value"),
        ("((A:1,B:1):1e308,(C:1,D:1):1e308);", "Pooled unrooted branch"),
    ]):
        case = root / f"invalid-{i}"
        case.mkdir()
        result = run(case, tree)
        assert result.returncode != 0, result.stdout + result.stderr
        assert diagnostic in result.stderr, result.stderr
        assert not (case / "result.run_manifest.json").exists()
        assert not (case / "result.matrix_no_fuse.txt").exists()

    for name, a, b in [
        ("zero", "-0", "+0e999"),
        ("subnormal", "5e-324", "1"),
        ("negative", "-0.1", ".2"),
        ("large", "1e308", "1e308"),
    ]:
        case = root / name
        case.mkdir()
        result = run(case, f"((A:{a},B:{b}):2,(C:1,D:1):3);")
        assert result.returncode == 0, result.stderr
        with (case / "result.matrix_no_fuse.txt").open() as handle:
            row = next(csv.DictReader(handle, delimiter="\t"))
        assert row["B1"] == a and row["B2"] == b
        assert row["B5"] == "5"
    print("PASS: finite input range, pooled overflow, missing publication, zero, subnormal and negative values")
