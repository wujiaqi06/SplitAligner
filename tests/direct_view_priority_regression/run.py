#!/usr/bin/env python3
"""Check direct-view priority in Perl and optionally a native Linux Fast CLI.

From the repository root:
  python3 tests/direct_view_priority_regression/run.py --fast-cli fast/build/saf
An optional --r-evidence CSV from the corresponding R fixture is compared too.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import tempfile

REPO = Path(__file__).resolve().parents[2]
REFERENCE = "((A,(B,C)),(D,E));"
GENES = {
    "del": "((A:1,B:1),(D:1,E:1));",
    "overflow": "((A:1,(B:1e308,C:3):1.5e308):5,(D:6,E:7):8);",
    "missing": "((A:1,(B,C:3):1.5e308):5,(D:6,E:7):8);",
    "zero": "((A:1,(B:0,C:3):0):5,(D:6,E:7):8);",
    "finite": "((A:1,(B:2,C:3):3):5,(D:6,E:7):8);",
    "absent": "(A:1,B:2,C:3,D:6,E:7);",
    "deleted": "(A:1,C:3,(D:6,E:7):8);",
    "active_missing": "((A:1,B),(D:1,E:1));",
}
FIELDS = ("reference_state", "recovery_status", "numeric_status", "construction_mode")
EXPECTED = {
    "del": ("active", "recovered", "available", "single_edge_lookup", 1),
    "overflow": ("decomposable", "recovered", "nonfinite_sum", "selected_component_edge_sum", None),
    "missing": ("decomposable", "recovered", "missing_length", "selected_component_edge_sum", None),
    "zero": ("decomposable", "recovered", "available", "selected_component_edge_sum", 0),
    "finite": ("decomposable", "recovered", "available", "selected_component_edge_sum", 5),
    "absent": ("decomposable", "not_recovered", "not_recovered", "selected_component_edge_sum", None),
    "deleted": ("unavailable", "not_applicable", "not_applicable", "none", None),
    "active_missing": ("active", "recovered", "missing_length", "single_edge_lookup", None),
}


def tsv(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def jsonl(path):
    return [json.loads(line) for line in path.read_text().splitlines()]


def call(command, root, success=True):
    number = len(list(root.glob("command-*.json"))) + 1
    log = root / f"command-{number:04d}"
    log.with_suffix(".json").write_text(json.dumps(command, indent=2) + "\n")
    run = subprocess.run(command, cwd=root, capture_output=True, text=True,
                         env={**os.environ, "LC_ALL": "C"}, timeout=30)
    log.with_suffix(".stdout.log").write_text(run.stdout)
    log.with_suffix(".stderr.log").write_text(run.stderr)
    assert (run.returncode == 0) == success, run.stdout + run.stderr
    return run


def check(cells, view, value_field):
    actual = {}
    for gene, row in cells.items():
        expected = EXPECTED[gene]
        if view == "direct" and expected[:2] == ("decomposable", "recovered"):
            expected = (*expected[:2], "direct_only", expected[3], None)
        value = row[value_field]
        value = None if value is None or value in ("NA", "") else float(value)
        actual[gene] = (*[row[field] for field in FIELDS], value)
        assert actual[gene] == expected, (view, gene, actual[gene], expected)
    assert set(actual) == set(GENES)
    return actual


def verify(root, fast=None, r_evidence=None):
    (root / "reference").write_text(REFERENCE + "\n")
    (root / "genes").write_text("".join(gene + tree + "\n" for gene, tree in GENES.items()))
    (root / "genes.tsv").write_text("".join(gene + "\t" + tree + "\n" for gene, tree in GENES.items()))
    all_results = {}
    for view in ("whole", "direct"):
        call(["perl", str(REPO / "SplitAligner.pl"), "--mode", "matrix", "--species", "reference",
              "--gene", "genes", "--label", view, "--fusion-view", view], root)
        cells = {row["gene"]: row for row in tsv(root / f"{view}.fusion_cells.tsv")
                 if row["column"] == "B2|B6"}
        all_results[view] = check(cells, view, "length")
        if fast:
            out = root / f"fast-{view}"
            call([str(fast), "--reference", str(root / "reference"), "--genes", str(root / "genes.tsv"),
                  "--gene-format", "tsv", "--scope", "full", "--role", "free", "--unmatched-splits", "on",
                  "--out", str(out), "--fusion-view", view], root)
            target, = [row["members"] for row in jsonl(out / "composites.jsonl") if row["gene_id"] == "del"]
            fast_cells = {row["gene_id"]: row for row in jsonl(out / "fusion_cells.jsonl")
                          if row["members"] == target}
            assert check(fast_cells, view, "value") == all_results[view]
            evidence = {gene: row["components"] for gene, row in fast_cells.items()}
            if view == "whole": fast_evidence = evidence
            else: assert evidence == fast_evidence
    for suffix in ("primitive_state.tsv", "matrix_no_fuse.txt"):
        assert (root / f"whole.{suffix}").read_bytes() == (root / f"direct.{suffix}").read_bytes()
    for view in ("whole", "direct"):
        meta = json.loads((root / f"{view}.run_manifest.json").read_text())
        counts = meta["projection_graph"]["primitive_counts"]
        if view == "whole": expected_counts = counts
        else: assert counts == expected_counts
    # Compact JSONL replaces per-gene component TSVs. Comparing every byte
    # retains component/observation/source-edge and input-identity invariance.
    ledgers = [root / f"{view}.compact.jsonl" for view in ("whole", "direct")]
    assert ledgers[0].read_bytes() == ledgers[1].read_bytes()
    records = jsonl(ledgers[0])
    assert len(records) == len(GENES)
    assert {r["gene_id"] for r in records} == set(GENES)
    for record in records:
        observations = record["observations"]
        for component in record["components"]:
            oid = component["observation_id"]
            if oid is not None:
                assert type(oid) is int and 1 <= oid <= len(observations)
                assert observations[oid - 1]["key"] == component["key"]
    if fast:
        for name in ("matrix.tsv", "raw_states.tsv", "support.tsv", "composites.jsonl", "unmatched.jsonl"):
            assert (root / "fast-whole" / name).read_bytes() == (root / "fast-direct" / name).read_bytes()
    bad = "(A:1e308,(B:1,C:1,D:1,E:1):1e308);"
    (root / "bad").write_text("bad" + bad + "\n")
    run = call(["perl", str(REPO / "SplitAligner.pl"), "--mode", "matrix", "--species", "reference",
                "--gene", "bad", "--label", "bad", "--fusion-view", "direct"], root, success=False)
    assert "Pooled unrooted branch" in run.stderr
    assert not (root / "bad.run_manifest.json").exists()
    if fast:
        (root / "bad.tsv").write_text("bad\t" + bad + "\n")
        run = call([str(fast), "--reference", str(root / "reference"), "--genes", str(root / "bad.tsv"),
                    "--gene-format", "tsv", "--scope", "full", "--role", "free", "--unmatched-splits", "on",
                    "--out", str(root / "fast-bad"), "--fusion-view", "direct"], root, success=False)
        assert "finite source aggregation overflow" in run.stderr
        assert not (root / "fast-bad").exists()
    if r_evidence:
        with r_evidence.open(encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle))
        for view in ("whole", "direct"):
            cells = {row["gene_id"]: row for row in rows if row["view"] == view}
            assert check(cells, view, "numeric_value") == all_results[view]
    print(json.dumps(all_results, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fast-cli", type=Path)
    parser.add_argument("--r-evidence", type=Path)
    parser.add_argument("--work-dir", type=Path, help="new directory to retain outputs and command logs")
    args = parser.parse_args()
    fast = args.fast_cli.resolve() if args.fast_cli else None
    evidence = args.r_evidence.resolve() if args.r_evidence else None
    if args.work_dir:
        root = args.work_dir.resolve()
        root.mkdir()
        verify(root, fast, evidence)
    else:
        with tempfile.TemporaryDirectory(prefix="direct-view-priority-") as temp:
            verify(Path(temp), fast, evidence)
    print("PASS: direct-only priority, complete component evidence, unchanged topology, fatal single-observation overflow")
