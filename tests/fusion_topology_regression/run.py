#!/usr/bin/env python3
"""Exercise the primitive and group outputs through both CLI stages."""
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile

CONTROLLER = Path(__file__).resolve().parents[2] / "SplitAligner.pl"
FIGURE4 = {
    "g1": "((A,B):0.31,D,(E,F):0.17);",
    "g2": "((A,B):0.26,E,(D,F));",
    "g3": "((A,D),B,(E,F):0.23);",
}
SPLIT = "A..B||D..E..F"


def cli(cwd, *args, ok=True):
    r = subprocess.run(["perl", str(CONTROLLER), *args], cwd=cwd,
                       env={**os.environ, "LC_ALL": "C"}, capture_output=True, text=True)
    assert (r.returncode == 0) == ok, r.stderr
    return r


def table(path):
    with path.open(encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def run(d):
    (d / "ref.nwk").write_text("((A,B),C,(D,(E,F)));\n")
    genes = {**FIGURE4, "zero": "((A,B):0,D,(E,F));",
             "no_length": "((A,B),D,(E,F));", "star": "(A,B,D,E,F);",
             "pooled": "((A,B):1,(D,(E,F)):2);", "mixed_root": "((A,B):1,(D,(E,F)));",
             "unique_topo": "((A,D),(B,C),(E,F));"}
    fixed = {gene: FIGURE4["g1"] for gene in genes}
    fixed["unique_topo"] = "((A,B):1,C,(D,(E,F)));"
    # The program must obtain the free/fixed role from the paired CLI arguments.
    for label, input_genes in [("alpha", {**genes, "free_only": FIGURE4["g3"]}), ("beta", fixed)]:
        path = d / (label + ".nwk")
        path.write_text("".join(g + tree + "\n" for g, tree in input_genes.items()))
        cli(d, "--mode", "matrix", "--species", "ref.nwk", "--gene", path.name, "--label", label)
    base = ("--mode", "finalize", "--free", "alpha.matrix_with_fuse.txt", "--fix", "beta.matrix_with_fuse.txt")
    cli(d, *base, "--final_label", "out")
    raw_groups_path = d / "alpha_split_branch_label/alpha.projected_groups.tsv"
    raw = table(raw_groups_path)
    assert all(r["status"] in {"recovered", "not_displayed"} for r in raw)
    groups = table(d / "out.free.projected_groups.tsv")
    records = {r["gene"]: r for r in groups if r["projected_split"] == SPLIT}
    assert len(records) == len(groups) == len(genes) - 1
    expected = {"g1": "0.31", "g2": "0.26", "g3": None, "zero": "0", "no_length": "NA",
                "star": None, "pooled": "3", "mixed_root": "NA"}
    for g, value in expected.items():
        r = records[g]
        assert r["status"] == ("NA_fuse_topo" if value is None else "recovered"), r
        assert r["length"] == (value if value is not None else "NA"), r
    members = set(records["g3"]["members"].split("|"))
    assert len(members) == 2
    for side in ["fix", "free"]:
        for row in table(d / ("out." + side + ".na_classified.txt")):
            if row["gene"] in expected:
                assert all(row[m] == "NA_fuse" for m in members), row
            assert "NA_fuse_topo" not in row.values()
    truth = table(d / "alpha_split_branch_label/alpha.primitive_truth.tsv")
    counts = {"Mapped": 0, "NA_fuse": 0, "NA_struct": 0, "not_displayed": 0}
    for row in truth:
        for k, v in row.items():
            if k != "gene":
                counts[v] += 1
    meta = json.loads((d / "alpha.run_manifest.json").read_text())
    assert counts == meta["projection_graph"]["primitive_counts"]
    assert all(n > 0 for n in counts.values())
    assert sum(meta["projection_graph"]["fused_group_counts"].values()) == len(raw)

    accounting = json.loads((d / "out.finalize_manifest.json").read_text())["free_projection_accounting"]
    by_member = accounting["by_reference_primitive"]
    by_coordinate = accounting["by_projected_coordinate_plus_structural_absences"]
    assert (by_member["NA_fuse"], by_member["NA_fuse_topo"]) == (12, 4)
    assert (by_coordinate["NA_fuse"], by_coordinate["NA_fuse_topo"]) == (6, 2)
    assert by_member["Total"] == len(genes) * (len(truth[0]) - 1)
    assert by_coordinate["Total"] == by_member["Total"] - len(groups)
    for tally in accounting.values():
        assert tally["Total"] == sum(v for category, v in tally.items() if category != "Total")

    # A v3 run has no projected-group evidence. Its U is left unresolved
    # when the fixed primitive only carried a fusion label.
    for label in ["alpha", "beta"]:
        manifest = json.loads((d / (label + ".run_manifest.json")).read_text())
        manifest["schema_version"] = "SplitAligner-run-manifest-v3"
        del manifest["projection_graph"]
        (d / (label + ".v3.json")).write_text(json.dumps(manifest))
    r = cli(d, *base, "--final_label", "legacy", "--free_manifest", "alpha.v3.json", "--fix_manifest", "beta.v3.json")
    assert "legacy v3" in r.stderr
    g3 = next(row for row in table(d / "legacy.free.na_classified.txt") if row["gene"] == "g3")
    assert all(g3[m] == "NA" for m in members)
    assert not (d / "legacy.free.projected_groups.tsv").exists()

    # Damaged evidence cannot change an existing run on a forced retry.
    before = (d / "out.finalize_manifest.json").read_bytes()
    original = raw_groups_path.read_text()
    raw_groups_path.write_text(original.replace("not_displayed", "recovered"))
    cli(d, *base, "--final_label", "out", "--force", ok=False)
    assert (d / "out.finalize_manifest.json").read_bytes() == before

    # Even a recomputed file checksum cannot hide a missing group's members.
    raw_groups_path.write_text(original.replace(next(l for l in original.splitlines(True) if l.startswith("g3\t")), ""))
    meta["projection_graph"]["projected_groups"]["sha256"] = hashlib.sha256(raw_groups_path.read_bytes()).hexdigest()
    (d / "incomplete.json").write_text(json.dumps(meta))
    cli(d, *base, "--final_label", "incomplete", "--free_manifest", "incomplete.json", ok=False)
    assert not (d / "incomplete.finalize_manifest.json").exists()


if __name__ == "__main__":
    with tempfile.TemporaryDirectory() as tmp:
        run(Path(tmp))
    print("PASS: paired projected groups, member truth, numeric availability, legacy and damaged inputs")
