#!/usr/bin/env python3
"""Trace a three-branch target through two ways of obtaining the same length."""
import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import tempfile

DRIVER = Path(__file__).resolve().parents[2] / "SplitAligner.pl"
REF = "(((((A,B),C),D),E),F,(G,H));"
GENES = {
    "independent": "(((((A,B):3,C):5,D):7,E),F,(G,H));",
    "pair": "((((A,B):8,D):7,E),F,(G,H));",
    "whole": "(((A,B):15,E),F,(G,H));",
    "pooled": "((A,B):4,(E,(F,(G,H))):11);",
    "outside": "((A,B):24,F,(G,H));",
    "deleted": "(((C,D):7,E),F,(G,H));",
    "partial": "(((A,B):8,(D,E):7),F,(G,H));",
    "lengthless": "((((A,B),D):7,E),F,(G,H));",
    "zero": "((((A,B):0,D):0,E),F,(G,H));",
    "absent": "((A,E),B,F,(G,H));",
    "precision": "(((((A,B):0.1,C):0.2,D):0,E),F,(G,H));",
    "integer_rounding": "(((((A,B):9007199254740993,C):1,D):0,E),F,(G,H));",
    # No other gene in this panel restricts the reference by E alone.
    "unique_miss": "((A,F),B,C,D,G,H);",
}


def run(cwd, *args):
    subprocess.run(["perl", str(DRIVER), *args], cwd=cwd, check=True,
                   env={**os.environ, "LC_ALL": "C"}, capture_output=True)


def rows(path):
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def verify(root):
    shared_graphs = []
    for label in ("whole", "direct"):
        (root / "ref.nwk").write_text(REF + "\n")
        (root / "genes.nwk").write_text("".join(gene + tree + "\n" for gene, tree in GENES.items()))
        args = ["--mode", "matrix", "--species", "ref.nwk", "--gene", "genes.nwk", "--label", label]
        run(root, *(args if label == "whole" else args + ["--fusion-view", label]))
        graph = root / (label + "_split_branch_label")
        groups = rows(graph / (label + ".projected_groups.tsv"))
        target = next(r["members"] for r in groups if r["gene"] == "whole")
        nested = next(r["members"] for r in groups if r["gene"] == "pair")
        miss = next(r["members"] for r in groups if r["gene"] == "unique_miss")
        assert [r["status"] for r in groups if r["members"] == miss] == ["not_displayed"]
        cells = {(r["gene"], r["column"]): r for r in rows(root / (label + ".fusion_cells.tsv"))}
        matrix = {r["gene"]: r for r in rows(root / (label + ".matrix_with_fuse.txt"))}
        assert matrix["unique_miss"][miss] == "NA"
        assert cells["unique_miss", miss]["reference_state"] == "active"
        assert len(target.split("|")) == 3 and len(nested.split("|")) == 2
        for gene, state, recovery, numeric, mode, value in [
            ("independent", "decomposable", "recovered", "available", "selected_component_edge_sum", "15"),
            ("pair", "decomposable", "recovered", "available", "selected_component_edge_sum", "15"),
            ("whole", "active", "recovered", "available", "single_edge_lookup", "15"),
            ("pooled", "active", "recovered", "available", "single_edge_lookup", "15"),
            ("outside", "unavailable", "not_applicable", "not_applicable", "none", "NA"),
            ("deleted", "unavailable", "not_applicable", "not_applicable", "none", "NA"),
            ("partial", "decomposable", "not_recovered", "not_recovered", "selected_component_edge_sum", "NA"),
            ("lengthless", "decomposable", "recovered", "missing_length", "selected_component_edge_sum", "NA"),
            ("zero", "decomposable", "recovered", "available", "selected_component_edge_sum", "0"),
            ("absent", "active", "not_recovered", "not_recovered", "single_edge_lookup", "NA"),
            ("integer_rounding", "decomposable", "recovered", "available", "selected_component_edge_sum", "9007199254740992"),
            ("precision", "decomposable", "recovered", "available", "selected_component_edge_sum", str(0.1 + 0.2)),
        ]:
            c = cells[gene, target]
            if recovery == "recovered" and mode == "selected_component_edge_sum" and label == "direct":
                numeric, value = "direct_only", "NA"
            assert (c["reference_state"], c["recovery_status"], c["numeric_status"],
                    c["construction_mode"], c["length"]) == (state, recovery, numeric, mode, value), c
            assert matrix[gene][target] == c["length"]
        assert nested in cells["pair", target]["components"].split(",")
        assert cells["outside", target]["cut_components"] in cells["outside", target]["components"].split(",")
        assert cells["deleted", target]["deleted_members"]
        # Follow the complete compact component to its normalized observation.
        # Decode the documented LSB-first bits independently of the Perl core.
        records = [json.loads(line) for line in
                   (root / (label + ".compact.jsonl")).read_text().splitlines()]
        pooled, = [r for r in records if r["gene_id"] == "pooled"]
        component, = [r for r in pooled["components"] if r["members"] == target]
        oid = component["observation_id"]
        assert type(oid) is int and 1 <= oid <= len(pooled["observations"])
        observation = pooled["observations"][oid - 1]
        assert observation["key"] == component["key"]
        dictionary = json.loads((root / (label + ".taxon_dictionary.json")).read_text())
        assert dictionary["taxa"] == list("ABCDEFGH")
        assert dictionary["bit_order"] == "vec-LSB-first"
        def members(bits):
            raw = bytes.fromhex(bits)
            return {taxon for i, taxon in enumerate(dictionary["taxa"])
                    if raw[i // 8] & (1 << (i % 8))}
        retained = members(pooled["tg"])
        side = members(observation["key"])
        assert retained == set("ABEFGH")
        assert {frozenset(side), frozenset(retained - side)} == {
            frozenset("AB"), frozenset("EFGH")}
        assert observation["length"] == "15"
        assert observation["source_nodes"] == [1, 4]
        manifest = json.loads((root / (label + ".run_manifest.json")).read_text())
        assert manifest["fusion_matrix"]["view"] == label
        shared_graphs.append([groups, manifest["projection_graph"]["primitive_counts"],
                              (root / (label + ".primitive_state.tsv")).read_bytes(),
                              (root / (label + ".matrix_no_fuse.txt")).read_bytes()])
    assert shared_graphs[0] == shared_graphs[1]
    saved = (root / "direct.run_manifest.json").read_bytes()
    try:
        run(root, "--mode", "matrix", "--species", "ref.nwk", "--gene", "genes.nwk",
            "--label", "direct", "--fusion-view", "whole", "--force")
    except subprocess.CalledProcessError as error:
        assert b"changing the fusion view" in error.stderr
    else:
        raise AssertionError("a different fusion view overwrote the direct run")
    assert (root / "direct.run_manifest.json").read_bytes() == saved


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", type=Path, help="new directory retaining all fixture outputs")
    args = parser.parse_args()
    if args.work_dir:
        args.work_dir.mkdir()
        verify(args.work_dir.resolve())
    else:
        with tempfile.TemporaryDirectory() as directory:
            verify(Path(directory))
    print("PASS: common fusion columns, whole components, direct view and source edges")
