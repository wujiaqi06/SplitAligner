"""Construct one three-axis fusion from complete components and check the source edges."""
import csv
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

CLI = sys.argv[1]
REFERENCE = "(((((A,B),C),D),E),F,(G,H));"
GENES = {
    "separate": "(((((A:1,B:2):3,C:4):5,D:6):7,E:8):9,F:10,(G:11,H:12):13);",
    "rounded": "(((((A:1,B:2):1e16,C:4):1,D:6):1,E:8):9,F:10,(G:11,H:12):13);",
    "nested": "((((A:1,B:2):8,D:6):7,E:8):9,F:10,(G:11,H:12):13);",
    "cut": "((A:1,B:2):24,F:10,(G:11,H:12):13);",
    "deleted": "(((C:4,D:6):7,E:8):9,F:10,(G:11,H:12):13);",
    "missing_split": "(((A:1,B:2):8,(D:6,E:8):7):9,F:10,(G:11,H:12):13);",
    "missing_length": "((((A:1,B:2),D:6):7,E:8):9,F:10,(G:11,H:12):13);",
    "zero": "((((A:1,B:2):0,D:6):0,E:8):9,F:10,(G:11,H:12):13);",
    "overflow": "((((A:1,B:2):1.7e308,D:6):1.7e308,E:8):9,F:10,(G:11,H:12):13);",
    "absent_group": "((A:1,E:8):9,B:2,F:10,(G:11,H:12):13);",
    "direct": "((A:1,B:2):4,(E:8,(F:10,(G:11,H:12):13):9):11);",
}


def records(path):
    return [json.loads(line) for line in path.read_bytes().split(b'\n') if line]


class FusionMatrix(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.ref = self.root / "ref"
        self.ref.write_text(REFERENCE)
        self.rows = "".join(f"{gene}\t{tree}\n" for gene, tree in GENES.items()).encode()

    def run_cli(self, name, *options, rows=None):
        out = self.root / name
        path = self.root / "genes.tsv"
        path.write_bytes(self.rows if rows is None else rows)
        subprocess.run([CLI, "--reference", str(self.ref), "--genes", str(path),
                        "--out", str(out), *options], check=True, stdout=subprocess.DEVNULL, timeout=10)
        return out

    def test_whole_direct_off(self):
        whole = self.run_cli("all")
        direct = self.run_cli("direct", "--fusion-view", "direct")
        off = self.run_cli("off", "--fusion-view", "off")
        for file in ("matrix.tsv", "raw_states.tsv", "support.tsv", "support.nwk", "composites.jsonl",
                     "unmatched.jsonl", "unmatched.index.tsv", "coverage.jsonl", "mapping.jsonl"):
            self.assertEqual((whole/file).read_bytes(), (direct/file).read_bytes())
            self.assertEqual((whole/file).read_bytes(), (off/file).read_bytes())
        direct_groups = records(whole/"composites.jsonl")
        group, = [r["members"] for r in direct_groups if r["gene_id"] == "direct" and len(r["members"]) == 3]
        context = json.loads((whole/"context.json").read_bytes())
        splits = {a["id"]: a["split_hex"] for a in context["axes"]}
        # Reference sides AB, ABC and ABCD are the three original coordinates.
        self.assertEqual({int(splits[m], 16) for m in group}, {3, 7, 15})

        label = "|".join(f"A{x}" for x in group)
        cells = {r["gene_id"]: r for r in records(whole/"fusion_cells.jsonl") if r["column"] == label}
        self.assertEqual(set(cells), set(GENES))
        whole_values = [15, 10000000000000002, 15, None, None, None, None, 0, None, None, 15]
        statuses = ["available", "available", "available", "not_applicable", "not_applicable",
                    "not_recovered", "missing_length", "available", "nonfinite_sum",
                    "not_recovered", "available"]
        self.assertEqual([cells[k]["value"] for k in GENES], whole_values)
        self.assertEqual([cells[k]["numeric_status"] for k in GENES], statuses)
        self.assertTrue(cells["cut"]["partial_component"])
        self.assertTrue(cells["deleted"]["deleted_member"])
        self.assertEqual(cells["absent_group"]["reference_state"], "active")
        self.assertEqual(cells["missing_split"]["reference_state"], "decomposable")
        self.assertEqual(cells["missing_length"]["recovery_status"], "recovered")
        self.assertEqual(cells["absent_group"]["recovery_status"], "not_recovered")
        self.assertEqual(cells["separate"]["construction_mode"], "selected_component_edge_sum")
        self.assertEqual(cells["direct"]["construction_mode"], "single_edge_lookup")
        self.assertEqual([len(c["members"]) for c in cells["nested"]["components"]], [2, 1])
        self.assertEqual(len(cells["direct"]["components"][0]["source_nodes"]), 2)

        mapping = {r["gene_id"]: r for r in records(whole/"mapping.jsonl")}
        for gene, cell in cells.items():
            edges = mapping[gene]["gene_edges"]
            for component in cell["components"]:
                edge = component["gene_edge_id"]
                if edge is None:
                    self.assertEqual(component["source_nodes"], [])
                else:
                    self.assertEqual(component["source_nodes"], edges[edge]["source_nodes"])
                    self.assertEqual(component["value"], edges[edge]["value"])
                    self.assertEqual(component["split_hex"], edges[edge]["split_hex"])
        for out, view in ((whole, "whole"), (direct, "direct"), (off, "off")):
            manifest = json.loads((out/"manifest.json").read_bytes())
            self.assertEqual(manifest["fusion_view"], view)
            self.assertEqual((out/"fusion_matrix.tsv").exists(), view != "off")
            if view == "off":
                continue
            cells_view = {r["gene_id"]: r for r in records(out/"fusion_cells.jsonl") if r["column"] == label}
            with (out/"fusion_matrix.tsv").open(encoding="utf-8") as stream:
                for row in csv.DictReader(stream, delimiter="\t"):
                    cell = cells_view[row["gene_id"]]
                    self.assertEqual(None if row[label] == "NA" else float(row[label]), cell["value"])
                    if view == "direct" and cell["recovery_status"] == "recovered" and cell["reference_state"] == "decomposable":
                        self.assertEqual(cell["numeric_status"], "direct_only")
                    self.assertEqual(cell["components"], cells[row["gene_id"]]["components"])

    def test_direct_does_not_hide_single_observation_overflow(self):
        rows = "bad\t(A:1e308,(B:2,C:3,D:4,E:5,F:6,G:7,H:8):1e308);\n"
        path = self.root / "root-overflow.tsv"
        path.write_text(rows)
        out = self.root / "root-overflow"
        run = subprocess.run([CLI, "--reference", str(self.ref), "--genes", str(path),
                              "--out", str(out), "--fusion-view", "direct"],
                             capture_output=True, text=True, timeout=10)
        self.assertNotEqual(run.returncode, 0)
        self.assertIn("finite source aggregation overflow", run.stderr)
        self.assertFalse(out.exists())

    def test_matrix_scope_and_degenerate_runs(self):
        full = self.run_cli("full")
        matrix = self.run_cli("matrix", "--scope", "matrix", "--role", "fixed")
        for file in ("fusion_matrix.tsv", "fusion_cells.jsonl"):
            self.assertEqual((full/file).read_bytes(), (matrix/file).read_bytes())
        for name, options in (("support", ["--scope", "support"]), ("states", ["--no-numeric"])):
            out = self.run_cli(name, *options)
            self.assertFalse((out/"fusion_cells.jsonl").exists())
        out = self.run_cli("no_groups", rows=self.rows.splitlines()[0] + b"\n")
        self.assertEqual((out/"fusion_matrix.tsv").read_text(), "gene_id\nseparate\n")
        self.assertEqual((out/"fusion_cells.jsonl").read_bytes(), b"")
        empty = self.run_cli("empty", rows=b"")
        self.assertEqual((empty/"fusion_matrix.tsv").read_text(), "gene_id\n")


if __name__ == "__main__":
    unittest.main(argv=[sys.argv[0]])
