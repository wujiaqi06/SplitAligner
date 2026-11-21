"""Cross-check the support tree against the TSV and hand-counted unrooted splits."""
import csv
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest


CLI = sys.argv[1]
TIP = "A's, sample"
QUOTED_TIP = "'A''s, sample'"
# Unique lengths on the two root-adjacent edges must stay on their own sides.
SHAPE = (
    f"((({QUOTED_TIP}:1,B:2){{ab}}:3,(C:4,D:5){{cd}}:6){{left}}:7,"
    "((E:8,F:9){ef}:10,(G:11,H:12){gh}:13){right}:14);"
)
REFERENCE = SHAPE.format(ab="old1", cd="old2", left="old3", ef="old4",
                         gh="old5", right="old6")
# Every tree lacks H, so GH projects onto the G terminal edge and is ineligible.
RECOVERED = f"(({QUOTED_TIP},B),C,D,((E,F),G));"
UNRESOLVED = f"(({QUOTED_TIP},B),C,D,(E,F),G);"


def tokens(newick):
    """Ignore optional quoting and whitespace while keeping Newick punctuation."""
    parts = re.findall(r"'(?:[^']|'')*'|[(),:;]|[^\s'(),:;]+|\s+|.", newick)
    if "'" in parts:
        raise ValueError("unterminated Newick label")
    return [p[1:-1].replace("''", "'") if p.startswith("'") else p
            for p in parts if not p.isspace()]


class SupportTree(unittest.TestCase):
    def setUp(self):
        tmp = TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name)
        self.ref = self.root / "reference.nwk"
        self.genes = self.root / "genes.tsv"
        self.ref.write_text(REFERENCE, encoding="utf-8")
        self.genes.write_text(f"g1\t{RECOVERED}\ng2\t{RECOVERED}\ng3\t{UNRESOLVED}\n",
                              encoding="utf-8")

    def run_cli(self, *options, truth=False):
        output = self.root / ("output-" + "-".join(options))
        if truth:
            taxa = self.root / "taxa.tsv"
            taxa.write_text("g1\tB\tC\tD\n", encoding="utf-8")
            inputs = ["--taxa", str(taxa)]
        else:
            inputs = ["--genes", str(self.genes)]
        subprocess.run([CLI, "--reference", str(self.ref), *inputs, "--out", str(output),
                        *options], check=True, timeout=10)
        return output

    def check_tree(self, out):
        with (out / "support.tsv").open(encoding="utf-8") as stream:
            rows = {int(r["axis_id"]): r for r in csv.DictReader(stream, delimiter="\t")}
        context = json.loads((out / "context.json").read_bytes())
        taxa = context["taxa"]
        axis_rows = {}
        for axis in context["axes"]:
            if not axis["terminal"]:
                bits = int(axis["split_hex"], 16)
                side = frozenset(t for i, t in enumerate(taxa) if (bits >> i) & 1)
                axis_rows[side] = axis_rows[frozenset(taxa) - side] = rows[axis["id"]]

        def label(side, recovered, decidable):
            row = axis_rows[frozenset(side)]
            self.assertEqual((int(row["n_rec"]), int(row["n_dec"])), (recovered, decidable))
            if decidable == 0:
                self.assertEqual(row["support_percent"], "NA")
                return ""
            self.assertAlmostEqual(float(row["support_percent"]), 100.0 * recovered / decidable)
            return row["support_percent"]

        middle = label((TIP, "B", "C", "D"), 2, 3)
        expected = SHAPE.format(ab=label((TIP, "B"), 3, 3), cd=label(("C", "D"), 0, 3),
                                gh=label(("G", "H"), 0, 0), ef=label(("E", "F"), 3, 3),
                                left=middle, right=middle)
        written = (out / "support.nwk").read_bytes()
        self.assertEqual(tokens(written.decode("utf-8")), tokens(expected))
        manifest = json.loads((out / "manifest.json").read_bytes())
        entries = {entry["path"]: entry for entry in manifest["files"]}
        for filename in ("support.tsv", "support.nwk"):
            data = (out / filename).read_bytes()
            self.assertEqual(entries[filename]["bytes"], len(data))
            self.assertEqual(entries[filename]["sha256"], hashlib.sha256(data).hexdigest())

    def test_default_and_all_gene_tree_modes(self):
        self.check_tree(self.run_cli())
        for scope in ("matrix", "full", "support"):
            for role in ("fixed", "free"):
                with self.subTest(scope=scope, role=role):
                    self.check_tree(self.run_cli("--scope", scope, "--role", role))

    def test_binary_display_root_on_terminal_edge(self):
        self.ref.write_text("(A_one:1,('B_two':2,(C:3,D:4):5):6);", encoding="utf-8")
        self.genes.write_text("g\t(A_one,'B_two',(C,D));\n", encoding="utf-8")
        out = self.run_cli()
        # The split opposite A shares its terminal axis. Only CD is internal.
        expected = "(A_one:1,('B_two':2,(C:3,D:4)100:5):6);"
        self.assertEqual((out / "support.nwk").read_text(encoding="utf-8").strip(), expected)

    def test_taxa_only_truth(self):
        out = self.run_cli("--scope", "truth", truth=True)
        self.assertFalse((out / "support.tsv").exists())
        self.assertFalse((out / "support.nwk").exists())


if __name__ == "__main__":
    unittest.main(argv=[sys.argv[0]])
