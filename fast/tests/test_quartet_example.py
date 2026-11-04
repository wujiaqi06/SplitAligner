"""Check the documented CLI and expected biological results on two quartets."""
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest


class QuartetExample(unittest.TestCase):
    def test_reference_and_changed_internal_edge(self):
        cli = sys.argv[1]
        example = Path(__file__).resolve().parents[1] / "examples" / "quickstart.sh"
        with TemporaryDirectory() as tmp:
            out = Path(tmp) / "quartet"
            subprocess.run([str(example), cli, str(out)], check=True)
            context = json.loads((out / "context.json").read_text())
            internal, = [axis["id"] for axis in context["axes"] if not axis["terminal"]]
            column = "A{}".format(internal)
            with (out / "matrix.tsv").open() as stream:
                same, different = csv.DictReader(stream, delimiter="\t")
            self.assertEqual((same["gene_id"], different["gene_id"]), ("same", "different"))
            self.assertEqual(float(same[column]), 9.0)
            self.assertEqual(different[column], "NA")
            with (out / "raw_states.tsv").open() as stream:
                same, different = csv.DictReader(stream, delimiter="\t")
            self.assertEqual([same["A{}".format(i)] for i in range(5)], ["D"] * 5)
            self.assertEqual(different[column], "U")
            self.assertEqual(list(different.values()).count("D"), 4)
            manifest = json.loads((out / "manifest.json").read_text())
            self.assertEqual(manifest["status"], "COMPLETED")
            self.assertEqual(manifest["gene_count"], 2)
            policy = manifest["unmatched_splits"]
            self.assertEqual((policy["n_reference_absent"], policy["n_gene_unmatched"]), (1, 1))
            self.assertTrue(policy["written"])
            for entry in manifest["files"]:
                data = (out / entry["path"]).read_bytes()
                self.assertEqual((len(data), hashlib.sha256(data).hexdigest()),
                                 (entry["bytes"], entry["sha256"]))


if __name__ == "__main__":
    unittest.main(argv=[sys.argv[0]])
