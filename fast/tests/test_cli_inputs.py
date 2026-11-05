"""Exercise provenance against pipes, replaced paths and byte-sensitive inputs."""
import codecs
import errno
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import time
import unittest


CLI = sys.argv[1]
TREE = b"(A:1,B:2);"


class InputContract(unittest.TestCase):
    def setUp(self):
        self.tmp = TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        self.ref = self.root / "reference"
        self.ref.write_bytes(TREE)
        self.out = self.root / "out"

    def args(self, rows, truth=False, out=None):
        return [CLI, "--reference", str(self.ref), "--taxa" if truth else "--genes", str(rows),
                "--scope", "truth" if truth else "full", "--out", str(out or self.out)]

    def check_digest(self, entry, data):
        self.assertEqual((entry["bytes"], entry["sha256"]),
                         (len(data), hashlib.sha256(data).hexdigest()))

    def start_fifo(self):
        fifo = self.root / "genes.fifo"
        os.mkfifo(fifo)
        proc = subprocess.Popen(self.args(fifo), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.addCleanup(proc.communicate)
        self.addCleanup(proc.kill)
        # Opening the writer proves the CLI has consumed the reference and opened the rows.
        until = time.monotonic() + 10
        while proc.poll() is None and time.monotonic() < until:
            try:
                fd = os.open(fifo, os.O_WRONLY | os.O_NONBLOCK)
                return proc, fifo, os.fdopen(fd, "wb")
            except OSError as exc:
                if exc.errno != errno.ENXIO:
                    raise
                time.sleep(.01)
        self.fail("CLI did not reach the FIFO")

    def test_replacement_records_bytes_from_open_inputs(self):
        proc, fifo, stream = self.start_fifo()
        original = b"\r\nfirst\t" + TREE + b"\r\nsecond\t" + TREE
        with stream:
            self.ref.write_bytes(b"(C,D);")
            replacement = self.root / "replacement"
            replacement.write_bytes(b"never_read\t(C,D);\n")
            replacement.replace(fifo)
            stream.write(original)
        _, err = proc.communicate(timeout=10)
        self.assertEqual(proc.returncode, 0, err)
        manifest = json.loads((self.out / "manifest.json").read_bytes())
        self.assertEqual(manifest["gene_count"], 2)
        self.check_digest(manifest["inputs"][0], TREE)
        self.check_digest(manifest["inputs"][1], original)

    def test_stdin_bom_line_separators_and_final_newline(self):
        self.ref.write_bytes(codecs.BOM_UTF8 + TREE)
        name = "g\u2028x\u2029y"
        rows = codecs.BOM_UTF8 + name.encode() + b"\t" + TREE + b"\r\n\n"
        subprocess.run(self.args("/dev/stdin", out=str(self.out) + "/"),
                       input=rows, check=True, timeout=10)
        manifest = json.loads((self.out / "manifest.json").read_bytes())
        self.check_digest(manifest["inputs"][0], self.ref.read_bytes())
        self.check_digest(manifest["inputs"][1], rows)
        record, = (self.out / "mapping.jsonl").read_text(encoding="utf-8").splitlines()
        self.assertEqual(json.loads(record)["gene_id"], name)

    def test_truth_raw_path_and_bom(self):
        data = codecs.BOM_UTF8 + b"\nfirst\tA\tB"
        raw = os.fsencode(self.root) + b"/taxa-\xff"
        with open(raw, "wb") as stream:
            stream.write(data)
        subprocess.run(self.args(os.fsdecode(raw), truth=True), check=True, timeout=10)
        entry = json.loads((self.out / "manifest.json").read_bytes())["inputs"][1]
        self.assertEqual(entry["path_bytes_hex"], raw.hex())
        self.check_digest(entry, data)
        self.assertIn("\tfirst\t", (self.out / "genes.tsv").read_text(encoding="utf-8"))

    def test_existing_publish_target_is_preserved_and_failure_is_locatable(self):
        proc, fifo, stream = self.start_fifo()
        with stream:
            self.out.mkdir()
            (self.out / "keep").write_bytes(b"another run")
            stream.write(b"gene\t" + TREE)
        _, err = proc.communicate(timeout=10)
        self.assertNotEqual(proc.returncode, 0)
        self.assertEqual([f.name for f in self.out.iterdir()], ["keep"])
        stage, = self.root.glob(".out.staging.*")
        self.assertIn(os.fsencode(stage), err)
        self.assertEqual(json.loads((stage / "FAILED.json").read_bytes())["status"], "FAILED")
        self.assertFalse((stage / "manifest.json").exists())

    def test_missing_output_parent(self):
        dest = self.root / "missing" / "out"
        result = subprocess.run(self.args("/dev/stdin", out=dest), input=b"", capture_output=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn(os.fsencode(dest.parent), result.stderr)
        self.assertIn(b"parent directory", result.stderr)


if __name__ == "__main__":
    unittest.main(argv=[sys.argv[0]])
