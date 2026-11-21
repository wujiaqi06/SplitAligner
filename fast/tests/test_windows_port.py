"""Portable CLI checks, including native Windows path and publication boundaries.

Run on every platform. A Mac/Linux pass is regression evidence, not evidence of
native Windows compatibility. Only the optional symlink case may skip.
"""
import codecs
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

CLI = str(Path(sys.argv[1]).resolve())
REFERENCE = "((甲,b),x,(c,d));"
GENES = {
    "全体": "((甲:1,b:2):3,x:4,(c:5,d:6):7);",
    "融合": "((甲:1,b:2):5,(c:3,d:4):6);",
    "零": "((甲:1,b:2):0,(c:3,d:4):0);",
    "缺长度": "((甲:1,b:2),(c:3,d:4):6);",
    "拓扑": "((甲:1,c:2):3,b:4,d:5);",
}


def records(path):
    return [json.loads(line) for line in path.read_bytes().split(b"\n") if line]


class PortableCLI(unittest.TestCase):
    def setUp(self):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        self.root = Path(tmp.name) / "Windows 路径 with spaces"
        self.root.mkdir()
        self.ref = self.root / "参考 树.nwk"
        self.ref_bytes = codecs.BOM_UTF8 + REFERENCE.encode("utf-8") + b"\r\n"
        self.ref.write_bytes(self.ref_bytes)
        self.genes = self.root / "基因 输入.tsv"
        # No final LF: the manifest must hash exact bytes, including CRLF/BOM.
        self.rows = codecs.BOM_UTF8 + "\r\n".join(
            f"{gene}\t{tree}" for gene, tree in GENES.items()).encode("utf-8")
        self.genes.write_bytes(self.rows)

    def invoke(self, name, *options, stdin=None, truth=False, target=None):
        out = target if target is not None else self.root / name
        source = "-" if stdin is not None else str(self.genes)
        argv = [CLI, "--reference", str(self.ref),
                "--taxa" if truth else "--genes", source,
                "--out", str(out), "--scope", "truth" if truth else "full",
                *options]
        run = subprocess.run(argv, input=stdin, capture_output=True, timeout=20)
        return out, run

    def successful(self, name, *options, **kwargs):
        out, run = self.invoke(name, *options, **kwargs)
        self.assertEqual(run.returncode, 0, run.stderr.decode("utf-8", "replace"))
        manifest = json.loads((out / "manifest.json").read_bytes())
        self.assertEqual(manifest["status"], "COMPLETED")
        self.assertFalse(list(out.rglob(".fusion-spool")))
        self.assertEqual(set(p.name for p in out.iterdir()),
                         {entry["path"] for entry in manifest["files"]} | {"manifest.json"})
        for entry in manifest["files"]:
            data = (out / entry["path"]).read_bytes()
            self.assertEqual((entry["bytes"], entry["sha256"]),
                             (len(data), hashlib.sha256(data).hexdigest()))
            self.assertNotIn(b"\r\n", data, entry["path"])
        self.assertEqual(manifest["binary_sha256"], hashlib.sha256(Path(CLI).read_bytes()).hexdigest())
        return out, manifest

    def check_input(self, entry, content):
        self.assertEqual((entry["bytes"], entry["sha256"]),
                         (len(content), hashlib.sha256(content).hexdigest()))

    def test_unicode_crlf_views_preserve_science_and_sources(self):
        whole, wm = self.successful("完整 输出", "--fusion-view", "whole")
        direct, dm = self.successful("直接 输出", "--fusion-view", "direct")
        off, om = self.successful("关闭 输出", "--fusion-view", "off")
        for manifest, view in ((wm, "whole"), (dm, "direct"), (om, "off")):
            self.assertEqual(manifest["fusion_view"], view)
            self.check_input(manifest["inputs"][0], self.ref_bytes)
            self.check_input(manifest["inputs"][1], self.rows)
            self.assertTrue(Path(manifest["inputs"][0]["path"]).samefile(self.ref))
            self.assertTrue(Path(manifest["inputs"][1]["path"]).samefile(self.genes))
            self.assertNotIn("\ufffd", manifest["inputs"][0]["path"])
        common = {entry["path"] for entry in om["files"]}
        for name in common:
            self.assertEqual((whole / name).read_bytes(), (direct / name).read_bytes(), name)
            self.assertEqual((whole / name).read_bytes(), (off / name).read_bytes(), name)
        self.assertFalse((off / "fusion_matrix.tsv").exists())
        self.assertFalse((off / "fusion_cells.jsonl").exists())
        wc = {r["gene_id"]: r for r in records(whole / "fusion_cells.jsonl")}
        dc = {r["gene_id"]: r for r in records(direct / "fusion_cells.jsonl")}
        self.assertEqual(set(wc), set(GENES))
        self.assertEqual({g: r["value"] for g, r in wc.items()},
                         {"全体": 10, "融合": 11, "零": 0, "缺长度": None, "拓扑": None})
        self.assertEqual(wc["全体"]["reference_state"], "decomposable")
        self.assertEqual(dc["全体"]["numeric_status"], "direct_only")
        self.assertIsNone(dc["全体"]["value"])
        # A degree-two display root contributes two raw edges to one normalized
        # observation. Direct therefore retains 5 + 6 = 11, not NA or one side.
        self.assertEqual(dc["融合"]["value"], 11)
        self.assertEqual(len(dc["融合"]["components"]), 1)
        self.assertEqual(len(dc["融合"]["components"][0]["source_nodes"]), 2)
        for gene in GENES:
            self.assertEqual(wc[gene]["components"], dc[gene]["components"])
            self.assertEqual(wc[gene]["reference_state"], dc[gene]["reference_state"])
            self.assertEqual(wc[gene]["recovery_status"], dc[gene]["recovery_status"])

    def test_stdin_genes_and_taxa_hash_exact_bytes(self):
        file_out, fm = self.successful("文件")
        pipe_out, pm = self.successful("标准 输入", stdin=self.rows)
        self.assertEqual(pm["inputs"][1]["path"], "-")
        self.check_input(pm["inputs"][1], self.rows)
        for entry in fm["files"]:
            self.assertEqual((file_out / entry["path"]).read_bytes(),
                             (pipe_out / entry["path"]).read_bytes())
        taxa = codecs.BOM_UTF8 + "真值\t甲\tb\tc\td\r\n".encode("utf-8")
        _, tm = self.successful("真值 输入", stdin=taxa, truth=True)
        self.assertEqual(tm["gene_count"], 1)
        self.assertEqual(tm["inputs"][1]["path"], "-")
        self.check_input(tm["inputs"][1], taxa)

    def test_existing_files_and_directories_are_preserved(self):
        for kind in ("file", "directory"):
            with self.subTest(kind=kind):
                dest = self.root / ("已存在 " + kind)
                keep = dest if kind == "file" else dest / "keep"
                if kind == "directory":
                    dest.mkdir()
                keep.write_bytes(b"leave this untouched")
                _, run = self.invoke("unused", target=dest)
                self.assertNotEqual(run.returncode, 0)
                self.assertIn(b"output already exists", run.stderr)
                self.assertEqual(keep.read_bytes(), b"leave this untouched")
        self.assertFalse(list(self.root.glob("*.staging.*")))

    def test_failures_leave_receipt_without_published_or_leaked_spool(self):
        for injection in ("read", "write", "close"):
            with self.subTest(injection=injection):
                name = "失败 " + injection
                out, run = self.invoke(name, "--inject-io-failure", injection)
                self.assertNotEqual(run.returncode, 0)
                self.assertFalse(out.exists())
                stages = list(self.root.glob("." + name + ".staging.*"))
                self.assertEqual(len(stages), 1)
                self.assertEqual(json.loads((stages[0] / "FAILED.json").read_bytes())["status"], "FAILED")
                self.assertFalse((stages[0] / "manifest.json").exists())
                self.assertFalse((stages[0] / ".fusion-spool").exists())

    def test_dangling_symlink_target_is_preserved_when_available(self):
        dest = self.root / "输出 链接"
        absent = self.root / "never-created"
        try:
            dest.symlink_to(absent, target_is_directory=True)
        except OSError as exc:
            if os.name == "nt" and getattr(exc, "winerror", None) == 1314:
                self.skipTest("Windows symlink privilege unavailable: " + str(exc))
            raise
        _, run = self.invoke("unused", target=dest)
        self.assertNotEqual(run.returncode, 0)
        self.assertTrue(dest.is_symlink())
        self.assertFalse(absent.exists())


if __name__ == "__main__":
    unittest.main(argv=[sys.argv[0]])
