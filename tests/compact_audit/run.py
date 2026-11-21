#!/usr/bin/env python3
"""Bounded audit CLI regression, including honestly rebound semantic faults.

Every case copies one generated three-gene package. Audit is the real command,
not a helper import. Mutation and resealing are confined to those fixture copies;
source trees and the generated seed are never edited. No benchmark is run.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import stat
import subprocess
import sys
import tempfile
import traceback


ROOT = Path(__file__).resolve().parents[2]
LABEL = "panel"
FULL = "z_full"
FUSED = "a_missing"
UNICODE = "中/基因"
GENES = [FULL, FUSED, UNICODE]
TABLES = {"no", "with", "state", "truth", "groups", "cells", "map"}


def source_binding():
    paths = [ROOT / "SplitAligner.pl", ROOT / "scripts/audit_compact_gene.pl",
             *sorted((ROOT / "lib/SplitAligner").glob("*.pm")), Path(__file__).resolve()]
    return {str(path.relative_to(ROOT)): sha(path.read_bytes()) for path in paths}


def sha(data):
    return hashlib.sha256(data).hexdigest()


def dump(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, sort_keys=True, indent=2) + "\n", encoding="utf-8")


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def snapshot(directory):
    result = {}
    for path in sorted(directory.rglob("*")):
        info = path.lstat()
        key = path.relative_to(directory).as_posix()
        base = [info.st_mode, info.st_dev, info.st_ino, info.st_mtime_ns]
        if stat.S_ISREG(info.st_mode):
            result[key] = base + [sha(path.read_bytes())]
        elif stat.S_ISLNK(info.st_mode):
            result[key] = base + [os.readlink(path)]
        else:
            result[key] = base
    return result


RESEAL = r'''
use strict; use warnings;
use File::Spec;
use SplitAligner::Provenance qw(read_json_file write_json_file);
use SplitAligner::OutputOwnership qw(attach_output_ownership);
my ($base, $manifest_path) = @ARGV;
my $manifest = read_json_file($manifest_path);
my $ownership = delete $manifest->{ownership};
my @items = map { +{
    source => $_->{relative_path},
    destination => File::Spec->catfile($base, $_->{relative_path}),
    role => $_->{role}, type => $_->{type}, shared => $_->{shared},
} } @{$ownership->{outputs}};
push @items, {
    source => $ownership->{manifest_relative_path}, destination => $manifest_path,
    role => $ownership->{manifest_role}, type => 'file', owner_manifest => 1,
};
attach_output_ownership($manifest, mode => $ownership->{owner_mode},
    owner_label => $ownership->{owner_label}, source_root => $base,
    owner_manifest_path => $manifest_path, items => \@items);
write_json_file($manifest_path, $manifest);
'''


class Harness:
    def __init__(self, work):
        self.work = work
        self.logs = work / "logs"
        self.logs.mkdir()
        self.results = []
        self.seed = work / "seed"
        self.env = dict(os.environ, LC_ALL="C", PERL_HASH_SEED="731", PERL_PERTURB_KEYS="2")

    def command(self, name, argv, cwd):
        log = self.logs / name
        log.mkdir()
        before = source_binding()
        result = subprocess.run(argv, cwd=cwd, env=self.env, capture_output=True, text=True,
                                encoding="utf-8", errors="strict", timeout=60)
        after = source_binding()
        dump(log / "command.json", {"argv": [str(a) for a in argv], "cwd": str(cwd),
                                    "exit_code": result.returncode, "LC_ALL": self.env["LC_ALL"]})
        dump(log / "source_binding.json", {"before": before, "after": after,
                                           "unchanged_during_command": before == after})
        (log / "stdout.txt").write_text(result.stdout, encoding="utf-8")
        (log / "stderr.txt").write_text(result.stderr, encoding="utf-8")
        require(before == after, name + ": candidate source changed while this command was running")
        return result

    def generate(self):
        (self.seed / "inputs").mkdir(parents=True)
        (self.seed / "inputs/species.nwk").write_text("((A,(B,C)),(D,E));\n", encoding="utf-8")
        # Input order differs from character order; CRLF, blanks, Unicode and no
        # final newline exercise raw-byte locations without adding a BOM.
        rows = ["", FULL + "((A:1,(B:2,C:3):4):5,(D:6,E:7):8);", "",
                FUSED + "((A:1,B:6):5,(D:6,E:7):8);",
                UNICODE + "((A:1,D:6):5,(B:2,E:7):8);"]
        (self.seed / "inputs/genes.nwk").write_bytes("\r\n".join(rows).encode("utf-8"))
        result = self.command("generate", ["perl", str(ROOT / "SplitAligner.pl"), "--mode", "matrix",
            "--species", "inputs/species.nwk", "--gene", "inputs/genes.nwk", "--label", LABEL,
            "--fusion-view", "whole"], self.seed)
        require(result.returncode == 0, "Seed generation failed: " + result.stderr)
        dump(self.work / "seed_binding.json", snapshot(self.seed))

    def copy(self, name):
        directory = self.work / "cases" / name
        shutil.copytree(self.seed, directory)
        return directory

    def audit(self, name, directory, *, genes=None, success=False, diagnostic=None,
              gene_input="inputs/genes.nwk", species_input="inputs/species.nwk", explicit=True):
        selected = genes if genes is not None else [FULL]
        argv = ["perl", str(ROOT / "scripts/audit_compact_gene.pl"),
                "--manifest", str(directory / f"{LABEL}.run_manifest.json")]
        for gene in selected:
            argv += ["--audit-gene", gene]
        if explicit:
            argv += ["--gene-input", str(directory / gene_input),
                     "--species-input", str(directory / species_input)]
        before = snapshot(directory)
        result = self.command(name, argv, directory)
        require(snapshot(directory) == before, name + ": audit modified the source package or inputs")
        if success:
            require(result.returncode == 0, name + ": " + result.stderr)
            payload = json.loads(result.stdout)
            require(payload.get("status") == "PASS", name + ": missing PASS JSON")
            require(sorted(r["gene_id"] for r in payload["genes"]) == sorted(set(selected)),
                    name + ": wrong selected-gene records")
            require(set(payload.get("scientific_tables", [])) == TABLES,
                    name + ": seven scientific tables were not verified")
            require("scientific" in payload.get("scope", "").lower(), name + ": incomplete audit scope")
        else:
            require(result.returncode != 0, name + ": corruption incorrectly passed")
            require(diagnostic is not None and re.search(diagnostic, result.stderr, re.I),
                    name + ": unexpected failure (possibly a hash-only rejection): " + result.stderr)
            require('"status" : "PASS"' not in result.stdout and '"status":"PASS"' not in result.stdout,
                    name + ": printed success despite failure")
        self.results.append({"case": name, "status": "PASS", "audit_exit": result.returncode,
                             "package_unchanged": True, "selected_genes": selected})
        print("PASS", name, flush=True)

    def paths(self, directory):
        manifest_path = directory / f"{LABEL}.run_manifest.json"
        meta = read_json(manifest_path)
        files = {k: directory / v["filename"] for k, v in meta["compact_provenance"]["files"].items()}
        return manifest_path, meta, files

    def rebind(self, name, directory):
        manifest_path, meta, files = self.paths(directory)
        index = read_json(files["index"])
        index["ledger_sha256"] = sha(files["ledger"].read_bytes())
        index["dictionary_sha256"] = sha(files["dictionary"].read_bytes())
        dump(files["index"], index)

        def refresh(node):
            if isinstance(node, dict):
                if "filename" in node and "sha256" in node:
                    path = directory / node["filename"]
                    require(path.resolve().is_relative_to(directory.resolve()), "Unsafe fixture path")
                    require(path.is_file(), "Missing fixture artifact " + str(path))
                    node["sha256"] = sha(path.read_bytes())
                for value in node.values():
                    refresh(value)
            elif isinstance(node, list):
                for value in node:
                    refresh(value)
        refresh(meta)
        meta["coordinate_state"]["matrix_with_fuse_sha256"] = meta["outputs"]["matrix_with_fuse"]["sha256"]
        dump(manifest_path, meta)
        result = self.command(name + "_reseal", ["perl", "-I" + str(ROOT / "lib"), "-e", RESEAL,
                                                str(directory), str(manifest_path)], directory)
        require(result.returncode == 0, "Fixture reseal failed: " + result.stderr)

    def mutate_record(self, name, directory, gene, mutate):
        _, _, files = self.paths(directory)
        records = [json.loads(raw) for raw in files["ledger"].read_bytes().splitlines()]
        target = next(r for r in records if r["gene_id"] == gene)
        mutate(target)
        index = read_json(files["index"])
        by_gene = {e["gene_id"]: e for e in index["genes"]}
        chunks = []
        offset = 0
        for record in records:
            raw = (json.dumps(record, ensure_ascii=False, sort_keys=True, separators=(",", ":")) + "\n").encode("utf-8")
            entry = by_gene[record["gene_id"]]
            entry.update(ledger_offset=offset, ledger_bytes=len(raw), sha256=sha(raw),
                         observations=len(record["observations"]), components=len(record["components"]))
            chunks.append(raw)
            offset += len(raw)
        files["ledger"].write_bytes(b"".join(chunks))
        dump(files["index"], index)
        self.rebind(name, directory)

    def run(self):
        self.generate()
        self.audit("valid_single", self.copy("valid_single"), success=True)
        self.audit("valid_multi_unicode", self.copy("valid_multi_unicode"), genes=[UNICODE, FULL, FUSED], success=True)
        self.audit("valid_defaults", self.copy("valid_defaults"), genes=[FUSED], success=True, explicit=False)
        p = self.copy("valid_relocated_inputs")
        (p / "relocated").mkdir()
        for filename in ("genes.nwk", "species.nwk"):
            (p / "inputs" / filename).rename(p / "relocated" / filename)
        self.audit("valid_relocated_inputs", p, success=True, genes=[UNICODE],
                   gene_input="relocated/genes.nwk", species_input="relocated/species.nwk")
        for kind, filename in [("gene", "genes.nwk"), ("species", "species.nwk")]:
            name = "changed_" + kind + "_input"
            p = self.copy(name)
            input_path = p / "inputs" / filename
            input_path.write_bytes(input_path.read_bytes() + b"\n")
            self.audit(name, p, diagnostic="hash differs")
            name = "missing_" + kind + "_input"
            p = self.copy(name)
            (p / "inputs" / filename).unlink()
            self.audit(name, p, diagnostic="Cannot access original audit input")
        name = "ledger_unbound_tamper"
        p = self.copy(name)
        _, _, files = self.paths(p)
        original = files["ledger"].read_bytes()
        altered = original.replace(FULL.encode("utf-8"), b"y_full", 1)
        require(original != altered, "Could not mutate a saved ledger block")
        files["ledger"].write_bytes(altered)
        self.audit(name, p, diagnostic="hash differs")

        name = "index_coordinates_rebound"
        p = self.copy(name)
        _, _, files = self.paths(p)
        index = read_json(files["index"])
        next(e for e in index["genes"] if e["gene_id"] == FULL)["input_offset"] += 1
        dump(files["index"], index)
        self.rebind(name, p)
        self.audit(name, p, diagnostic="Saved input coordinates differ")

        name = "ledger_number_rebound"
        p = self.copy(name)
        def number(record):
            obs = next(o for o in record["observations"] if o["length"] != "NA")
            obs["length"] = "987.25"
        self.mutate_record(name, p, FULL, number)
        self.audit(name, p, diagnostic=r"Audit mismatch.*field 'observations'")

        name = "ledger_members_rebound"
        p = self.copy(name)
        def member(record):
            first, second = record["components"][:2]
            first["members"], second["members"] = second["members"], first["members"]
        self.mutate_record(name, p, FUSED, member)
        self.audit(name, p, genes=[FUSED], diagnostic=r"Audit mismatch.*field 'components'")

        name = "ledger_sources_rebound"
        p = self.copy(name)
        def source(record):
            a, b = record["observations"][:2]
            a["source_nodes"], b["source_nodes"] = b["source_nodes"], a["source_nodes"]
        self.mutate_record(name, p, FULL, source)
        self.audit(name, p, diagnostic=r"Audit mismatch.*field 'observations'")

        table_paths = {
            "no": f"{LABEL}.matrix_no_fuse.txt", "with": f"{LABEL}.matrix_with_fuse.txt",
            "state": f"{LABEL}.primitive_state.tsv", "map": f"{LABEL}.gene_id_map.tsv",
            "truth": f"{LABEL}_split_branch_label/{LABEL}.primitive_truth.tsv",
            "groups": f"{LABEL}_split_branch_label/{LABEL}.projected_groups.tsv",
            "cells": f"{LABEL}.fusion_cells.tsv",
        }
        for kind, relative in table_paths.items():
            name = "scientific_" + kind + "_rebound"
            p = self.copy(name)
            path = p / relative
            rows = [line.split("\t") for line in path.read_text(encoding="utf-8").splitlines()]
            selected = FUSED if kind == "groups" else FULL
            gene_column = 1 if kind == "map" else 0
            row = next(r for r in rows[1:] if r[gene_column] == selected)
            if kind == "map":
                row[0] = "g99999999"
            elif kind == "state":
                row[1] = "S" if row[1] != "S" else "D"
            elif kind == "truth":
                row[1] = "NA_struct" if row[1] != "NA_struct" else "Mapped"
            elif kind in {"groups", "cells"}:
                row[-1] = "987.25" if row[-1] != "987.25" else "123.5"
            else:
                row[1] = "987.25" if row[1] != "987.25" else "123.5"
            path.write_text("".join("\t".join(r) + "\n" for r in rows), encoding="utf-8")
            self.rebind(name, p)
            self.audit(name, p, genes=[selected], diagnostic="Audit scientific row mismatch")
        require(snapshot(self.seed) == read_json(self.work / "seed_binding.json"), "Seed package was altered")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work-dir", type=Path, help="New evidence directory, kept for review")
    args = parser.parse_args()
    work = args.work_dir.resolve() if args.work_dir else Path(tempfile.mkdtemp(prefix="splitaligner-audit-test-"))
    if work.exists():
        require(not any(work.iterdir()), "Refusing to overwrite evidence directory " + str(work))
    work.mkdir(parents=True, exist_ok=True)
    h = Harness(work)
    binding = source_binding()
    dump(work / "source_binding.json", binding)
    result = {"scope": "Real audit CLI; three-gene fixture; no benchmark or publication acceptance",
              "source": str(ROOT), "cases": h.results, "status": "RUNNING"}
    exit_code = 0
    try:
        h.run()
        result["status"] = "PASS"
    except Exception as exc:
        result["status"] = "FAIL"
        result["failure"] = str(exc)
        (work / "failure.txt").write_text(traceback.format_exc(), encoding="utf-8")
        print(traceback.format_exc(), file=sys.stderr)
        exit_code = 1
    result["passed"] = len(h.results)
    result["source_unchanged_across_suite"] = source_binding() == binding
    if result["status"] == "PASS" and not result["source_unchanged_across_suite"]:
        result["status"] = "FAIL"
        result["failure"] = "Candidate source changed across suite; per-command bindings retained."
        exit_code = 1
    dump(work / "RESULTS.json", result)
    print(json.dumps({"status": result["status"], "passed": len(h.results), "evidence": str(work)}, ensure_ascii=False))
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
