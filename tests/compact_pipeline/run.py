#!/usr/bin/env python3
"""Independent, bounded v1/v2 matrix and provenance acceptance harness.

This script never edits either source tree.  --prepare-only runs the frozen
controller on small fixtures; a later invocation uses those pinned outputs.
All subprocess working directories and logs belong to --work-dir.  Existing
baseline evidence is verified, never silently regenerated.
"""
import argparse
import ast
import base64
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import re
import shutil
import stat
import subprocess
import sys
import traceback


SCHEMA = "SplitAligner-compact-cross-version-harness-v1"
SEED = {"PERL_HASH_SEED": "1", "PERL_PERTURB_KEYS": "2"}
SEED_VARIANTS = ({"PERL_HASH_SEED": "0", "PERL_PERTURB_KEYS": "0"},
                 {"PERL_HASH_SEED": "7351", "PERL_PERTURB_KEYS": "1"})
COMMON = ["species_tree.forSplit.nwk", "species_tree.FigTree.tre",
          "species_tree.splits.txt", "species_tree.branch_map.txt",
          "species_tree.primitive_axis.tsv"]
SUFFIXES = ["matrix_no_fuse.txt", "matrix_with_fuse.txt", "fusion_cells.tsv",
            "primitive_state.tsv", "gene_id_map.tsv", "primitive_axis.tsv"]
AGGREGATES = ["primitive_truth.tsv", "projected_groups.tsv", "branch_patterns.txt"]


def sha(data):
    return hashlib.sha256(data).hexdigest()


def dump(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8")


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def encoded(data):
    return base64.b64encode(data if isinstance(data, bytes) else data.encode("utf-8")).decode("ascii")


def decoded(data):
    return base64.b64decode(data, validate=True)


def literal(source, name):
    """Read fixture constants without executing a baseline test module."""
    for node in ast.parse(source.read_text(encoding="utf-8")).body:
        if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == name for t in node.targets):
            return ast.literal_eval(node.value)
    raise AssertionError(f"Missing literal {name} in {source}")


def tree_input(genes, ending="\n", blanks=False, final_newline=True):
    rows = [gene + tree for gene, tree in genes.items()]
    if blanks:
        rows = [""] + [part for row in rows for part in (row, "")]
        rows.pop()
    return ending.join(rows) + (ending if final_newline else "")


def plan(frozen):
    tests = frozen / "tests"
    whole = tests / "whole_component_fusion_regression/run.py"
    priority = tests / "direct_view_priority_regression/run.py"
    figure = tests / "fusion_topology_regression/run.py"
    cases = []

    def positive(name, ref, genes, view="whole", paired=False, seeds=False):
        panels = genes if paired else {"x": genes}
        cases.append({"name": name, "species": encoded(ref + "\n"),
                      "panels": {label: encoded(tree_input(trees)) for label, trees in panels.items()},
                      "view": view, "success": True, "paired": paired, "seeds": seeds})

    for view in ("whole", "direct"):
        positive("whole_components_" + view, literal(whole, "REF"), literal(whole, "GENES"),
                 view, seeds=view == "whole")
        positive("priority_" + view, literal(priority, "REFERENCE"), literal(priority, "GENES"), view)
        fg = literal(figure, "FIGURE4")
        genes = {**fg, "zero": "((A,B):0,D,(E,F));", "no_length": "((A,B),D,(E,F));",
                 "star": "(A,B,D,E,F);", "pooled": "((A,B):1,(D,(E,F)):2);",
                 "mixed_root": "((A,B):1,(D,(E,F)));", "unique_topo": "((A,D),(B,C),(E,F));"}
        fixed = {gene: fg["g1"] for gene in genes}
        fixed["unique_topo"] = "((A,B):1,C,(D,(E,F)));"
        positive("figure4_" + view, "((A,B),C,(D,(E,F)));",
                 {"alpha": {**genes, "free_only": fg["g3"]}, "beta": fixed}, view, paired=True)

    numeric = {
        "z_token": "((A:1.,B:+1e0):-2.5E-1,(C:.5,D:0.5):3.);",
        "a_zero": "((A:-0,B:+0e999):2,(C:1,D:1):3);",
        "subnormal": "((A:5e-324,B:1):2,(C:1,D:1):3);",
        "negative": "((A:-0.1,B:.2):2,(C:1,D:1):3);",
        "root_missing": "((A:1,B:1):2,(C:1,D:1));",
        "star": "(A:1,B:2,C:3,D:4);",
        "two_taxa": "(A:2,C:3);",
        "three_taxa": "(A:1,B:2,C:3);",
        "g/a": "((A:1,B:2):3,(C:4,D:5):6);",
        "g_a": "((A:7,C:8):9,(B:10,D:11):12);",
    }
    positive("numeric_shapes", "((A,B),(C,D));", numeric)
    unicode_genes = {"末尾/基因": "((猫:1,A|B:2):3,(.C:4,雪:5):6);",
                     "alpha_é": "(猫:1,A|B:2,雪:3);",
                     "g/a": "((猫:0,A|B:0):0,(.C:1,雪:2):0);",
                     "g_a": "((猫:1,.C:2):3,(A|B:4,雪:5):6);"}
    positive("unicode_crlf", "((猫,A|B),(.C,雪));", unicode_genes, seeds=True)
    cases[-1]["species"] = encoded("((猫,A|B),(.C,雪));\r\n")
    cases[-1]["panels"]["x"] = encoded(tree_input(unicode_genes, "\r\n", blanks=True, final_newline=False))
    reversed_genes = dict(reversed(list(unicode_genes.items())))
    positive("unicode_reversed", "((猫,A|B),(.C,雪));", reversed_genes)
    positive("numeric_taxa", "((1,2),(3,4));", {"gene9": "((1:1,2:2):5,(3:3,4:4):6);"})
    for name, taxon in [("double_dot", "A..B"), ("double_pipe", "A||B"),
                        ("trailing_dot", "NC.045."), ("ordinary_period", "NC_045512.2")]:
        positive("codec_" + name, f"((A,B),({taxon},C));",
                 {"g": f"((A:1,B:2):10,({taxon}:3,C:4):20);"})

    normal = "g((A:1,B:2):3,(C:4,D:5):6);\n"
    rejected = [
        ("empty", b"", "No non-empty|Empty"),
        ("duplicate_gene", normal + normal, "Duplicate gene ID"),
        ("empty_gene", "((A,B),(C,D));\n", "Empty gene ID"),
        ("unknown", "g((A,B),(C,Z));\n", "Unknown taxon|absent from the reference"),
        ("duplicate_taxon", "g((A,A),(C,D));\n", "Duplicate terminal taxon"),
        ("unbalanced", "g((A,B),(C,D);\n", "Expected|Unexpected"),
        ("two_trees", "g((A,B),(C,D));((A,C),(B,D));\n", "Extra content after"),
        ("quoted", "g(('A',B),(C,D));\n", "Quoted labels"),
        ("comments", "g((A[comment],B),(C,D));\n", "Bracket comments"),
        ("unary", "g(((A),B),(C,D));\n", "at least two children"),
        ("tab_gene", "g\tbad((A,B),(C,D));\n", "Gene ID contains a tab"),
        ("bom_gene", b"\xef\xbb\xbf" + normal.encode(), "BOM"),
        ("invalid_utf8", b"g\xff((A,B),(C,D));\n", "Invalid UTF-8"),
        ("bad_value", "g((A:bad,B:1),(C:1,D:1));\n", "Invalid branch value"),
        ("space_taxon", "g((A B,B),(C,D));\n", "Expected|Unexpected|at least two children"),
        ("overflow", "g((A:1e999,B:1):2,(C:1,D:1):3);\n", "not finite"),
        ("underflow", "g((A:1e-999,B:1):2,(C:1,D:1):3);\n", "underflows to zero"),
        ("unicode_digit", "g((A:١,B:1):2,(C:1,D:1):3);\n", "Invalid branch value"),
    ]
    for name, data, diagnostic in rejected:
        cases.append({"name": "reject_" + name, "species": encoded("((A,B),(C,D));\n"),
                      "panels": {"x": encoded(data)}, "view": "whole", "success": False,
                      "diagnostic": diagnostic, "paired": False, "seeds": False})
    for view in ("whole", "direct"):
        cases.append({"name": "reject_pooled_" + view, "species": encoded("((A,B),(C,D));\n"),
                      "panels": {"x": encoded("g((A:1,B:1):1e308,(C:1,D:1):1e308);\n")},
                      "view": view, "success": False, "diagnostic": "Pooled unrooted branch",
                      "paired": False, "seeds": False})
    cases.append({"name": "reject_bom_species", "species": encoded(b"\xef\xbb\xbf((A,B),(C,D));\n"),
                  "panels": {"x": encoded(normal)}, "view": "whole", "success": False,
                  "diagnostic": "BOM", "paired": False, "seeds": False})
    return {"schema": SCHEMA, "fixture_sources": {str(p.relative_to(frozen)): sha(p.read_bytes())
            for p in (whole, priority, figure)}, "cases": cases}


def source_identity(root):
    paths = [root / "SplitAligner.pl"]
    paths += sorted((root / "lib").rglob("*.pm")) + sorted((root / "scripts").glob("*.pl"))
    return {str(p.relative_to(root)): sha(p.read_bytes()) for p in paths}


def table(path):
    return list(csv.DictReader(io.StringIO(path.read_text(encoding="utf-8")), delimiter="\t"))


def scientific(label):
    return COMMON + [f"{label}.{suffix}" for suffix in SUFFIXES] + [
        f"{label}_split_branch_label/{label}.{suffix}" for suffix in AGGREGATES]


def snapshot(directory):
    result = {}
    for p in sorted(directory.rglob("*")):
        s = p.lstat()
        result[str(p.relative_to(directory))] = {"mode": stat.S_IMODE(s.st_mode),
            "type": stat.S_IFMT(s.st_mode), "inode": s.st_ino,
            "data": os.readlink(p) if p.is_symlink() else sha(p.read_bytes()) if p.is_file() else None}
    return result


class Harness:
    def __init__(self, args):
        self.args = args
        self.work = args.work_dir.resolve()
        self.frozen = args.frozen.resolve()
        self.candidate = args.candidate.resolve() if args.candidate else None
        self.checks = []
        self.commands = 0
        self.log_root = self.work / ("prepare-commands" if args.prepare_only else args.run_name + "-commands")
        self.log_root.mkdir(parents=True, exist_ok=False)

    def check(self, condition, name, detail=None):
        self.checks.append({"name": name, "pass": bool(condition), "detail": detail})
        require(condition, f"{name}: {detail}")

    def call(self, driver_root, cwd, args, seed=SEED):
        self.commands += 1
        stem = self.log_root / f"{self.commands:04d}"
        command = ["perl", str(driver_root / "SplitAligner.pl"), *args]
        env = {**os.environ, "LC_ALL": "C", "LANG": "C", **seed}
        # Prevent ambient Perl startup hooks from injecting behavior into either implementation.
        for key in ("PERL5OPT", "PERL5LIB", "PERLLIB"):
            env.pop(key, None)
        dump(stem.with_suffix(".json"), {"command": command, "cwd": str(cwd), "environment": seed})
        result = subprocess.run(command, cwd=cwd, env=env, capture_output=True, timeout=120)
        stem.with_suffix(".stdout").write_bytes(result.stdout)
        stem.with_suffix(".stderr").write_bytes(result.stderr)
        dump(stem.with_suffix(".exit.json"), {"returncode": result.returncode})
        return result

    def run_case(self, root, source, case, seed=SEED):
        root.mkdir(parents=True, exist_ok=False)
        (root / "reference.nwk").write_bytes(decoded(case["species"]))
        for label, data in case["panels"].items():
            (root / (label + ".input")).write_bytes(decoded(data))
        before = snapshot(root)
        for label in case["panels"]:
            result = self.call(source, root, ["--mode", "matrix", "--species", "reference.nwk",
                               "--gene", label + ".input", "--label", label,
                               "--fusion-view", case["view"]], seed)
            self.check((result.returncode == 0) == case["success"], case["name"] + ":exit:" + label,
                       result.stderr.decode("utf-8", "replace")[-3000:])
            if not case["success"]:
                self.check(bool(re.search(case["diagnostic"], result.stderr.decode("utf-8", "replace"))),
                           case["name"] + ":diagnostic")
                if case["name"] == "reject_unknown":
                    error = result.stderr.decode("utf-8", "replace")
                    self.check("'Z'" in error and re.search(r"line\s+1\s+gene\s+g\b", error),
                               case["name"] + ":original-taxon-gene-and-line")
                self.check(snapshot(root) == before, case["name"] + ":no-publication-or-residue")
            else:
                self.check(not list(root.glob(".splitaligner-*")), case["name"] + ":no-staging-residue")
        return root

    def prepare(self):
        self.check(not (self.work / "baseline.json").exists(), "baseline-does-not-already-exist")
        p = plan(self.frozen)
        dump(self.work / "plan.json", p)
        identity = source_identity(self.frozen)
        expected = {"schema": SCHEMA, "frozen_source": str(self.frozen), "source_identity": identity,
                    "plan_sha256": sha((self.work / "plan.json").read_bytes()), "files": {}}
        for case in p["cases"]:
            print("Preparing", case["name"], flush=True)
            d = self.run_case(self.work / "baseline" / case["name"], self.frozen, case)
            if case["success"]:
                for label in case["panels"]:
                    for name in scientific(label):
                        expected["files"][str((d / name).relative_to(self.work))] = sha((d / name).read_bytes())
                if case["paired"]:
                    for mode in ("finalize", "finalize_fix"):
                        final = self.finalize(self.work / "baseline-final" / case["name"] / mode, d, mode)
                        for name in final:
                            path = self.work / "baseline-final" / case["name"] / mode / name
                            expected["files"][str(path.relative_to(self.work))] = sha(path.read_bytes())
                # Pin old component and observation evidence as well as scientific files.
                for label in case["panels"]:
                    for directory in (d / (label + "_splits"), d / (label + "_split_branch_label")):
                        for file in directory.iterdir():
                            if file.is_file():
                                expected["files"][str(file.relative_to(self.work))] = sha(file.read_bytes())
        self.check(source_identity(self.frozen) == identity, "frozen-source-unchanged-after-prepare")
        dump(self.work / "baseline.json", expected)
        dump(self.work / "prepare-results.json", {"status": "PASS", "commands": self.commands, "checks": self.checks})

    def finalize(self, cwd, matrix_dir, mode):
        cwd.mkdir(parents=True, exist_ok=False)
        # The old CLI publishes *.na_fuse.txt beside the supplied matrices.
        # Give each mode an intact, byte-identical result-directory copy so
        # paired/fix-only outputs cannot collide or modify the pinned baseline.
        inputs = cwd / "inputs"
        shutil.copytree(matrix_dir, inputs)
        args = ["--mode", mode, "--fix", str(inputs / "beta.matrix_with_fuse.txt"), "--final_label", "final"]
        if mode == "finalize":
            args += ["--free", str(inputs / "alpha.matrix_with_fuse.txt"),
                     "--species_tree", str(inputs / "species_tree.forSplit.nwk")]
        # Deliberately the untouched old controller, even when matrix_dir is v2.
        result = self.call(self.frozen, cwd, args)
        self.check(result.returncode == 0, str(cwd.relative_to(self.work)) + ":old-controller-" + mode,
                   result.stderr.decode("utf-8", "replace")[-4000:])
        meta = read_json(cwd / "final.finalize_manifest.json")
        return sorted(final_relative(name) for name in meta["outputs"])

    def compare_files(self, old, new, names, prefix):
        for name in names:
            self.check((new / name).is_file(), prefix + ":exists:" + name)
            self.check((old / name).read_bytes() == (new / name).read_bytes(), prefix + ":bytes:" + name,
                       {"baseline": sha((old / name).read_bytes()), "candidate": sha((new / name).read_bytes())})

    def validate_ledger(self, old, new, label, prefix):
        dictionary = read_json(new / (label + ".taxon_dictionary.json"))
        index = read_json(new / (label + ".compact.index.json"))
        ledger_raw = (new / (label + ".compact.jsonl")).read_bytes()
        input_raw = (new / (label + ".input")).read_bytes()
        self.check(index["schema_version"] == "SplitAligner-compact-ledger-v1", prefix + ":ledger-schema")
        self.check(index["ledger_file"] == label + ".compact.jsonl", prefix + ":ledger-path")
        self.check(index["ledger_sha256"] == sha(ledger_raw), prefix + ":ledger-hash")
        self.check(index["dictionary_file"] == label + ".taxon_dictionary.json", prefix + ":dictionary-path")
        self.check(index["dictionary_sha256"] == sha((new / index["dictionary_file"]).read_bytes()), prefix + ":dictionary-hash")
        reference_raw = (new / "reference.nwk").read_bytes()
        self.check(index["input_sha256"] == {"gene_trees": sha(input_raw), "species_tree": sha(reference_raw)} and
                   index["input_bytes"] == {"gene_trees": len(input_raw), "species_tree": len(reference_raw)},
                   prefix + ":input-binding")
        # The frozen axis ledger is headerless, unlike the scientific matrices.
        source_axis = [dict(zip(("branch_id", "canonical_split"), line.split("\t")))
                       for line in (old / (label + ".primitive_axis.tsv")).read_text(encoding="utf-8").splitlines()]
        first_split = source_axis[0]["canonical_split"]
        left, right = split_sides(first_split)
        taxa = sorted(left | right)
        self.check(dictionary["taxa"] == taxa, prefix + ":dictionary-order")
        axis_digest = read_json(old / (label + ".run_manifest.json"))["axis"]["sha256"]
        self.check(dictionary["axis_sha256"] == index["axis_sha256"] == axis_digest, prefix + ":axis-identity-binding")
        self.check(dictionary["bit_order"] == "vec-LSB-first" and dictionary["bitset_bytes"] == (len(taxa) + 7) // 8,
                   prefix + ":dictionary-bit-contract")
        gene_map = table(old / (label + ".gene_id_map.tsv"))
        self.check([e["gene_id"] for e in index["genes"]] == [g["gene_id"] for g in gene_map], prefix + ":sorted-index")
        all_records = [json.loads(line) for line in ledger_raw.splitlines()]
        self.check(len(all_records) == len(gene_map), prefix + ":one-ledger-record-per-gene")
        input_lines = input_raw.splitlines(keepends=True)
        input_positions = []
        offset = 0
        for number, line in enumerate(input_lines, 1):
            text = line.decode("utf-8").rstrip("\r\n")
            if text.strip():
                gene = text[:text.index("(")].strip()
                input_positions.append((gene, number, offset, len(line)))
            offset += len(line)
        self.check([r["gene_id"] for r in all_records] == [r[0] for r in input_positions], prefix + ":ledger-input-order")
        positions = {gene: (i + 1, line, offset, count) for i, (gene, line, offset, count) in enumerate(input_positions)}
        truth = {r["gene"]: r for r in table(old / f"{label}_split_branch_label/{label}.primitive_truth.tsv")}
        states = {r["gene"]: r for r in table(old / f"{label}.primitive_state.tsv")}
        primitive = {r["gene"]: r for r in table(old / f"{label}.matrix_no_fuse.txt")}
        axis_ids = [row["branch_id"] for row in source_axis]
        covered = []
        for entry, identity in zip(index["genes"], gene_map):
            gene, key = identity["gene_id"], identity["storage_key"]
            detail = prefix + ":gene:" + gene
            self.check(entry["storage_key"] == key, detail + ":storage-key")
            self.check((entry["input_record_id"], entry["input_line"], entry["input_offset"], entry["input_bytes"]) == positions[gene],
                       detail + ":input-byte-location")
            start, count = entry["ledger_offset"], entry["ledger_bytes"]
            raw = ledger_raw[start:start + count]
            covered.append((start, start + count))
            self.check(raw.endswith(b"\n") and raw.count(b"\n") == 1 and sha(raw) == entry["sha256"], detail + ":ledger-block")
            record = json.loads(raw)
            self.check(record["schema"] == "splitaligner-compact-gene-v1", detail + ":record-schema")
            self.check((record["gene_id"], record["input_record_id"], record["input_line"]) == (gene, positions[gene][0], positions[gene][1]),
                       detail + ":stable-identity")
            old_splits = [line.split("\t") for line in (old / f"{label}_splits/{key}.split.txt").read_text().splitlines()]
            old_nodes = [[int(n) for n in line.split(",")] for line in
                         (old / f"{label}_splits/{key}.source_nodes.txt").read_text().splitlines()]
            self.check(len(old_splits) == len(old_nodes), detail + ":baseline-source-row-count")
            expected = {}
            for (split, length), nodes in zip(old_splits, old_nodes):
                left, right = split_sides(split)
                expected[bit_key(left, right, taxa)] = {"length": length, "source_nodes": nodes}
            observations = record["observations"]
            actual = {r["key"]: {"length": r["length"], "source_nodes": r["source_nodes"]} for r in observations}
            self.check(len(actual) == len(observations) == entry["observations"], detail + ":unique-observation-keys")
            self.check(actual == expected, detail + ":all-observations-and-source-edges", {"expected": expected, "actual": actual} if actual != expected else None)
            left, right = split_sides(old_splits[0][0])
            self.check(record["tg"] == set_bits(left | right, taxa), detail + ":tg")
            nodes = sorted(node for observation in observations for node in observation["source_nodes"])
            self.check(nodes == list(range(1, record["source_node_count"] + 1)), detail + ":each-input-edge-exactly-once")
            old_components = table(old / f"{label}_split_branch_label/{key}.components.tsv")
            expected_components = {}
            for c in old_components:
                row = int(c["matched_gene_split_row"]) - 1 if c["matched_gene_split_row"] else None
                split = old_splits[row][0] if row is not None else c["unmatched_projected_split"]
                a, b = split_sides(split)
                expected_components[c["members"]] = {"key": bit_key(a, b, taxa), "recovered": row is not None,
                                                       "length": c["length"]}
            components = record["components"]
            actual_components = {}
            for component in components:
                oid = component["observation_id"]
                self.check(oid is None or type(oid) is int and 1 <= oid <= len(observations), detail + ":observation-reference")
                obs = observations[oid - 1] if oid is not None else None
                if obs is not None:
                    self.check(obs["key"] == component["key"], detail + ":matched-key")
                actual_components[component["members"]] = {"key": component["key"], "recovered": obs is not None,
                                                           "length": obs["length"] if obs is not None else "NA"}
            self.check(len(actual_components) == len(components) == entry["components"], detail + ":unique-complete-components")
            self.check(actual_components == expected_components, detail + ":component-evidence", {"expected": expected_components, "actual": actual_components} if actual_components != expected_components else None)
            for field, source in [("truth", truth), ("state", states), ("primitive", primitive)]:
                self.check(record[field] == [source[gene][branch] for branch in axis_ids], detail + ":" + field)
        covered.sort()
        self.check(covered[0][0] == 0 and covered[-1][1] == len(ledger_raw) and
                   all(a[1] == b[0] for a, b in zip(covered, covered[1:])), prefix + ":index-covers-ledger-exactly")
        self.check(not (new / (label + "_splits")).exists(), prefix + ":old-split-directory-absent")
        self.check(not list((new / (label + "_split_branch_label")).glob("g*.components.tsv")) and
                   not list((new / (label + "_split_branch_label")).glob("g*.split.txt")), prefix + ":old-per-gene-text-absent")

    def execute(self):
        expected = read_json(self.work / "baseline.json")
        self.check(expected["schema"] == SCHEMA, "baseline-schema")
        self.check(source_identity(self.frozen) == expected["source_identity"], "frozen-source-still-bound")
        self.check(sha((self.work / "plan.json").read_bytes()) == expected["plan_sha256"], "fixture-plan-still-bound")
        for path, digest in expected["files"].items():
            self.check(sha((self.work / path).read_bytes()) == digest, "baseline-evidence:" + path)
        p = read_json(self.work / "plan.json")
        initial_candidate = source_identity(self.candidate)
        out = self.work / self.args.run_name
        out.mkdir(exist_ok=False)
        for case in p["cases"]:
            print("Checking", case["name"], flush=True)
            base = self.work / "baseline" / case["name"]
            d = self.run_case(out / case["name"], self.candidate, case)
            if not case["success"]:
                continue
            for label in case["panels"]:
                prefix = case["name"] + ":" + label
                self.compare_files(base, d, scientific(label), prefix)
                self.validate_ledger(base, d, label, prefix)
                baseline_meta = read_json(base / (label + ".run_manifest.json"))
                candidate_meta = read_json(d / (label + ".run_manifest.json"))
                for field in ("axis", "gene_identity", "coordinate_state", "projection_graph"):
                    self.check(candidate_meta[field] == baseline_meta[field], prefix + ":manifest-science:" + field)
                self.check(candidate_meta["fusion_matrix"]["view"] == case["view"], prefix + ":manifest-view")
                self.check(candidate_meta["fusion_matrix"]["schema_version"] != "whole-components-v1", prefix + ":new-storage-schema")
            if case["paired"]:
                for mode in ("finalize", "finalize_fix"):
                    target = out / (case["name"] + "-old-" + mode)
                    names = self.finalize(target, d, mode)
                    baseline_final = self.work / "baseline-final" / case["name"] / mode
                    self.check(names == sorted(final_relative(name) for name in read_json(baseline_final / "final.finalize_manifest.json")["outputs"]),
                               case["name"] + ":old-" + mode + ":output-set")
                    self.compare_files(baseline_final, target, names, case["name"] + ":old-" + mode)
                    for field in ("free_projection_accounting", "primitive_fusion_source"):
                        old_meta = read_json(baseline_final / "final.finalize_manifest.json")
                        if field in old_meta:
                            self.check(read_json(target / "final.finalize_manifest.json")[field] == old_meta[field],
                                       case["name"] + ":old-" + mode + ":" + field)
            if case["seeds"]:
                for i, seed in enumerate(SEED_VARIANTS):
                    alt = self.run_case(out / (case["name"] + f"-seed-{i}"), self.candidate, case, seed)
                    for label in case["panels"]:
                        self.compare_files(d, alt, scientific(label) + [label + ".compact.jsonl", label + ".compact.index.json", label + ".taxon_dictionary.json"],
                                           case["name"] + f":seed-{i}")
                        self.validate_ledger(base, alt, label, case["name"] + f":seed-{i}:" + label)
        # Independent assertions shared by whole and direct must still agree.
        for name in ("whole_components", "priority"):
            whole, direct = out / (name + "_whole"), out / (name + "_direct")
            self.compare_files(whole, direct, ["x.compact.jsonl", "x.matrix_no_fuse.txt", "x.primitive_state.tsv",
                "x_split_branch_label/x.primitive_truth.tsv", "x_split_branch_label/x.projected_groups.tsv"], name + ":view-invariant-evidence")
        # Record layout follows input order; scientific row order follows gene identity.
        self.compare_files(out / "unicode_crlf", out / "unicode_reversed",
                           [n for n in scientific("x") if n != "x.gene_id_map.tsv"], "gene-order-metamorphism")
        self.force_boundary(out)
        self.check(source_identity(self.frozen) == expected["source_identity"], "frozen-source-unchanged-after-tests")
        self.check(source_identity(self.candidate) == initial_candidate, "candidate-source-stable-during-tests")
        dump(self.work / (self.args.run_name + "-results.json"), {"status": "PASS", "commands": self.commands,
             "candidate_source": str(self.candidate), "candidate_identity": initial_candidate, "checks": self.checks})

    def force_boundary(self, out):
        old = out / "force-old-layout"
        shutil.copytree(self.work / "baseline/numeric_shapes", old)
        before = snapshot(old)
        result = self.call(self.candidate, old, ["--mode", "matrix", "--species", "reference.nwk",
                           "--gene", "x.input", "--label", "x", "--force"])
        self.check(result.returncode != 0, "old-layout-force-refused")
        diagnostic = result.stderr.decode("utf-8", "replace")
        self.check(bool(re.search("layout|incompatible|legacy|compact|new label|new directory", diagnostic, re.I)),
                   "old-layout-explicit-diagnostic", diagnostic[-3000:])
        self.check(snapshot(old) == before, "old-layout-force-preserves-all-files-and-modes")
        # Current compact replacement: exact gene set, no stale g2 evidence.
        case = {"name": "force-current", "species": encoded("((A,B),(C,D));\n"),
                "panels": {"x": encoded("g2((A:1,B:2):3,(C:4,D:5):6);\ng1((A:7,B:8):9,(C:10,D:11):12);\n")},
                "view": "whole", "success": True}
        d = self.run_case(out / "force-current", self.candidate, case)
        (d / "replacement.input").write_text("g1((A:2,B:3):4,(C:5,D:6):7);\n")
        result = self.call(self.candidate, d, ["--mode", "matrix", "--species", "reference.nwk",
                           "--gene", "replacement.input", "--label", "x", "--force"])
        self.check(result.returncode == 0, "current-force-replacement-succeeds", result.stderr.decode("utf-8", "replace")[-3000:])
        self.check([r["gene_id"] for r in read_json(d / "x.compact.index.json")["genes"]] == ["g1"], "current-force-no-stale-index")
        self.check([json.loads(line)["gene_id"] for line in (d / "x.compact.jsonl").read_text().splitlines()] == ["g1"], "current-force-no-stale-ledger")
        (d / "invalid.input").write_text("g1((A:1,B:2):3,(C:4,Z:5):6);\n")
        before = snapshot(d)
        result = self.call(self.candidate, d, ["--mode", "matrix", "--species", "reference.nwk",
                           "--gene", "invalid.input", "--label", "x", "--force"])
        self.check(result.returncode != 0 and snapshot(d) == before, "invalid-force-preserves-current-complete-run")


def split_sides(key):
    if key.startswith("HX1:"):
        sides = key[4:].split("||")
        require(len(sides) == 2 and all(side.startswith("HS1:") for side in sides), "bad structured baseline split")
        return tuple({bytes.fromhex(token).decode("utf-8") for token in side[4:].split(",") if token} for side in sides)
    sides = key.split("||")
    require(len(sides) == 2, "bad baseline split")
    return tuple(set(side.split("..")) - {""} for side in sides)


def final_relative(name):
    return "inputs/" + name if name.endswith(".matrix_with_fuse.na_fuse.txt") else name


def set_bits(members, taxa):
    data = bytearray((len(taxa) + 7) // 8)
    for i, taxon in enumerate(taxa):
        if taxon in members:
            data[i // 8] |= 1 << (i % 8)
    require(set(members) <= set(taxa), "unknown taxon in baseline split")
    return data.hex()


def bit_key(left, right, taxa):
    minimum = min(left | right)
    return set_bits(right if minimum in left else left, taxa)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--frozen", type=Path, required=True)
    parser.add_argument("--candidate", type=Path)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--run-name", default="candidate-001")
    args = parser.parse_args()
    require(args.prepare_only or args.candidate, "--candidate is required after baseline preparation")
    require(re.fullmatch(r"[A-Za-z0-9_-]+", args.run_name), "run-name must be a simple basename")
    args.work_dir.mkdir(parents=True, exist_ok=True)
    h = None
    try:
        h = Harness(args)
        h.prepare() if args.prepare_only else h.execute()
    except Exception as error:
        if h:
            dump(h.work / ("prepare-failure.json" if args.prepare_only else args.run_name + "-failure.json"),
                 {"status": "FAIL", "error": str(error), "traceback": traceback.format_exc(),
                  "commands": h.commands, "checks": h.checks})
        raise
    print(f"PASS: {len(h.checks)} checks; {h.commands} bounded commands", flush=True)


if __name__ == "__main__":
    main()
