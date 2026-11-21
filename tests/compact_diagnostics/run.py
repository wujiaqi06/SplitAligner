#!/usr/bin/env python3
"""Exercise real-CLI gene diagnostics without changing frozen input handling.

Each invocation requires a new --work-dir and retains commands, inputs, stdout,
stderr and a machine-readable receipt, including failures. This is a small
diagnostic regression, not a performance or large-panel acceptance test.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess


REFERENCE = "((A,B),(C,D));\n"
GOOD = "g1((A:1,B:2):3,(C:4,D:5):6);\n"
CASES = (
    ("unknown_taxon", GOOD + "g2(A:1,B:2,Z:3);\n", 2, "g2", "Unknown taxon 'Z'"),
    ("parse_error", GOOD + "g2((A:1,B:2),C:3;\n", 2, "g2", "Expected"),
    ("duplicate_tip", GOOD + "g2(A:1,A:2,C:3);\n", 2, "g2", "Duplicate terminal taxon 'A'"),
    ("quoted_tip", GOOD + "g2('A':1,B:2,C:3);\n", 2, "g2", "Quoted labels"),
    ("nonfinite_length", GOOD + "g2(A:1e999,B:2,C:3);\n", 2, "g2", "not finite"),
    ("normalized_overflow", GOOD + "g2((A:1,B:2):1e308,(C:3,D:4):1e308);\n", 2, "g2", "Pooled unrooted branch length"),
    ("duplicate_gene", GOOD + GOOD, 2, "g1", "Duplicate gene ID"),
    ("missing_tree", GOOD + "g2 no tree\n", 2, None, "No '(' found"),
    ("blank_input", "\n \n", None, None, "No non-empty gene-tree records"),
    ("bom", "\ufeff" + GOOD, None, None, "BOM is not permitted"),
    ("bad_utf8", GOOD.encode() + b"g2(A:1,B:2,\xff:3);\n", 2, None, "Invalid UTF-8"),
    ("unicode_gene", GOOD + "基因二(A:1,B:2,Z:3);\n", 2, "基因二", "Unknown taxon 'Z'"),
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--work-dir", required=True, type=Path)
    args = parser.parse_args()
    candidate = args.candidate.resolve()
    root = args.work_dir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    records = []
    for name, data, line, gene, reason in CASES:
        case = root / name
        output = case / "output"
        inputs = case / "input space"
        inputs.mkdir(parents=True)
        output.mkdir()
        species = inputs / "reference.nwk"
        genes = inputs / "genes [original].txt"
        species.write_text(REFERENCE, encoding="utf-8")
        raw = data if isinstance(data, bytes) else data.encode("utf-8")
        genes.write_bytes(raw)
        # Use relative caller spelling: the controller changes into private
        # staging, but diagnostics must retain exactly this supplied label.
        user_path = "../input space/genes [original].txt"
        command = ["perl", str(candidate / "SplitAligner.pl"), "--mode", "matrix",
                   "--species", str(species), "--gene", user_path,
                   "--label", "diag", "--fusion-view", "whole"]
        proc = subprocess.run(command, cwd=output, capture_output=True,
                              env={**os.environ, "LC_ALL": "C", "LANG": "C"}, timeout=30)
        (case / "stdout.log").write_bytes(proc.stdout)
        (case / "stderr.log").write_bytes(proc.stderr)
        (case / "COMMAND.json").write_text(json.dumps(command, indent=2) + "\n")
        text = proc.stderr.decode("utf-8")
        context = user_path + (f" line {line}" if line else "") + (f" gene {gene}" if gene else "")
        checks = {
            "rejected": proc.returncode != 0,
            "original_source_context": f"[ERROR][{context}" in text,
            "single_error_marker": text.count("[ERROR]") == 1,
            "reason_preserved": reason in text,
            "original_input_unchanged": genes.read_bytes() == raw,
            "no_partial_publication_or_staging": not list(output.iterdir()),
        }
        records.append({"case": name, "exit_code": proc.returncode, "checks": checks,
                        "pass": all(checks.values()), "input_sha256": hashlib.sha256(raw).hexdigest()})
    receipt = {"candidate": str(candidate), "cases": records,
               "cases_passed": sum(x["pass"] for x in records), "cases_total": len(records),
               "checks_passed": sum(sum(x["checks"].values()) for x in records),
               "checks_total": sum(len(x["checks"]) for x in records)}
    receipt["status"] = "PASS" if receipt["cases_passed"] == len(records) else "FAIL"
    (root / "RESULTS.json").write_text(json.dumps(receipt, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(receipt, indent=2, ensure_ascii=False))
    return 0 if receipt["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
