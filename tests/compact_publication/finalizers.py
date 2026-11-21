#!/usr/bin/env python3
"""Run old/new finalizer CLIs on independent copies of new compact matrices.

Only a three-gene toy panel is generated, using the candidate. No old matrix run
or benchmark is performed. Finalizer science is compared byte-for-byte and each
result's ownership is checked against its files and by read-only force preflight.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import traceback


ROOT = Path(__file__).resolve().parents[2]
REF = "((A,(B,C)),(D,E));\n"
FIXED = {
    "z_full": "((A:1,(B:2,C:3):4):5,(D:6,E:7):8);",
    "a_drop": "((A:1,B:6):5,(D:6,E:7):8);",
    "m_topology": "((A:2,(B:3,C:4):5):6,(D:7,E:8):9);",
}
FREE = {
    "z_full": "((A:1,(B:2,C:3):4):5,(D:6,E:7):8);",
    "a_drop": "((A:1,D:6):5,(B:2,E:7):8);",
    "m_topology": "((A:2,(B:3,D:4):5):6,(C:7,E:8):9);",
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def dump(path, data):
    path.write_text(json.dumps(data, ensure_ascii=False, sort_keys=True, indent=2) + "\n", encoding="utf-8")


def read_json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def require(ok, message):
    if not ok:
        raise AssertionError(message)


def binding(root):
    paths = [root / "SplitAligner.pl", *sorted((root / "lib/SplitAligner").glob("*.pm")),
             *sorted((root / "scripts").glob("*.pl"))]
    return {str(p.relative_to(root)): sha(p) for p in paths}


def content_snapshot(root):
    return {p.relative_to(root).as_posix(): sha(p) for p in sorted(root.rglob("*")) if p.is_file()}


VERIFY_OWNER = r'''
use strict; use warnings;
use File::Spec; use File::Basename qw(basename);
use SplitAligner::Provenance qw(read_json_file);
use SplitAligner::OutputOwnership qw(preflight_output_ownership);
my ($base, $path) = @ARGV;
my $meta = read_json_file($path);
my $owner = $meta->{ownership};
my @items = map { +{
    source => basename($_->{relative_path}),
    destination => File::Spec->catfile($base, $_->{relative_path}),
    role => $_->{role}, type => $_->{type}, shared => $_->{shared},
} } @{$owner->{outputs}};
push @items, {source => basename($path), destination => $path,
    role => $owner->{manifest_role}, type => 'file', owner_manifest => 1};
my $auth = preflight_output_ownership(mode => $owner->{owner_mode},
    owner_label => $owner->{owner_label}, owner_manifest_path => $path,
    items => \@items, force => 1);
die "Prior owner was not verified\n" unless $auth->{has_prior_owner};
print "PASS read-only ownership preflight\n";
'''


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--candidate", type=Path, default=ROOT)
    ap.add_argument("--frozen", type=Path, required=True)
    ap.add_argument("--work-dir", type=Path, required=True)
    args = ap.parse_args()
    candidate, frozen, work = (p.resolve() for p in (args.candidate, args.frozen, args.work_dir))
    require(not work.exists() or not any(work.iterdir()), "Evidence directory must be new or empty")
    work.mkdir(parents=True, exist_ok=True)
    logs = work / "logs"
    logs.mkdir()
    env = dict(os.environ, LC_ALL="C", PERL_HASH_SEED="741", PERL_PERTURB_KEYS="2")
    initial = {"candidate": binding(candidate), "frozen": binding(frozen), "test_driver": sha(Path(__file__))}
    dump(work / "source_binding.json", initial)
    receipt = {"status": "RUNNING", "scope": "Three-gene real finalize/finalize_fix CLI comparison; no benchmark",
               "candidate": str(candidate), "frozen": str(frozen), "generation": [], "comparisons": []}

    def run(name, argv, cwd, source):
        directory = logs / name
        directory.mkdir()
        before = binding(source)
        result = subprocess.run(argv, cwd=cwd, env=env, capture_output=True, text=True,
                                encoding="utf-8", timeout=60)
        after = binding(source)
        dump(directory / "command.json", {"argv": [str(a) for a in argv], "cwd": str(cwd),
                                           "exit_code": result.returncode, "LC_ALL": "C"})
        dump(directory / "source_binding.json", {"before": before, "after": after, "unchanged": before == after})
        (directory / "stdout.txt").write_text(result.stdout, encoding="utf-8")
        (directory / "stderr.txt").write_text(result.stderr, encoding="utf-8")
        require(before == after, name + ": source changed during command")
        require(result.returncode == 0, name + ": " + result.stderr)
        return result

    try:
        seed = work / "seed"
        inputs = seed / "inputs"
        inputs.mkdir(parents=True)
        (inputs / "reference.nwk").write_text(REF, encoding="utf-8")
        for label, genes in [("free", FREE), ("fix", FIXED)]:
            (inputs / (label + ".genes")).write_text("".join(g + t + "\n" for g, t in genes.items()), encoding="utf-8")
            target = seed / "matrices" / label
            target.mkdir(parents=True)
            run("generate_" + label, ["perl", str(candidate / "SplitAligner.pl"), "--mode", "matrix",
                "--species", str(inputs / "reference.nwk"), "--gene", str(inputs / (label + ".genes")),
                "--label", label, "--fusion-view", "whole"], target, candidate)
            meta = read_json(target / (label + ".run_manifest.json"))
            require(meta["coordinate_state"]["gene_record_count"] == 3, "Wrong fixture gene count")
            receipt["generation"].append({"label": label, "genes": 3,
                "manifest_sha256": sha(target / (label + ".run_manifest.json"))})
        seed_before = content_snapshot(seed)
        results = {}
        for mode in ("finalize", "finalize_fix"):
            for name, source in [("frozen", frozen), ("candidate", candidate)]:
                case = work / "cases" / mode / name
                shutil.copytree(seed, case)
                cwd = case / "output"
                cwd.mkdir()
                matrices = case / "matrices"
                argv = ["perl", str(source / "SplitAligner.pl"), "--mode", mode,
                        "--fix", str(matrices / "fix/fix.matrix_with_fuse.txt"), "--final_label", "result"]
                if mode == "finalize":
                    argv += ["--free", str(matrices / "free/free.matrix_with_fuse.txt"),
                             "--species_tree", str(matrices / "free/species_tree.forSplit.nwk")]
                original = content_snapshot(case)
                run(mode + "_" + name, argv, cwd, source)
                after = content_snapshot(case)
                require(all(after.get(path) == value for path, value in original.items()),
                        mode + "/" + name + ": original matrix/input artifacts were modified")
                manifest_path = cwd / "result.finalize_manifest.json"
                meta = read_json(manifest_path)
                science = {}
                for entry in meta["ownership"]["outputs"]:
                    require(entry["type"] == "file", "Unexpected directory in finalizer fixture output")
                    relative = entry["relative_path"]
                    path = (cwd / relative).resolve()
                    require(path.is_relative_to(case), "Finalizer output escaped its independent package")
                    info = path.lstat()
                    fingerprint = entry["fingerprint"]
                    require(sha(path) == fingerprint["sha256"], "Ownership hash differs: " + str(path))
                    require(info.st_size == fingerprint["size"] and info.st_nlink == fingerprint["nlink"],
                            "Ownership size/link count differs: " + str(path))
                    require(format(info.st_mode & 0o7777, "04o") == fingerprint["mode"],
                            "Ownership file mode differs: " + str(path))
                    science[relative] = path.read_bytes()
                before_verify = content_snapshot(case)
                run(mode + "_" + name + "_ownership", ["perl", "-I" + str(source / "lib"), "-e", VERIFY_OWNER,
                    str(cwd), str(manifest_path)], cwd, source)
                require(content_snapshot(case) == before_verify, "Ownership verification modified results")
                results[(mode, name)] = (science, meta)
            old, old_meta = results[(mode, "frozen")]
            new, new_meta = results[(mode, "candidate")]
            require(old.keys() == new.keys(), mode + ": scientific output sets differ")
            for relative in sorted(old):
                require(old[relative] == new[relative], mode + ": scientific bytes differ: " + relative)
            for field in ("outputs", "primitive_axis_sha256", "primitive_fusion_source", "free_projection_accounting", "inputs"):
                require(old_meta.get(field) == new_meta.get(field), mode + ": scientific manifest field differs: " + field)
            receipt["comparisons"].append({"mode": mode, "status": "PASS", "ownership_verified": ["frozen", "candidate"],
                "scientific_files": {relative: hashlib.sha256(new[relative]).hexdigest() for relative in sorted(new)},
                "note": "Manifest program identity is intentionally different; science and input bindings compared."})
            print("PASS", mode, len(new), "scientific files; both ownership records verified", flush=True)
        require(content_snapshot(seed) == seed_before, "Shared seed was changed")
        require(binding(candidate) == initial["candidate"] and binding(frozen) == initial["frozen"], "Sources changed across suite")
        receipt["status"] = "PASS"
        receipt["source_unchanged_across_suite"] = True
    except Exception as exc:
        receipt["status"] = "FAIL"
        receipt["failure"] = str(exc)
        (work / "failure.txt").write_text(traceback.format_exc(), encoding="utf-8")
        print(traceback.format_exc(), file=sys.stderr)
    dump(work / "RESULTS.json", receipt)
    print(json.dumps({"status": receipt["status"], "comparisons": len(receipt["comparisons"]), "evidence": str(work)}))
    return 0 if receipt["status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
