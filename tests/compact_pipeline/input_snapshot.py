#!/usr/bin/env python3
"""Bounded CLI snapshot/use/publication tests with process-local Perl injection.

No source files are changed. A generated PERL5OPT preload wraps the compact
stage after the controller froze both inputs. Evidence belongs to --work-dir;
existing evidence is never overwritten. This is not a performance benchmark.
"""
import argparse
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path
import stat
import subprocess
import sys
import traceback
from zoneinfo import ZoneInfo


REFERENCE = b"((A,B),((C,D),(E,F)));\n"
GENES = b"g1(A:1,C:2,E:3);\ng2((A:1.00,B:2):3,(C:4,D:5,E:6,F:7):8);\n"
REPLACEMENTS = {
    "gene": b"g1(A:101,C:202,E:303);\ng2((A:9,B:8):7,(C:6,D:5,E:4,F:3):2);\n",
    "species": b"((A,C),((B,D),(E,F)));\n",
}
CAPTURE = [
    "species_tree.forSplit.nwk", "species_tree.FigTree.tre", "species_tree.splits.txt",
    "species_tree.branch_map.txt", "species_tree.primitive_axis.tsv",
    "panel.compact.jsonl", "panel.compact.index.json", "panel.taxon_dictionary.json",
    "panel.gene_id_map.tsv", "panel.primitive_state.tsv", "panel.matrix_no_fuse.txt",
    "panel.matrix_with_fuse.txt", "panel.fusion_cells.tsv",
    "panel_split_branch_label/panel.primitive_truth.tsv",
    "panel_split_branch_label/panel.projected_groups.tsv",
    "panel_split_branch_label/panel.branch_patterns.txt",
]
HOOK = r'''
package CompactSnapshotProbe;
use strict;
use warnings;
use File::Spec;
use File::Basename qw(dirname);
use File::Path qw(make_path);
use File::Copy qw(copy);
use Cwd qw(getcwd);
use JSON::PP;
use SplitAligner::Provenance qw(file_sha256);

sub replace_original {
    my ($receipt) = @_;
    my $source = $ENV{SA_COMPACT_MUTATE_PATH};
    $receipt->{original_before_mutation} = file_sha256($source);
    open my $in, '<:raw', $ENV{SA_COMPACT_REPLACEMENT} or die $!;
    local $/; my $new = <$in>; close $in or die $!;
    open my $out, '>:raw', $source or die $!;
    print {$out} $new or die $!; close $out or die $!;
    $receipt->{original_after_mutation} = file_sha256($source);
    $receipt->{mutation_calls}++;
}

sub import {
    return unless $0 =~ /SplitAligner\.pl\z/;
    require SplitAligner::CompactPipeline;
    my $real = \&SplitAligner::CompactPipeline::run_compact_matrix;
    no warnings 'redefine';
    *SplitAligner::CompactPipeline::run_compact_matrix = sub {
        my %args = @_;
        my $receipt = {
            phase => $ENV{SA_COMPACT_PHASE}, hook_calls => 1,
            stage_directory => getcwd(), snapshot_gene_path => $args{gene_path},
            snapshot_gene_sha256_before => file_sha256($args{gene_path}),
            input_sha256_arguments => $args{input_sha256},
        };
        replace_original($receipt) if $ENV{SA_COMPACT_PHASE} eq 'before';
        my $result = $real->(@_);
        $receipt->{stage_returned} = 1;
        $receipt->{snapshot_gene_sha256_after} = file_sha256($args{gene_path});
        replace_original($receipt) if $ENV{SA_COMPACT_PHASE} eq 'after';
        open my $list, '<:raw', $ENV{SA_COMPACT_CAPTURE_LIST} or die $!;
        while (my $relative = <$list>) {
            chomp $relative;
            my $dest = File::Spec->catfile($ENV{SA_COMPACT_EVIDENCE}, 'stage', $relative);
            make_path(dirname($dest));
            copy($relative, $dest) or die "Cannot capture $relative: $!";
            $receipt->{captured_sha256}{$relative} = file_sha256($dest);
        }
        close $list or die $!;
        my $path = File::Spec->catfile($ENV{SA_COMPACT_EVIDENCE}, 'hook.json');
        open my $out, '>:raw', $path or die $!;
        print {$out} JSON::PP->new->canonical(1)->pretty(1)->utf8(1)->encode($receipt) or die $!;
        close $out or die $!;
        return $result;
    };
}
1;
'''


def sha(data):
    return hashlib.sha256(data).hexdigest()


def dump(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def inventory(directory):
    result = {}
    for path in sorted(directory.rglob("*")):
        st = path.lstat()
        item = {"mode": stat.S_IMODE(st.st_mode), "inode": st.st_ino, "device": st.st_dev}
        if path.is_symlink():
            item.update(type="symlink", target=os.readlink(path))
        elif path.is_file():
            item.update(type="file", bytes=st.st_size, sha256=sha(path.read_bytes()))
        elif path.is_dir():
            item.update(type="directory")
        else:
            item.update(type="other")
        result[str(path.relative_to(directory))] = item
    return result


def source_binding(source):
    files = [source / "SplitAligner.pl"] + sorted((source / "lib/SplitAligner").glob("*.pm"))
    files += [source / "scripts/label_species_tree.pl", source / "scripts/tree_to_splits.pl"]
    return {str(p.relative_to(source)): sha(p.read_bytes()) for p in files}


def initialize_case(root):
    root.mkdir()
    for name in ("input", "public", "tmp", "evidence"):
        (root / name).mkdir()
    (root / "input/species.nwk").write_bytes(REFERENCE)
    (root / "input/genes.nwk").write_bytes(GENES)


def run_cli(source, root, *, injected=False, phase=None, target=None, force=False, name="run"):
    command = ["/usr/bin/perl", str(source / "SplitAligner.pl"), "--mode", "matrix",
               "--species", str(root / "input/species.nwk"), "--gene", str(root / "input/genes.nwk"),
               "--label", "panel", "--fusion-view", "whole"]
    if force:
        command.append("--force")
    env = dict(os.environ)
    # Avoid inheriting another test's injection; source lookup remains explicit.
    env.pop("PERL5OPT", None)
    env.pop("PERL5LIB", None)
    env.update(LC_ALL="C", TMPDIR=str(root / "tmp"), PERL_HASH_SEED="7351", PERL_PERTURB_KEYS="1")
    if injected:
        hook = root / "hook"
        hook.mkdir()
        (hook / "CompactSnapshotProbe.pm").write_text(HOOK, encoding="utf-8")
        (root / "replacement").write_bytes(REPLACEMENTS[target])
        (root / "capture-list").write_text("\n".join(CAPTURE) + "\n")
        target_path = root / "input" / ("genes.nwk" if target == "gene" else "species.nwk")
        env.update(PERL5OPT="-MCompactSnapshotProbe", PERL5LIB=os.pathsep.join([str(hook), str(source / "lib")]),
                   SA_COMPACT_MUTATE_PATH=str(target_path), SA_COMPACT_REPLACEMENT=str(root / "replacement"),
                   SA_COMPACT_PHASE=phase, SA_COMPACT_EVIDENCE=str(root / "evidence"),
                   SA_COMPACT_CAPTURE_LIST=str(root / "capture-list"))
    completed = subprocess.run(command, cwd=root / "public", env=env, capture_output=True, timeout=120)
    (root / "evidence" / f"{name}.stdout.txt").write_bytes(completed.stdout)
    (root / "evidence" / f"{name}.stderr.txt").write_bytes(completed.stderr)
    record = {"command": command, "cwd": str(root / "public"), "exit_code": completed.returncode,
              "environment": {k: v for k, v in env.items() if k in ("LC_ALL", "TMPDIR", "PERL_HASH_SEED", "PERL_PERTURB_KEYS", "PERL5OPT", "PERL5LIB") or k.startswith("SA_COMPACT_")},
              "stdout_sha256": sha(completed.stdout), "stderr_sha256": sha(completed.stderr)}
    dump(root / "evidence" / f"{name}.command.json", record)
    return completed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--work-dir", type=Path, required=True)
    args = parser.parse_args()
    source = args.source.resolve()
    work = args.work_dir.resolve()
    work.mkdir(parents=True, exist_ok=False)
    receipt = {"schema": "SplitAligner-compact-input-snapshot-tests-v1", "status": "RUNNING",
               "date": datetime.now(ZoneInfo("Asia/Tokyo")).isoformat(), "source": str(source),
               "source_sha256_before": source_binding(source), "cases": []}
    dump(work / "RECEIPT.json", receipt)
    try:
        baseline = work / "unchanged-control"
        initialize_case(baseline)
        control = run_cli(source, baseline)
        require(control.returncode == 0, "unchanged-input control failed")
        control_hashes = {name: sha((baseline / "public" / name).read_bytes()) for name in CAPTURE}
        dump(work / "CONTROL_STAGE_SHA256.json", control_hashes)
        for phase in ("before", "after"):
            for target in ("gene", "species"):
                for overwrite in (False, True):
                    name = f"{phase}-{target}-{'force-existing' if overwrite else 'new-output'}"
                    root = work / name
                    initialize_case(root)
                    if overwrite:
                        setup = run_cli(source, root, name="setup")
                        require(setup.returncode == 0, f"{name}: old candidate setup failed")
                    before = inventory(root / "public")
                    dump(root / "evidence/public-before.json", before)
                    completed = run_cli(source, root, injected=True, phase=phase, target=target, force=overwrite)
                    text = (completed.stdout + completed.stderr).decode("utf-8", errors="replace")
                    require(completed.returncode != 0, f"{name}: mutation must reject publication")
                    require("[Input immutability]" in text, f"{name}: failure is not the intended input guard: {text}")
                    require("completed successfully" not in text, f"{name}: incorrect completion message")
                    hook = json.loads((root / "evidence/hook.json").read_text())
                    require(hook["hook_calls"] == hook["mutation_calls"] == 1, f"{name}: injection did not run exactly once")
                    require(hook["stage_returned"] == 1, f"{name}: compact computation did not finish")
                    require(hook["snapshot_gene_sha256_before"] == hook["snapshot_gene_sha256_after"] == sha(GENES),
                            f"{name}: computation snapshot changed")
                    require(hook["input_sha256_arguments"] == {"gene_trees": sha(GENES), "species_tree": sha(REFERENCE)},
                            f"{name}: provenance arguments no longer identify frozen inputs")
                    require(Path(hook["snapshot_gene_path"]) != root / "input/genes.nwk", f"{name}: original path passed as snapshot")
                    require(hook["captured_sha256"] == control_hashes, f"{name}: intermediate outputs differ from frozen-input control")
                    original = GENES if target == "gene" else REFERENCE
                    require(hook["original_before_mutation"] == sha(original) and
                            hook["original_after_mutation"] == sha(REPLACEMENTS[target]), f"{name}: mutation evidence inconsistent")
                    changed_path = root / "input" / ("genes.nwk" if target == "gene" else "species.nwk")
                    require(changed_path.read_bytes() == REPLACEMENTS[target], f"{name}: test mutation was unexpectedly undone")
                    after = inventory(root / "public")
                    dump(root / "evidence/public-after.json", after)
                    require(after == before, f"{name}: published outputs changed despite failure")
                    require(not list((root / "tmp").iterdir()), f"{name}: temporary input snapshots leaked")
                    receipt["cases"].append({"name": name, "status": "PASS", "exit_code": completed.returncode,
                                             "intermediate_files_compared": len(control_hashes), "public_tree_unchanged": True,
                                             "source_mutation_persisted": True, "input_snapshots_cleaned": True})
                    dump(work / "RECEIPT.json", receipt)
        receipt["source_sha256_after"] = source_binding(source)
        require(receipt["source_sha256_after"] == receipt["source_sha256_before"], "source files changed during checks")
        receipt["status"] = "PASS"
        receipt["scope"] = "8 bounded CLI injections; 16 intermediate artifacts each; new publication rejected and old candidate outputs unchanged"
        dump(work / "RECEIPT.json", receipt)
        print(f"PASS: {len(receipt['cases'])} snapshot/publication cases; evidence: {work}")
        return 0
    except Exception as error:
        receipt["status"] = "FAIL"
        receipt["error"] = str(error)
        receipt["traceback"] = traceback.format_exc()
        receipt["source_sha256_after"] = source_binding(source)
        dump(work / "RECEIPT.json", receipt)
        print(receipt["traceback"], file=sys.stderr)
        print(f"FAIL evidence retained: {work}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
