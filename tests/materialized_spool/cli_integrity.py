#!/usr/bin/env python3
"""Bounded CLI injection into private spool/artifact boundaries; no source hooks."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import stat
import subprocess
import sys
import traceback

REFERENCE = b"(((A,B),C),(D,E));\n"
GENES = b"z((A:1.00,B:0):5,(D:6,E:7):8);\na((A:1,D:2):5,(B:6,E:7):8);\nm(((A:1,B:2):4,C):5,(D:-0,E:7):8);\n"
MODES = ["spool_byte", "spool_offset", "spool_identity", "ledger_during_read",
         "index_after", "dictionary_after", "ledger_after"]
HOOK = r'''
package SpoolIntegrityProbe;
use strict;
use warnings;
use JSON::PP;
use Cwd qw(getcwd);
use File::Spec;
use SplitAligner::Provenance qw(file_sha256);
our ($ACTIVE, $RECEIPT);

sub flip_file {
    my ($path) = @_;
    $RECEIPT->{mutated_path} = $path;
    $RECEIPT->{artifact_before} = file_sha256($path);
    open my $fh, '+<:raw', $path or die $!;
    read($fh, my $byte, 1) == 1 or die 'probe cannot read byte';
    seek($fh, 0, 0) or die $!;
    print {$fh} chr(ord($byte) ^ 1) or die $!;
    close $fh or die $!;
    $RECEIPT->{artifact_after} = file_sha256($path);
    ++$RECEIPT->{mutation_calls};
}

sub append_space {
    my ($path) = @_;
    $RECEIPT->{mutated_path} = $path;
    $RECEIPT->{artifact_before} = file_sha256($path);
    open my $fh, '>>:raw', $path or die $!;
    print {$fh} ' ' or die $!;
    close $fh or die $!;
    $RECEIPT->{artifact_after} = file_sha256($path);
    ++$RECEIPT->{mutation_calls};
}

sub import {
    return unless $0 =~ /SplitAligner\.pl\z/;
    $ACTIVE = 1;
    my $mode = $ENV{SA_SPOOL_PROBE_MODE};
    $RECEIPT = { mode => $mode, mutation_calls => 0, spool_read_calls => 0,
        materialize_calls => 0, public_materialize_calls => 0, json_decode_calls => 0 };
    require SplitAligner::CompactPipeline;
    require SplitAligner::CompactKernel;
    my $real_run = \&SplitAligner::CompactPipeline::run_compact_matrix;
    my $real_read = \&SplitAligner::CompactPipeline::read_spool_block;
    my $real_private = \&SplitAligner::CompactKernel::_materialize;
    my $real_public = \&SplitAligner::CompactKernel::materialize_gene;
    my $real_decode = \&JSON::PP::decode;
    no warnings 'redefine';
    *SplitAligner::CompactKernel::_materialize = sub {
        ++$RECEIPT->{materialize_calls};
        return $real_private->(@_);
    };
    my $public = sub {
        ++$RECEIPT->{public_materialize_calls};
        return $real_public->(@_);
    };
    *SplitAligner::CompactKernel::materialize_gene = $public;
    *SplitAligner::CompactPipeline::materialize_gene = $public
        if SplitAligner::CompactPipeline->can('materialize_gene');
    *JSON::PP::decode = sub {
        ++$RECEIPT->{json_decode_calls};
        die '[Probe] JSON decode called in fresh matrix path' if $mode eq 'no_decode';
        return $real_decode->(@_);
    };
    *SplitAligner::CompactPipeline::read_spool_block = sub {
        my ($fh, $entry, $ctx, $limits) = @_;
        ++$RECEIPT->{spool_read_calls};
        $RECEIPT->{spool_path} = File::Spec->rel2abs($fh->filename, getcwd());
        if ($RECEIPT->{spool_read_calls} == 1) {
            if ($mode eq 'spool_byte') {
                seek($fh, $entry->{offset}, 0) or die $!;
                read($fh, my $byte, 1) == 1 or die 'probe cannot read spool byte';
                seek($fh, $entry->{offset}, 0) or die $!;
                print {$fh} chr(ord($byte) ^ 1) or die $!;
                ++$RECEIPT->{mutation_calls};
            } elsif ($mode eq 'spool_offset' || $mode eq 'spool_identity') {
                my %wrong = %{$entry};
                if ($mode eq 'spool_offset') { ++$wrong{offset}; }
                else { ++$wrong{input_record_id}; }
                ++$RECEIPT->{mutation_calls};
                return $real_read->($fh, \%wrong, $ctx, $limits);
            } elsif ($mode eq 'ledger_during_read') {
                flip_file($RECEIPT->{ledger_path});
            }
        }
        return $real_read->(@_);
    };
    *SplitAligner::CompactPipeline::run_compact_matrix = sub {
        my %args = @_;
        $RECEIPT->{stage_directory} = getcwd();
        $RECEIPT->{snapshot_gene_path} = $args{gene_path};
        $RECEIPT->{ledger_path} = "$args{label}.compact.jsonl";
        my $result = $real_run->(@_);
        $RECEIPT->{run_returned} = 1;
        if ($mode =~ /\A(index|dictionary|ledger)_after\z/) {
            append_space($result->{$1});
        }
        return $result;
    };
}

END {
    if ($ACTIVE) {
        my $status = $?;
        my $path = $ENV{SA_SPOOL_PROBE_RECEIPT};
        open my $out, '>:raw', $path or die "Cannot write probe receipt: $!";
        print {$out} JSON::PP->new->canonical(1)->pretty(1)->utf8(1)->encode($RECEIPT) or die $!;
        close $out or die $!;
        $? = $status;
    }
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


def inventory(root):
    result = {}
    for path in sorted(root.rglob("*")):
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
        result[str(path.relative_to(root))] = item
    return result


def source_binding(source):
    paths = [source / "SplitAligner.pl"] + sorted((source / "lib/SplitAligner").glob("*.pm"))
    paths += sorted((source / "scripts").glob("*.pl"))
    return {str(p.relative_to(source)): sha(p.read_bytes()) for p in paths}


def initialize_case(root):
    root.mkdir()
    for name in ["input", "public", "tmp", "evidence", "hook"]:
        (root / name).mkdir()
    (root / "input/reference.nwk").write_bytes(REFERENCE)
    (root / "input/genes.nwk").write_bytes(GENES)
    (root / "hook/SpoolIntegrityProbe.pm").write_text(HOOK, encoding="utf-8")


def run_cli(source, root, mode=None, force=False, name="run"):
    command = ["/usr/bin/perl", str(source / "SplitAligner.pl"), "--mode", "matrix",
               "--species", str(root / "input/reference.nwk"), "--gene", str(root / "input/genes.nwk"),
               "--label", "panel", "--fusion-view", "whole"]
    if force:
        command.append("--force")
    env = dict(os.environ)
    env.pop("PERL5OPT", None)
    env.pop("PERL5LIB", None)
    env.update(LC_ALL="C", TMPDIR=str(root / "tmp"), PERL_HASH_SEED="4179", PERL_PERTURB_KEYS="1")
    if mode:
        env.update(PERL5OPT="-MSpoolIntegrityProbe", PERL5LIB=os.pathsep.join([str(root / "hook"), str(source / "lib")]),
                   SA_SPOOL_PROBE_MODE=mode, SA_SPOOL_PROBE_RECEIPT=str(root / "evidence/hook.json"))
    completed = subprocess.run(command, cwd=root / "public", env=env, capture_output=True, timeout=120)
    (root / "evidence" / f"{name}.stdout.log").write_bytes(completed.stdout)
    (root / "evidence" / f"{name}.stderr.log").write_bytes(completed.stderr)
    dump(root / "evidence" / f"{name}.command.json", {"command": command, "cwd": str(root / "public"),
         "exit_code": completed.returncode, "environment": {k: v for k, v in env.items()
         if k in ["LC_ALL", "TMPDIR", "PERL_HASH_SEED", "PERL_PERTURB_KEYS", "PERL5OPT", "PERL5LIB"] or k.startswith("SA_SPOOL_")},
         "stdout_sha256": sha(completed.stdout), "stderr_sha256": sha(completed.stderr)})
    return completed


def assert_clean(root, hook):
    require(not list((root / "tmp").iterdir()), "temporary input snapshot residue")
    require(not Path(hook["stage_directory"]).exists(), "private stage residue")
    require(not Path(hook["snapshot_gene_path"]).exists(), "gene snapshot residue")
    require(not Path(hook["spool_path"]).exists(), "materialized spool residue")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--work-dir", type=Path, required=True)
    args = parser.parse_args()
    source, work = args.source.resolve(), args.work_dir.resolve()
    work.mkdir(parents=True, exist_ok=False)
    receipt = {"status": "RUNNING", "source": str(source), "started_utc": datetime.now(timezone.utc).isoformat(),
               "source_sha256_before": source_binding(source), "cases": []}
    dump(work / "RESULTS.json", receipt)
    try:
        root = work / "fresh-no-json-decode"
        initialize_case(root)
        before = inventory(root / "input")
        result = run_cli(source, root, "no_decode")
        require(result.returncode == 0, "fresh no-decode run failed: " + result.stderr.decode(errors="replace"))
        hook = json.loads((root / "evidence/hook.json").read_text())
        require(hook["json_decode_calls"] == 0, "normal fresh path used JSON decode")
        require(hook["materialize_calls"] == 3, "expected exactly one materialization per gene")
        require(hook["public_materialize_calls"] == 0, "second materialize_gene pass occurred")
        require(hook["spool_read_calls"] == 3 and hook["run_returned"] == 1, "spool path did not finish")
        require(inventory(root / "input") == before, "positive case changed inputs")
        assert_clean(root, hook)
        receipt["cases"].append({"name": root.name, "status": "PASS", "json_decodes": 0,
                                "materializations": 3, "genes": 3, "spool_reads": 3})
        dump(work / "RESULTS.json", receipt)
        for mode in MODES:
            for force in [False, True]:
                root = work / f"{mode}-{'force-existing' if force else 'fresh'}"
                initialize_case(root)
                if force:
                    setup = run_cli(source, root, name="setup")
                    require(setup.returncode == 0, f"{root.name} setup failed")
                public_before = inventory(root / "public")
                input_before = inventory(root / "input")
                dump(root / "evidence/public_before.json", public_before)
                result = run_cli(source, root, mode, force)
                error = result.stderr.decode(errors="replace")
                require(result.returncode != 0, f"{root.name} unexpectedly succeeded")
                require("completed successfully" not in error, f"{root.name} claimed completion")
                hook = json.loads((root / "evidence/hook.json").read_text())
                require(hook["mutation_calls"] == 1, f"{root.name} injection did not run exactly once")
                require(hook["materialize_calls"] == 3 and hook["public_materialize_calls"] == 0,
                        f"{root.name} materialization counts changed")
                if mode.startswith("spool_"):
                    require("[Materialized spool]" in error, f"{root.name} wrong rejection: {error}")
                    require(hook["spool_read_calls"] == 1 and not hook.get("run_returned"), "spool fault was not rejected at first read")
                elif mode == "ledger_during_read":
                    require("Written compact ledger differs from its accumulated fingerprint" in error,
                            f"{root.name} wrong rejection: {error}")
                    require(hook["spool_read_calls"] == 3 and not hook.get("run_returned"), "ledger guard did not follow completed spool reads")
                else:
                    require("Compact artifact differs from its write-time seal" in error,
                            f"{root.name} wrong rejection: {error}")
                    require(hook.get("run_returned") == 1, "after-stage injection did not follow successful computation")
                if "artifact_before" in hook:
                    require(hook["artifact_before"] != hook["artifact_after"], "artifact mutation did not change bytes")
                require(inventory(root / "public") == public_before, f"{root.name} changed published data or inodes")
                require(inventory(root / "input") == input_before, f"{root.name} changed original inputs")
                assert_clean(root, hook)
                dump(root / "evidence/public_after.json", inventory(root / "public"))
                receipt["cases"].append({"name": root.name, "status": "PASS", "exit_code": result.returncode,
                    "original_inputs_unchanged": True, "published_tree_unchanged": True, "temporary_files_cleaned": True})
                dump(work / "RESULTS.json", receipt)
        receipt["source_sha256_after"] = source_binding(source)
        require(receipt["source_sha256_after"] == receipt["source_sha256_before"], "source changed during suite")
        receipt["status"] = "PASS"
        receipt["scope"] = "15 bounded full-CLI cases; no production hooks or source edits; force ownership may legitimately decode existing manifests"
        receipt["finished_utc"] = datetime.now(timezone.utc).isoformat()
        dump(work / "RESULTS.json", receipt)
        print(f"PASS {len(receipt['cases'])} CLI integrity cases: {work}")
        return 0
    except Exception as error:
        receipt.update(status="FAIL", error=str(error), traceback=traceback.format_exc(), source_sha256_after=source_binding(source))
        dump(work / "RESULTS.json", receipt)
        print(receipt["traceback"], file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
