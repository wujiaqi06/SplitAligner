"""Count free-tree outcomes with an explicit unit; lengths do not enter the sum."""
import argparse
import csv
import json
from itertools import groupby
from pathlib import Path
import sys


def counts(directory, unit):
    ctx = json.loads((directory / "context.json").read_text(encoding="utf-8"))
    if ctx["role"] != "free" or ctx["scope"] not in ("full", "support"):
        raise ValueError("use the output of a free full or support run")
    manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
    if (manifest["status"] != "COMPLETED" or manifest["role"] != ctx["role"] or
            manifest["scope"] != ctx["scope"] or manifest["axis_count"] != len(ctx["axes"])):
        raise ValueError("expected a completed run with matching scope and primitive axis")
    with (directory / "mapping.jsonl").open(encoding="utf-8") as mappings, \
         (directory / "composites.jsonl").open(encoding="utf-8") as composites:
        groups = groupby(map(json.loads, composites), key=lambda r: r["gene_idx"])
        batch = next(groups, None)
        n_genes = 0
        for idx, line in enumerate(mappings):
            n_genes += 1
            record = json.loads(line)
            states = record["raw_states"]
            if (idx != record["gene_idx"] or len(states) != len(ctx["axes"]) or
                    set(states) - set("DFSU")):
                raise ValueError("invalid mapping row")
            members = {i for i, state in enumerate(states) if state == "F"}
            expected_members = len(members)
            fuse = joint = n_groups = 0
            if batch is not None and batch[0] == idx:
                for group in batch[1]:
                    axes = set(group["members"])
                    if (group["gene_id"] != record["gene_id"] or
                            group["coverage_id"] != record["coverage_id"] or
                            len(axes) < 2 or len(axes) != len(group["members"]) or
                            not axes <= members):
                        raise ValueError("inconsistent group membership")
                    members -= axes
                    n_groups += 1
                    weight = 1 if unit == "coordinate" else len(axes)
                    if group["group_state"] == "Mapped" and group["recovered"] is True:
                        fuse += weight
                    elif group["group_state"] == "NA_fuse_topo" and group["recovered"] is False:
                        joint += weight
                    else:
                        raise ValueError("inconsistent group recovery")
                batch = next(groups, None)
            if members:
                raise ValueError("missing group for a fused member")
            total = len(states) if unit == "primitive" else len(states) - expected_members + n_groups
            terms = (states.count("D"), states.count("S"), fuse, states.count("U"), joint)
            if total != sum(terms):
                raise ValueError("topology accounting mismatch")
            yield (record["gene_id"], unit, total, *terms)
        if batch is not None or n_genes != manifest["gene_count"]:
            raise ValueError("group or record count differs from mapping run")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="completed Fast free full/support directory")
    parser.add_argument("--unit", choices=("coordinate", "primitive"), required=True)
    args = parser.parse_args()
    writer = csv.writer(sys.stdout, delimiter="\t", lineterminator="\n")
    writer.writerow(("gene_id", "unit", "Total", "Mapped", "NA_struct", "NA_fuse", "NA_topo", "NA_fuse_topo"))
    writer.writerows(counts(args.output, args.unit))
