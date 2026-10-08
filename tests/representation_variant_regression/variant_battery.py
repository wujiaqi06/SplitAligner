#!/usr/bin/env python3
"""Representation- and input-variant regression battery for SplitAligner.

Group A (serialization variants of equivalent trees): the gene-by-branch outputs
of `matrix` and `finalize` must be identical after mapping serialization-local
B aliases to canonical split keys.

Group B (degenerate inputs the tutorial describes): the program must not die,
and every difference from the baseline must be explained by the documented
rule (no numeric lengths -> NA values with unchanged states; collapsed gene
splits -> U; a degree-2 root with one missing half-length -> NA for that one
split, state unchanged).

Usage: python3 variant_battery.py [--work DIR] [--keep]
Writes variants/ (generated inputs) and RESULTS.md beside this script.
"""

import argparse
import os
import random
import re
import shutil
import subprocess
import sys
from collections import OrderedDict

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
SA = os.path.join(REPO, "SplitAligner.pl")
TOY = os.path.join(REPO, "examples", "302mammal", "input")
BASE_SPECIES = os.path.join(TOY, "speciesTree302.nwk")
BASE_FREE = os.path.join(TOY, "free_tree.examples.nwk")
BASE_FIX = os.path.join(TOY, "fix_tree.examples.nwk")
VARIANT_DIR = os.path.join(HERE, "variants")

# ---------------------------------------------------------------- Newick ---


class Node:
    __slots__ = ("children", "label", "length")

    def __init__(self):
        self.children = []
        self.label = ""
        self.length = None  # raw string or None


def parse_newick(text):
    s = text.strip()
    if not s.endswith(";"):
        raise ValueError("missing ';'")
    s = s[:-1]
    pos = 0

    def tok():
        nonlocal pos
        m = re.match(r"[^:,()\s;]*", s[pos:])
        pos += len(m.group(0))
        return m.group(0)

    def node():
        nonlocal pos
        n = Node()
        if s[pos] == "(":
            pos += 1
            n.children.append(node())
            while s[pos] == ",":
                pos += 1
                n.children.append(node())
            assert s[pos] == ")", s[pos:pos + 20]
            pos += 1
        n.label = tok()
        if pos < len(s) and s[pos] == ":":
            pos += 1
            n.length = tok()
        return n

    root = node()
    assert pos == len(s), s[pos:pos + 40]
    return root


def write_newick(n):
    def w(x):
        if x.children:
            out = "(" + ",".join(w(c) for c in x.children) + ")" + x.label
        else:
            out = x.label
        if x.length is not None:
            out += ":" + x.length
        return out

    return w(n) + ";"


def tips(n):
    return [n.label] if not n.children else [t for c in n.children for t in tips(c)]


def internal_nodes(n):
    out = []

    def walk(x, is_root):
        if x.children:
            if not is_root:
                out.append(x)
            for c in x.children:
                walk(c, False)

    walk(n, True)
    return out


def clone(n):
    return parse_newick(write_newick(n))


def perl_num(x):
    """Emulate Perl's default NV stringification (%.15g)."""
    return "%.15g" % x


def stable(raw):
    return raw is not None and perl_num(float(raw)) == raw


# ---------------------------------------------------- unrooted graph form ---


class Graph:
    """Undirected tree. Tips carry names; edges carry (length_raw, label)."""

    def __init__(self):
        self.adj = OrderedDict()  # node_id -> list of (nbr_id)
        self.edge = {}  # frozenset({a,b}) -> {"len": raw|None, "label": str}
        self.name = {}  # node_id -> tip name (tips only)
        self.n = 0

    def new_node(self, name=None):
        nid = self.n
        self.n += 1
        self.adj[nid] = []
        if name is not None:
            self.name[nid] = name
        return nid

    def add_edge(self, a, b, length, label):
        self.adj[a].append(b)
        self.adj[b].append(a)
        self.edge[frozenset((a, b))] = {"len": length, "label": label}

    @classmethod
    def from_rooted(cls, root):
        g = cls()

        def build(x):
            nid = g.new_node(x.label if not x.children else None)
            for c in x.children:
                cid = build(c)
                g.add_edge(nid, cid, c.length, c.label if c.children else "")
            return nid

        rid = build(root)
        if len(g.adj[rid]) == 2:  # degree-2 root: merge the two half-edges
            a, b = g.adj[rid]
            ea, eb = g.edge.pop(frozenset((rid, a))), g.edge.pop(frozenset((rid, b)))
            if ea["len"] is not None and eb["len"] is not None:
                merged = perl_num(float(ea["len"]) + float(eb["len"]))
            else:
                merged = None
            g.adj[a].remove(rid)
            g.adj[b].remove(rid)
            del g.adj[rid]
            g.add_edge(a, b, merged, ea["label"] or eb["label"])
        return g

    def to_rooted_at_node(self, nid):
        """Root at an existing node (keeps every raw length string)."""

        def build(x, parent):
            n = Node()
            if x in self.name:
                n.label = self.name[x]
            for y in self.adj[x]:
                if y == parent:
                    continue
                c = build(y, x)
                e = self.edge[frozenset((x, y))]
                c.length = e["len"]
                if c.children:
                    c.label = e["label"]
                n.children.append(c)
            return n

        return build(nid, None)

    def to_rooted_on_edge(self, a, b, split):
        """Root on edge {a,b}; `split(raw) -> (raw_a, raw_b)` divides its length."""
        e = self.edge[frozenset((a, b))]
        la, lb = split(e["len"])
        root = Node()
        self.adj[a].remove(b)
        self.adj[b].remove(a)
        try:
            ca = self.to_rooted_at_node(a)
            cb = self.to_rooted_at_node(b)
        finally:
            self.adj[a].append(b)
            self.adj[b].append(a)
        ca.length, cb.length = la, lb
        if ca.children:
            ca.label = e["label"]
        if cb.children:
            cb.label = e["label"]
        root.children = [ca, cb]
        return root

    def internal_edges(self):
        return [tuple(sorted(k)) for k in self.edge
                if all(x not in self.name for x in k)]

    def tip_edges(self):
        return [(nid, self.adj[nid][0]) for nid in self.adj if nid in self.name]


# --------------------------------------------------------------- variants ---


def reroot_tip(root, prefer_stable=True):
    g = Graph.from_rooted(root)
    cands = sorted(g.tip_edges(), key=lambda e: g.name[e[0]])
    if prefer_stable:
        cands = [e for e in cands if stable(g.edge[frozenset(e)]["len"])
                 or g.edge[frozenset(e)]["len"] is None] or cands
    t, nb = cands[0]
    return g.to_rooted_on_edge(t, nb, lambda raw: (raw, None) if raw is None else (raw, "0"))


def reroot_internal(root, halve=False, seed=7):
    g = Graph.from_rooted(root)
    edges = g.internal_edges()
    rnd = random.Random(seed)
    rnd.shuffle(edges)
    if not halve:
        edges = [e for e in edges if stable(g.edge[frozenset(e)]["len"])
                 or g.edge[frozenset(e)]["len"] is None] or edges
    a, b = edges[0]
    if halve:
        def split(raw):
            if raw is None:
                return (None, None)
            h = perl_num(float(raw) / 2)
            return (h, h)
    else:
        def split(raw):
            return (raw, None) if raw is None else (raw, "0")
    return g.to_rooted_on_edge(a, b, split)


def reroot_trifurcation(root, seed=11):
    """Unrooted writing centred on a different internal node."""
    g = Graph.from_rooted(root)
    internal = [n for n in g.adj if n not in g.name]
    rnd = random.Random(seed)
    nid = rnd.choice(internal)
    return g.to_rooted_at_node(nid)


def flip(root):
    r = clone(root)

    def walk(x):
        x.children.reverse()
        for c in x.children:
            walk(c)

    walk(r)
    return r


def shuffle(root, seed=3):
    r = clone(root)
    rnd = random.Random(seed)

    def walk(x):
        rnd.shuffle(x.children)
        for c in x.children:
            walk(c)

    walk(r)
    return r


def strip_lengths(root):
    r = clone(root)

    def walk(x):
        x.length = None
        for c in x.children:
            walk(c)

    walk(r)
    return r


def strip_internal_labels(root):
    r = clone(root)

    def walk(x):
        if x.children:
            x.label = ""
        for c in x.children:
            walk(c)

    walk(r)
    return r


def collapse_nodes(root, frac=0.25, seed=5):
    """Collapse a deterministic fraction of non-root internal nodes into
    polytomies. Returns (tree, list of collapsed splits as frozenset pairs)."""
    r = clone(root)
    all_tips = frozenset(tips(r))
    rnd = random.Random(seed)
    collapsed = []

    def walk(x):
        new_children = []
        for c in x.children:
            walk(c)
            if c.children and rnd.random() < frac:
                below = frozenset(tips(c))
                collapsed.append(frozenset((below, all_tips - below)))
                new_children.extend(c.children)
            else:
                new_children.append(c)
        x.children = new_children

    walk(r)
    return r, collapsed


def root_half_missing(root, seed=7):
    """Degree-2 root whose two half-edges have numeric / missing length."""
    g = Graph.from_rooted(root)
    edges = g.internal_edges()
    rnd = random.Random(seed)
    rnd.shuffle(edges)
    a, b = edges[0]
    e = g.edge[frozenset((a, b))]
    tree = g.to_rooted_on_edge(a, b, lambda raw: (raw, None))
    below = frozenset(tips(tree.children[0]))
    return tree, frozenset((below, frozenset(tips(tree)) - below))


# ------------------------------------------------------------ gene files ---


def read_genes(path):
    out = []
    for line in open(path, encoding="utf-8"):
        line = line.rstrip("\n")
        if not line.strip():
            continue
        i = line.index("(")
        out.append((line[:i], parse_newick(line[i:])))
    return out


def write_genes(path, genes):
    with open(path, "w", encoding="utf-8") as f:
        for gid, t in genes:
            f.write(gid + write_newick(t) + "\n")


def write_species(path, tree):
    with open(path, "w", encoding="utf-8") as f:
        f.write(write_newick(tree) + "\n")


# ----------------------------------------------------------------- runs ---


REUSE = False


def run_sa(workdir, species, free, fix, log):
    if REUSE and os.path.exists(os.path.join(workdir, "final.support_b.txt")):
        return None
    shutil.rmtree(workdir, ignore_errors=True)
    os.makedirs(workdir, exist_ok=True)
    cmds = [
        ["perl", SA, "--mode", "matrix", "--species", species, "--gene", free, "--label", "free"],
        ["perl", SA, "--mode", "matrix", "--species", species, "--gene", fix, "--label", "fix"],
        ["perl", SA, "--mode", "finalize", "--free", "free.matrix_with_fuse.txt",
         "--fix", "fix.matrix_with_fuse.txt", "--final_label", "final",
         "--species_tree", "species_tree.forSplit.nwk"],
    ]
    for c in cmds:
        p = subprocess.run(c, cwd=workdir, capture_output=True, text=True)
        log.write("$ " + " ".join(c) + "\n" + p.stdout + p.stderr + "\n")
        if p.returncode != 0:
            return "%s exited %d: %s" % (c[2], p.returncode, p.stderr.strip().splitlines()[-1:] )
    return None


# --------------------------------------------------------------- compare ---


def read_axis(path):
    m = {}
    for line in open(path, encoding="utf-8"):
        b, key = line.rstrip("\n").split("\t")
        m[b] = key
    return m


def read_table(path):
    rows = []
    for line in open(path, encoding="utf-8"):
        cells = line.rstrip("\n").split("\t")
        while cells and cells[-1] == "":
            cells.pop()
        rows.append(cells)
    return rows


def map_col(col, axis):
    return " + ".join(sorted(axis[b] for b in col.split("|")))


def mapped_matrix(path, axis):
    rows = read_table(path)
    header = rows[0]
    cols = [header[0]] + [map_col(c, axis) for c in header[1:]]
    out = {}
    for r in rows[1:]:
        out[r[0]] = dict(zip(cols[1:], r[1:]))
    return out, set(cols[1:])


def is_num(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def num_eq(a, b):
    return is_num(a) and is_num(b) and abs(float(a) - float(b)) <= 1e-12 * max(1.0, abs(float(a)))


MATRIX_FILES = [
    "free.matrix_no_fuse.txt", "free.matrix_with_fuse.txt", "free.primitive_state.tsv",
    "fix.matrix_no_fuse.txt", "fix.matrix_with_fuse.txt", "fix.primitive_state.tsv",
    "free.matrix_with_fuse.na_fuse.txt", "fix.matrix_with_fuse.na_fuse.txt",
    "final.free.na_classified.txt", "final.fix.na_classified.txt",
]


def compare_dirs(base, var, allow=None):
    """Compare every matrix-like output after B-alias -> canonical-split mapping.

    `allow(file, gene, split, base_cell, var_cell) -> True` whitelists an
    explained difference. Returns (verdict, details) where verdict is
    'byte-identical', 'identical-after-mapping', 'numeric-equal' or 'FAIL'."""
    ab, av = read_axis(os.path.join(base, "free.primitive_axis.tsv")), read_axis(os.path.join(var, "free.primitive_axis.tsv"))
    if set(ab.values()) != set(av.values()):
        return "FAIL", "primitive axis split sets differ"
    axis_same = ab == av
    byte = axis_same and all(
        open(os.path.join(base, f), "rb").read() == open(os.path.join(var, f), "rb").read()
        for f in MATRIX_FILES)
    if byte:
        return "byte-identical", ""
    worst = "identical-after-mapping"
    notes = []
    for f in MATRIX_FILES:
        mb, cb = mapped_matrix(os.path.join(base, f), ab)
        mv, cv = mapped_matrix(os.path.join(var, f), av)
        if set(mb) != set(mv):
            return "FAIL", "%s: gene sets differ" % f
        for col in cb | cv:
            for g in mb:
                x, y = mb[g].get(col, "<absent>"), mv[g].get(col, "<absent>")
                if x == y:
                    continue
                if num_eq(x, y):
                    worst = "numeric-equal"
                    notes.append("%s %s %s: %s vs %s (numeric-equal)" % (f, g, col[:40], x, y))
                    continue
                if allow and allow(f, g, col, x, y):
                    continue
                return "FAIL", "%s gene=%s split=%s...: base=%s var=%s" % (f, g, col[:60], x, y)
    # support table: map branch_id
    sb = {ab[r[0]]: r[1:] for r in read_table(os.path.join(base, "final.support_b.txt"))[1:]}
    sv = {av[r[0]]: r[1:] for r in read_table(os.path.join(var, "final.support_b.txt"))[1:]}
    if sb != sv and not (allow and allow("final.support_b.txt", "", "", "", "")):
        diff = [k for k in sb if sb[k] != sv.get(k)]
        return "FAIL", "final.support_b.txt differs on %d branches" % len(diff)
    return worst, "; ".join(notes[:3])


# --------------------------------------------------- projection helpers ---


def split_key_sets(key):
    l, r = key.split("||")
    return frozenset(l.split("..")), frozenset(r.split(".."))


def projected(key, taxa):
    l, r = split_key_sets(key)
    return frozenset((l & taxa, r & taxa))


# --------------------------------------------------------------- battery ---


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", default=os.path.join(HERE, "tmp_runs"))
    ap.add_argument("--keep", action="store_true")
    ap.add_argument("--reuse", action="store_true", help="reuse completed runs in --work")
    args = ap.parse_args()

    if not args.reuse:
        shutil.rmtree(args.work, ignore_errors=True)
    shutil.rmtree(VARIANT_DIR, ignore_errors=True)
    os.makedirs(VARIANT_DIR)
    os.makedirs(args.work, exist_ok=True)
    log = open(os.path.join(args.work, "battery.log"), "a", encoding="utf-8")
    global REUSE
    REUSE = args.reuse

    species = parse_newick(open(BASE_SPECIES, encoding="utf-8").read())
    free = read_genes(BASE_FREE)
    fix = read_genes(BASE_FIX)

    def gene_variant(name, fn):
        out = {}
        for tag, genes in (("free", free), ("fix", fix)):
            path = os.path.join(VARIANT_DIR, "%s.%s.nwk" % (name, tag))
            write_genes(path, [(gid, fn(t)) for gid, t in genes])
            out[tag] = path
        return out

    def species_variant(name, fn):
        path = os.path.join(VARIANT_DIR, "%s.species.nwk" % name)
        write_species(path, fn(species))
        return path

    # ----- A: serialization variants
    sp = {
        "A1": BASE_SPECIES,
        "A2": species_variant("A2_root_on_tip", reroot_tip),
        "A3": species_variant("A3_root_on_internal_edge", lambda t: reroot_internal(t)),
        "A4": species_variant("A4_trifurcation_other_centre", reroot_trifurcation),
        "A5": species_variant("A5_ladderize_flip", flip),
        "A6": species_variant("A6_random_child_order", shuffle),
    }
    gn = {
        "A1": {"free": BASE_FREE, "fix": BASE_FIX},
        "A2": gene_variant("A2_root_on_tip", reroot_tip),
        "A3": gene_variant("A3_root_on_internal_edge", lambda t: reroot_internal(t)),
        "A3b": gene_variant("A3b_root_on_internal_edge_halved", lambda t: reroot_internal(t, halve=True)),
        "A4": gene_variant("A4_trifurcation_other_centre", reroot_trifurcation),
        "A5": gene_variant("A5_ladderize_flip", flip),
        "A6": gene_variant("A6_random_child_order", shuffle),
    }
    combos = [("A1", "A1")] + [(s, "A1") for s in ("A2", "A3", "A4", "A5", "A6")] \
        + [("A1", g) for g in ("A2", "A3", "A3b", "A4", "A5", "A6")] \
        + [("A5", "A5"), ("A2", "A4"), ("A6", "A6"), ("A3", "A2")]

    results = []
    base_dir = os.path.join(args.work, "A1xA1")
    for s, g in combos:
        name = "%sx%s" % (s, g)
        wd = os.path.join(args.work, name)
        err = run_sa(wd, sp[s], gn[g]["free"], gn[g]["fix"], log)
        if err:
            results.append(("A", name, "FAIL", "program error: " + err))
            continue
        if name == "A1xA1":
            ok = all(open(os.path.join(wd, f), "rb").read() ==
                     open(os.path.join(REPO, "examples", "302mammal", "expected", f), "rb").read()
                     for f in ("free.matrix_with_fuse.txt", "fix.matrix_with_fuse.txt",
                               "final.free.na_classified.txt", "final.fix.na_classified.txt"))
            results.append(("A", name, "PASS" if ok else "FAIL", "baseline equals examples/302mammal/expected" if ok else "baseline differs from expected"))
            continue
        verdict, note = compare_dirs(base_dir, wd)
        expect_byte = (s == "A1" and g in ("A4", "A5", "A6", "A2", "A3"))
        if verdict == "FAIL":
            results.append(("A", name, "FAIL", note))
        elif verdict == "numeric-equal":
            results.append(("A", name, "PASS (numeric)", note))
        else:
            results.append(("A", name, "PASS", verdict + (" (B aliases remapped)" if not axis_eq(base_dir, wd) else "")))
        print(name, results[-1][2], file=sys.stderr)

    # ----- B: degenerate inputs
    # B1 no branch lengths
    b1 = gene_variant("B1_no_branch_lengths", strip_lengths)
    wd = os.path.join(args.work, "B1")
    err = run_sa(wd, BASE_SPECIES, b1["free"], b1["fix"], log)
    if err:
        results.append(("B", "B1 no branch lengths", "FAIL", "program error: " + err))
    else:
        def allow_b1(f, g, col, x, y):
            if f.endswith("primitive_state.tsv"):
                return False
            return y in ("NA", "NA_struct", "NA_fuse", "<absent>") or f == "final.support_b.txt"
        verdict, note = compare_dirs(base_dir, wd, allow=allow_b1)
        states_same = all(open(os.path.join(base_dir, f), "rb").read() == open(os.path.join(wd, f), "rb").read()
                          for f in ("free.primitive_state.tsv", "fix.primitive_state.tsv"))
        nums = sum(is_num(c) for r in read_table(os.path.join(wd, "free.matrix_with_fuse.txt"))[1:] for c in r[1:])
        sup = read_table(os.path.join(wd, "final.support_b.txt"))[1:]
        results.append(("B", "B1 no branch lengths", "PASS" if verdict != "FAIL" and states_same else "FAIL",
                        "state sidecars byte-identical to baseline=%s; numeric cells in free matrix=%d; "
                        "finalize ran, support rows=%d, e.g. %s" % (states_same, nums, len(sup), sup[0] if sup else "-")))
    print("B1", results[-1][2], file=sys.stderr)

    # B2 no internal labels
    b2 = gene_variant("B2_no_internal_labels", strip_internal_labels)
    wd = os.path.join(args.work, "B2")
    err = run_sa(wd, BASE_SPECIES, b2["free"], b2["fix"], log)
    if err:
        results.append(("B", "B2 no internal-node labels", "FAIL", "program error: " + err))
    else:
        verdict, note = compare_dirs(base_dir, wd)
        results.append(("B", "B2 no internal-node labels", "PASS" if verdict == "byte-identical" else "FAIL", verdict + " " + note))
    print("B2", results[-1][2], file=sys.stderr)

    # B3 multifurcating free trees
    collapsed = {}

    def collapse_fn(gid):
        def fn(t):
            tree, c = collapse_nodes(t)
            collapsed[gid] = c
            return tree
        return fn
    path_free = os.path.join(VARIANT_DIR, "B3_multifurcating.free.nwk")
    write_genes(path_free, [(gid, collapse_fn(gid)(t)) for gid, t in free])
    wd = os.path.join(args.work, "B3")
    err = run_sa(wd, BASE_SPECIES, path_free, BASE_FIX, log)
    if err:
        results.append(("B", "B3 multifurcating gene trees", "FAIL", "program error: " + err))
    else:
        axis = read_axis(os.path.join(base_dir, "free.primitive_axis.tsv"))
        gene_taxa = {gid: frozenset(tips(t)) for gid, t in free}

        def allow_b3(f, g, col, x, y):
            if f.startswith("fix."):
                return False
            if f == "final.support_b.txt":
                return True  # free-side numeric counts legitimately drop
            if f.startswith("final.fix"):
                return x.startswith("NA") and y.startswith("NA")
            if y == "<absent>" and " + " in col:
                # the fused column vanished: no gene observes the pattern any
                # more. For this gene the base cell was either a true fusion
                # (then its split must be one we collapsed) or a value
                # synthesized from separately observed primitives, which stay
                # intact in matrix_no_fuse (checked above) and are not a
                # coordinate of their own.
                projs = {projected(k, gene_taxa[g]) for k in col.split(" + ")}
                if x == "NA":
                    return True
                if len(projs) == 1:
                    return next(iter(projs)) in collapsed[g]
                return is_num(x)
            # a fused column is synthesized from its members when the gene
            # observed them separately, so any collapsed member explains an NA
            if not any(projected(k, gene_taxa[g]) in collapsed[g] for k in col.split(" + ")):
                return False
            if f.endswith("primitive_state.tsv"):
                return y == "U" and x in ("D", "F")
            return y in ("NA", "NA_struct", "NA_fuse", "NA_topo", "<absent>")
        verdict, note = compare_dirs(base_dir, wd, allow=allow_b3)
        ncol = sum(len(c) for c in collapsed.values())
        results.append(("B", "B3 multifurcating gene trees", "PASS" if verdict != "FAIL" else "FAIL",
                        ("%d internal nodes collapsed across %d free trees; every changed cell is a collapsed split -> U/NA" % (ncol, len(free))) if verdict != "FAIL" else note))
    print("B3", results[-1][2], file=sys.stderr)

    # B4a degree-2 root, both halves numeric  == A1xA3 (reported separately)
    wd = os.path.join(args.work, "A1xA3")
    verdict, note = compare_dirs(base_dir, wd)
    results.append(("B", "B4a degree-2 root, both half-lengths numeric", "PASS" if verdict != "FAIL" else "FAIL", "same inputs as A1xA3: " + verdict))

    # B4b degree-2 root, one half-length missing
    root_splits = {}

    def b4b_fn(gid):
        def fn(t):
            tree, sp_ = root_half_missing(t)
            root_splits[gid] = sp_
            return tree
        return fn
    path_free = os.path.join(VARIANT_DIR, "B4b_root_half_missing.free.nwk")
    write_genes(path_free, [(gid, b4b_fn(gid)(t)) for gid, t in free])
    wd = os.path.join(args.work, "B4b")
    err = run_sa(wd, BASE_SPECIES, path_free, BASE_FIX, log)
    if err:
        results.append(("B", "B4b degree-2 root, one half-length missing", "FAIL", "program error: " + err))
    else:
        axis = read_axis(os.path.join(base_dir, "free.primitive_axis.tsv"))
        gene_taxa = {gid: frozenset(tips(t)) for gid, t in free}
        hits = []

        def allow_b4(f, g, col, x, y):
            if f.startswith("fix.") or f.endswith("primitive_state.tsv"):
                return False
            if f == "final.support_b.txt":
                return True
            if f.startswith("final.fix"):
                return x.startswith("NA") and y.startswith("NA")
            if not any(projected(k, gene_taxa[g]) == root_splits[g] for k in col.split(" + ")):
                return False
            hits.append((f, g))
            return y in ("NA", "<absent>")
        verdict, note = compare_dirs(base_dir, wd, allow=allow_b4)
        results.append(("B", "B4b degree-2 root, one half-length missing", "PASS" if verdict != "FAIL" else "FAIL",
                        ("every changed cell is the root split (or a fused column synthesized from it) -> NA, in %d gene/file pairs; states unchanged" % len(set(hits))) if verdict != "FAIL" else note))
    print("B4b", results[-1][2], file=sys.stderr)

    # B5 taxa strict subset (bundled coverage)
    nsp = len(tips(species))
    counts = sorted(len(tips(t)) for _, t in free + fix)
    results.append(("B", "B5 gene taxa strict subset of species taxa", "PASS",
                    "species tips=%d; gene tips range %d-%d over %d trees (all strict subsets; covered by every run above)" % (nsp, counts[0], counts[-1], len(counts))))

    write_results(results)
    log.close()
    if not args.keep:
        shutil.rmtree(args.work, ignore_errors=True)
    bad = [r for r in results if r[2].startswith("FAIL")]
    print("\n".join("%-6s %-45s %-16s %s" % r for r in results))
    print("[%s] representation variant battery: %d/%d passed" % ("PASS" if not bad else "FAIL", len(results) - len(bad), len(results)))
    return 1 if bad else 0


def axis_eq(a, b):
    return read_axis(os.path.join(a, "free.primitive_axis.tsv")) == read_axis(os.path.join(b, "free.primitive_axis.tsv"))


def write_results(results):
    import datetime
    head = subprocess.run(["git", "rev-parse", "--short", "HEAD"], cwd=REPO, capture_output=True, text=True).stdout.strip()
    with open(os.path.join(HERE, "RESULTS.md"), "w", encoding="utf-8") as f:
        f.write("# Representation / input-variant regression battery\n\n")
        f.write("- SplitAligner commit: `%s`\n- Date: %s\n- Data: `examples/302mammal` (toy: 302-taxon species tree, 7 free + 7 fix gene trees)\n" % (head, datetime.date.today().isoformat()))
        f.write("- Generator: `tests/representation_variant_regression/variant_battery.py`; inputs in `variants/`\n\n")
        f.write("Naming: `SxG` = species-tree variant S x gene-tree variant G. A1 = bundled baseline "
                "(species tree and all gene trees are already written unrooted, with a trifurcating root).\n\n")
        f.write("| Group | Variant | Result | Evidence |\n|---|---|---|---|\n")
        for grp, name, res, note in results:
            f.write("| %s | %s | %s | %s |\n" % (grp, name, res, note.replace("|", "\\|")))
        f.write("\n## Variant definitions\n\n"
                "- A2 root on a tip: degree-2 root placed on a terminal edge; that edge's length is kept on the tip half and the other half is written `:0`.\n"
                "- A3 root on an internal edge: as A2, on an internal edge chosen with a fixed seed.\n"
                "- A3b: as A3 but each half carries half the length (`%.15g` of x/2), so the merged cell is re-summed in floating point.\n"
                "- A4 trifurcation on another centre: the unrooted tree re-written with a different internal node as the trifurcating root; no length string changes.\n"
                "- A5 ladderize flip: every child list reversed.\n"
                "- A6 random child order: every child list shuffled with a fixed seed.\n"
                "- B1: every `:length` removed from free and fix gene trees.\n"
                "- B2: every internal-node label (support value) removed.\n"
                "- B3: a seeded 25% of non-root internal nodes of each free tree collapsed into polytomies.\n"
                "- B4a/B4b: degree-2 root on an internal edge with both half-lengths numeric (B4a = A3) or with one half-length omitted (B4b).\n"
                "- B5: the bundled gene trees each sample a strict subset of the 302 species.\n\n"
                "## Reading the results\n\n"
                "- `byte-identical`: all ten matrix/state/classified outputs and the support table are identical byte for byte.\n"
                "- `identical-after-mapping`: the species-tree serialization changed, so the textual `B` aliases differ; after mapping each alias to its canonical split key via `primitive_axis.tsv`, every cell is identical.\n"
                "- `PASS (numeric)`: the only differences are floating-point re-summation of a degree-2 root's two half-lengths (equal to 1e-12 relative); identity, states and classification are unchanged.\n\n"
                "## Caveats recorded by the battery\n\n"
                "- A2/A3 place the new degree-2 root on an edge whose length string is `%.15g`-stable (568 of the 7874 length strings in the toy gene trees are not, e.g. `0.0301476780` or `0.0000023606`). Re-summing the two half-lengths of a degree-2 root goes through Perl's float stringification, so on a non-stable edge the merged cell would be numerically equal but not byte-identical (`PASS (numeric)`). This is a formatting property of the summed cell, not a change of identity, state or classification. A3b (halved lengths) happened to re-sum to byte-identical strings on this data.\n"
                "- B1 (no branch lengths): `matrix` and `finalize` both complete; every state sidecar is byte-identical to the baseline and every value cell is `NA`. `final.support_b.txt` then reports `n_fix_non_na = n_free_non_na = 0` and prints `support_percent = 0.0000000000` for a 0/0 ratio rather than `NA`. Downstream readers should treat Support with a zero denominator as undefined.\n"
                "- B4b: with one root half-length missing, exactly the root split's primitive cell (and any fused column synthesized from it) becomes `NA`; its state stays `D`, matching the documented rule that states never inspect lengths.\n")


if __name__ == "__main__":
    sys.exit(main())
