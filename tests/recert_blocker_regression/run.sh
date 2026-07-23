#!/usr/bin/env bash
set -euo pipefail

TEST_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(cd "$TEST_DIR/../.." && pwd)"
SA="$ROOT_DIR/SplitAligner.pl"
TMP="$(mktemp -d "${TMPDIR:-/tmp}/splitaligner_recert_blockers_XXXXXX")"
trap 'rm -rf "$TMP"' EXIT
export LC_ALL=C

fail() { echo "[FAIL] $*" >&2; exit 1; }
pass() { echo "[PASS] $*"; }

run_matrix() {
    local dir="$1" label="$2" species="$3" genes="$4"
    (
        cd "$dir"
        perl "$SA" --mode matrix --species "$species" --gene "$genes" --label "$label"
    )
}

run_finalize() {
    local dir="$1" final="$2"
    (
        cd "$dir"
        perl "$SA" --mode finalize \
            --free free.matrix_with_fuse.txt \
            --fix fix.matrix_with_fuse.txt \
            --final_label "$final"
    )
}

expect_failure() {
    local name="$1"
    shift
    if "$@" >"$TMP/$name.stdout" 2>"$TMP/$name.stderr"; then
        fail "$name unexpectedly succeeded"
    fi
}

# RECERT-001-A/B/F/K: direct topology is independent of branch-value availability.
mkdir -p "$TMP/branchless"
cat > "$TMP/branchless/species.nwk" <<'EOF'
((A,B),(C,D));
EOF
cat > "$TMP/branchless/free.nwk" <<'EOF'
g((A,B),(C,D));
EOF
cp "$TMP/branchless/free.nwk" "$TMP/branchless/fix.nwk"
run_matrix "$TMP/branchless" free species.nwk free.nwk
run_matrix "$TMP/branchless" fix species.nwk fix.nwk

python3 - "$TMP/branchless" <<'PY'
from pathlib import Path
import json, sys
d = Path(sys.argv[1])
for label in ("free", "fix"):
    state = (d / f"{label}.primitive_state.tsv").read_text(encoding="utf-8").splitlines()
    assert state[1] == "g\tD\tD\tD\tD\tD", state
    matrix = (d / f"{label}.matrix_no_fuse.txt").read_text(encoding="utf-8").splitlines()
    assert matrix[1] == "g\tNA\tNA\tNA\tNA\tNA", matrix
    manifest = json.loads((d / f"{label}.run_manifest.json").read_text(encoding="utf-8"))
    assert manifest["schema_version"] == "SplitAligner-run-manifest-v3"
    assert manifest["coordinate_state"]["schema_version"] == "SplitAligner-primitive-state-v1"
    assert manifest["coordinate_state"]["gene_record_count"] == 1
    assert manifest["coordinate_state"]["cell_count"] == 5
PY
run_finalize "$TMP/branchless" final
grep -qx $'g\tNA\tNA\tNA\tNA\tNA' "$TMP/branchless/final.fix.na_classified.txt"
grep -qx $'g\tNA\tNA\tNA\tNA\tNA' "$TMP/branchless/final.free.na_classified.txt"
(
    cd "$TMP/branchless"
    perl "$SA" --mode finalize_fix \
        --fix fix.matrix_with_fuse.txt --final_label final_fix --force
)
grep -qx $'g\tNA\tNA\tNA\tNA\tNA' "$TMP/branchless/final_fix.fix.na_classified.txt"
pass "RECERT-001-A/B branchless direct coordinates remain residual NA"

# Partial length removal and numeric spelling must not alter coordinate state.
mkdir -p "$TMP/metamorphic"
cp "$TMP/branchless/species.nwk" "$TMP/metamorphic/species.nwk"
cat > "$TMP/metamorphic/full.nwk" <<'EOF'
g((A:1,B:2):3,(C:4,D:5):6);
EOF
cat > "$TMP/metamorphic/partial.nwk" <<'EOF'
g((A:1,B):3,(C:4,D));
EOF
cat > "$TMP/metamorphic/spelling.nwk" <<'EOF'
g((A:1.0,B:2.00e0):3.,(C:4E+0,D:5.000):6e0);
EOF
cat > "$TMP/metamorphic/none.nwk" <<'EOF'
g((A,B),(C,D));
EOF
for label in full partial spelling none; do
    run_matrix "$TMP/metamorphic" "$label" species.nwk "$label.nwk"
done
cmp "$TMP/metamorphic/full.primitive_state.tsv" "$TMP/metamorphic/partial.primitive_state.tsv"
cmp "$TMP/metamorphic/full.primitive_state.tsv" "$TMP/metamorphic/spelling.primitive_state.tsv"
cmp "$TMP/metamorphic/full.primitive_state.tsv" "$TMP/metamorphic/none.primitive_state.tsv"
pass "RECERT-001-C/K length and numeric spelling metamorphisms preserve state"

# Equivalent unrooted graph serializations and batch context must preserve state.
cat > "$TMP/metamorphic/degree2.nwk" <<'EOF'
g((A:1,B:2):3,(C:4,D:5):6);
EOF
cat > "$TMP/metamorphic/degree3.nwk" <<'EOF'
g(A:1,B:2,(C:4,D:5):9);
EOF
cat > "$TMP/metamorphic/child_order.nwk" <<'EOF'
g((D:5,C:4):6,(B:2,A:1):3);
EOF
cat > "$TMP/metamorphic/rerooted.nwk" <<'EOF'
g(D:5,C:4,(B:2,A:1):9);
EOF
cat > "$TMP/metamorphic/batch.nwk" <<'EOF'
g((A:1,B:2):3,(C:4,D:5):6);
extra((A:7,C:8):9,(B:10,D:11):12);
EOF
for label in degree2 degree3 child_order rerooted batch; do
    run_matrix "$TMP/metamorphic" "$label" species.nwk "$label.nwk"
done
for label in degree3 child_order rerooted; do
    cmp "$TMP/metamorphic/degree2.primitive_state.tsv" "$TMP/metamorphic/$label.primitive_state.tsv"
done
python3 - "$TMP/metamorphic/degree2.primitive_state.tsv" "$TMP/metamorphic/batch.primitive_state.tsv" <<'PY'
from pathlib import Path
import sys

def row(path, gene):
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    header = lines[0]
    found = [line for line in lines[1:] if line.split("\t", 1)[0] == gene]
    assert len(found) == 1, (path, gene, found)
    return header, found[0]

assert row(sys.argv[1], "g") == row(sys.argv[2], "g")
PY
pass "RECERT Gate 3 equivalent rooting, child order, and batch composition preserve state"

# RECERT-001-D/G/H: structural loss and numeric/nonnumeric fusion.
for kind in numeric nonnumeric; do
    mkdir -p "$TMP/fusion_$kind"
    cp "$TMP/branchless/species.nwk" "$TMP/fusion_$kind/species.nwk"
done
cat > "$TMP/fusion_numeric/fix.nwk" <<'EOF'
g((A:1,B:2):3,C:4);
EOF
cp "$TMP/fusion_numeric/fix.nwk" "$TMP/fusion_numeric/free.nwk"
cat > "$TMP/fusion_nonnumeric/fix.nwk" <<'EOF'
g((A,B),C);
EOF
cp "$TMP/fusion_nonnumeric/fix.nwk" "$TMP/fusion_nonnumeric/free.nwk"

for kind in numeric nonnumeric; do
    run_matrix "$TMP/fusion_$kind" free species.nwk free.nwk
    run_matrix "$TMP/fusion_$kind" fix species.nwk fix.nwk
    run_finalize "$TMP/fusion_$kind" final
    (
        cd "$TMP/fusion_$kind"
        perl "$SA" --mode finalize_fix \
            --fix fix.matrix_with_fuse.txt --final_label final_fix --force
    )
done

python3 - "$TMP/fusion_numeric" "$TMP/fusion_nonnumeric" <<'PY'
from pathlib import Path
import sys

def row(path):
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    header = lines[0].split("\t")[1:]
    values = lines[1].split("\t")[1:]
    return dict(zip(header, values))

for dname, numeric in ((sys.argv[1], True), (sys.argv[2], False)):
    d = Path(dname)
    state = row(d / "fix.primitive_state.tsv")
    final = row(d / "final.fix.na_classified.txt")
    final_fix = row(d / "final_fix.fix.na_classified.txt")
    assert "S" in state.values(), state
    assert list(state.values()).count("F") >= 2, state
    for branch, code in state.items():
        if code == "S":
            assert final[branch] == "NA_struct", (branch, final[branch])
            assert final_fix[branch] == "NA_struct", (branch, final_fix[branch])
        if code == "F":
            expected = "NA_fuse" if numeric else "NA"
            assert final[branch] == expected, (branch, final[branch], expected)
            assert final_fix[branch] == expected, (branch, final_fix[branch], expected)
PY
pass "RECERT-001-D/G/H structural and fused states classify from explicit evidence in paired and fix-only modes"

# RECERT-001-I: a matched endpoint-collapsed singleton is U, not D or S.
mkdir -p "$TMP/endpoint/input"
cat > "$TMP/endpoint/species.splits.txt" <<'EOF'
A..B||C..D	B5
EOF
cat > "$TMP/endpoint/axis.tsv" <<'EOF'
B5	A..B||C..D
EOF
cat > "$TMP/endpoint/map.tsv" <<'EOF'
storage_key	gene_id	input_line
g00000001	g	1
EOF
cat > "$TMP/endpoint/input/g00000001.split.txt" <<'EOF'
A||C..D	1
EOF
(
    cd "$TMP/endpoint"
    perl "$ROOT_DIR/scripts/split_branch_label.pl" \
        -i input -j species.splits.txt -a axis.tsv -g map.tsv -o endpoint
)
grep -qx $'g\tU' "$TMP/endpoint/endpoint.primitive_state.tsv"
[[ ! -s "$TMP/endpoint/endpoint_split_branch_label/g00000001.split.txt" ]]
pass "RECERT-001-I endpoint-collapsed singleton remains U"

# RECERT-001-E: numeric fixed D plus free U is NA_topo.
mkdir -p "$TMP/topology"
cp "$TMP/branchless/species.nwk" "$TMP/topology/species.nwk"
cat > "$TMP/topology/fix.nwk" <<'EOF'
g((A:1,B:2):3,(C:4,D:5):6);
EOF
cat > "$TMP/topology/free.nwk" <<'EOF'
g((A:1,C:2):3,(B:4,D:5):6);
EOF
run_matrix "$TMP/topology" free species.nwk free.nwk
run_matrix "$TMP/topology" fix species.nwk fix.nwk
run_finalize "$TMP/topology" final
python3 - "$TMP/topology" <<'PY'
from pathlib import Path
import sys
d = Path(sys.argv[1])
def row(name):
    lines=(d/name).read_text(encoding="utf-8").splitlines()
    return dict(zip(lines[0].split("\t")[1:], lines[1].split("\t")[1:]))
fs=row("fix.primitive_state.tsv")
fr=row("free.primitive_state.tsv")
fm=row("fix.matrix_no_fuse.txt")
out=row("final.free.na_classified.txt")
targets=[b for b in fs if fs[b]=="D" and fr[b]=="U" and fm[b] not in ("NA", "")]
assert targets, (fs, fr, fm)
assert all(out[b]=="NA_topo" for b in targets), (targets, out)
PY
pass "RECERT-001-E topology absence requires fixed D numeric and free U"

# RECERT-001-F: direct topology without a free-side value remains D/residual NA.
mkdir -p "$TMP/direct_nonnumeric"
cp "$TMP/branchless/species.nwk" "$TMP/direct_nonnumeric/species.nwk"
cat > "$TMP/direct_nonnumeric/fix.nwk" <<'EOF'
g((A:1,B:2):3,(C:4,D:5):6);
EOF
cat > "$TMP/direct_nonnumeric/free.nwk" <<'EOF'
g((A,B),(C,D));
EOF
run_matrix "$TMP/direct_nonnumeric" free species.nwk free.nwk
run_matrix "$TMP/direct_nonnumeric" fix species.nwk fix.nwk
run_finalize "$TMP/direct_nonnumeric" final
python3 - "$TMP/direct_nonnumeric" <<'PY'
from pathlib import Path
import sys
d=Path(sys.argv[1])
def row(name):
    lines=(d/name).read_text(encoding="utf-8").splitlines()
    return dict(zip(lines[0].split("\t")[1:], lines[1].split("\t")[1:]))
assert set(row("fix.primitive_state.tsv").values()) == {"D"}
assert set(row("free.primitive_state.tsv").values()) == {"D"}
assert set(row("fix.matrix_no_fuse.txt").values()) != {"NA"}
assert set(row("final.free.na_classified.txt").values()) == {"NA"}
PY
pass "RECERT-001-F direct nonnumeric free cells remain D and residual NA"

# Differing FREE/FIX taxon sets are classified from each side's own state.
mkdir -p "$TMP/asymmetric_taxa"
cat > "$TMP/asymmetric_taxa/species.nwk" <<'EOF'
((A,B),(C,(D,E)));
EOF
cat > "$TMP/asymmetric_taxa/fix.nwk" <<'EOF'
g((A:1,B:2):3,(C:4,D:5):6);
EOF
cat > "$TMP/asymmetric_taxa/free.nwk" <<'EOF'
g((A:1,B:2):3,(C:4,E:5):6);
EOF
run_matrix "$TMP/asymmetric_taxa" free species.nwk free.nwk
run_matrix "$TMP/asymmetric_taxa" fix species.nwk fix.nwk
run_finalize "$TMP/asymmetric_taxa" final
python3 - "$TMP/asymmetric_taxa" <<'PY'
from pathlib import Path
import sys
d=Path(sys.argv[1])
def row(name):
    lines=(d/name).read_text(encoding="utf-8").splitlines()
    return dict(zip(lines[0].split("\t")[1:], lines[1].split("\t")[1:]))
fs=row("fix.primitive_state.tsv"); fr=row("free.primitive_state.tsv")
fo=row("final.fix.na_classified.txt"); fro=row("final.free.na_classified.txt")
fix_s={b for b,v in fs.items() if v=="S"}; free_s={b for b,v in fr.items() if v=="S"}
assert fix_s and free_s and fix_s != free_s, (fs,fr)
assert all(fo[b]=="NA_struct" for b in fix_s), (fix_s,fo)
assert all(fro[b]=="NA_struct" for b in free_s), (free_s,fro)
assert all(fo[b]!="NA_struct" for b in fs if b not in fix_s), (fs,fo)
assert all(fro[b]!="NA_struct" for b in fr if b not in free_s), (fr,fro)
PY
pass "RECERT side-specific states handle differing FREE/FIX taxon sets"

# RECERT-001-J: tampering fails before transactional replacement.
python3 - "$TMP/branchless" <<'PY'
from pathlib import Path
import hashlib, json, sys
d=Path(sys.argv[1])
tracked=[d/"final.fix.na_classified.txt", d/"final.free.na_classified.txt", d/"final.finalize_manifest.json"]
(d/"before.sha256").write_text("\n".join(f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.name}" for p in tracked)+"\n")
PY
cp "$TMP/branchless/fix.primitive_state.tsv" "$TMP/branchless/fix.primitive_state.tsv.saved"
python3 - "$TMP/branchless/fix.primitive_state.tsv" <<'PY'
from pathlib import Path
import sys
p=Path(sys.argv[1]); p.write_text(p.read_text(encoding="utf-8").replace("\tD", "\tX", 1), encoding="utf-8")
PY
expect_failure state_tamper bash -c "cd '$TMP/branchless' && perl '$SA' --mode finalize --free free.matrix_with_fuse.txt --fix fix.matrix_with_fuse.txt --final_label final --force"
mv "$TMP/branchless/fix.primitive_state.tsv.saved" "$TMP/branchless/fix.primitive_state.tsv"
python3 - "$TMP/branchless" <<'PY'
from pathlib import Path
import hashlib, sys
d=Path(sys.argv[1])
before=(d/"before.sha256").read_text().splitlines()
for line in before:
    expected,name=line.split("  ",1)
    assert hashlib.sha256((d/name).read_bytes()).hexdigest()==expected, name
PY

for variant in missing bad_sha bad_axis bad_order bad_gene truncated duplicate unknown_code; do
    cp -R "$TMP/branchless" "$TMP/tamper_$variant"
done
rm "$TMP/tamper_missing/free.primitive_state.tsv"
python3 - "$TMP/tamper_bad_sha/free.run_manifest.json" <<'PY'
from pathlib import Path
import json,sys
p=Path(sys.argv[1]); x=json.loads(p.read_text()); x["coordinate_state"]["sha256"]="0"*64; x["outputs"]["primitive_state"]["sha256"]="0"*64; p.write_text(json.dumps(x,sort_keys=True,indent=3)+"\n")
PY
python3 - \
    "$TMP/tamper_bad_axis" "$TMP/tamper_bad_order" "$TMP/tamper_bad_gene" \
    "$TMP/tamper_truncated" "$TMP/tamper_duplicate" "$TMP/tamper_unknown_code" <<'PY'
from pathlib import Path
import hashlib,json,sys
def patch(d, transform):
    d=Path(d); s=d/"free.primitive_state.tsv"
    lines=s.read_text(encoding="utf-8").splitlines(); transform(lines)
    s.write_text("\n".join(lines)+"\n",encoding="utf-8")
    h=hashlib.sha256(s.read_bytes()).hexdigest(); m=d/"free.run_manifest.json"; x=json.loads(m.read_text())
    x["coordinate_state"]["sha256"]=h; x["outputs"]["primitive_state"]["sha256"]=h
    m.write_text(json.dumps(x,sort_keys=True,indent=3)+"\n")
patch(sys.argv[1], lambda x: x.__setitem__(0, x[0].replace("B1", "B99", 1)))
patch(sys.argv[2], lambda x: x.__setitem__(0, "\t".join([x[0].split("\t")[0], x[0].split("\t")[2], x[0].split("\t")[1], *x[0].split("\t")[3:]])))
patch(sys.argv[3], lambda x: x.__setitem__(1, x[1].replace("g\t", "other\t", 1)))
patch(sys.argv[4], lambda x: x.__setitem__(1, x[1].rsplit("\t", 1)[0]))
patch(sys.argv[5], lambda x: x.append(x[1]))
patch(sys.argv[6], lambda x: x.__setitem__(1, x[1].replace("\tD", "\tX", 1)))
PY
for variant in missing bad_sha bad_axis bad_order bad_gene truncated duplicate unknown_code; do
    expect_failure "tamper_$variant" bash -c "cd '$TMP/tamper_$variant' && perl '$SA' --mode finalize --free free.matrix_with_fuse.txt --fix fix.matrix_with_fuse.txt --final_label tampered"
    [[ ! -e "$TMP/tamper_$variant/tampered.free.na_classified.txt" ]]
done
pass "RECERT-001-J state provenance and transactional tamper checks"

# RECERT-002-A/B/C: decoded Unicode survives every text/provenance boundary.
mkdir -p "$TMP/unicode"
python3 - "$TMP/unicode" <<'PY'
from pathlib import Path
import sys
d=Path(sys.argv[1])
taxa=("猫|A", ".犬", "鳥..β", "𠮷")
d.joinpath("species.nwk").write_text(f"(({taxa[0]},{taxa[1]}),({taxa[2]},{taxa[3]}));\n",encoding="utf-8")
records=(
    f"基因_猫_é(({taxa[0]}:1,{taxa[1]}:2):3,({taxa[2]}:4,{taxa[3]}:5):6);",
    f"遺伝子_𠮷_かな(({taxa[0]}:1,{taxa[1]}:2):3,({taxa[2]}:4,{taxa[3]}:5):6);",
)
d.joinpath("genes.nwk").write_text("\n".join(records)+"\n",encoding="utf-8")
PY
run_matrix "$TMP/unicode" free species.nwk genes.nwk
run_matrix "$TMP/unicode" fix species.nwk genes.nwk
(
    cd "$TMP/unicode"
    perl "$SA" --mode finalize --free free.matrix_with_fuse.txt --fix fix.matrix_with_fuse.txt \
        --species_tree species_tree.forSplit.nwk --final_label final
)
python3 - "$TMP/unicode" <<'PY'
from pathlib import Path
import json,sys
d=Path(sys.argv[1]); genes=["基因_猫_é","遺伝子_𠮷_かな"]; taxa={"猫|A",".犬","鳥..β","𠮷"}
for label in ("free","fix"):
    matrix=(d/f"{label}.matrix_no_fuse.txt").read_text(encoding="utf-8").splitlines()
    state=(d/f"{label}.primitive_state.tsv").read_text(encoding="utf-8").splitlines()
    mapping=(d/f"{label}.gene_id_map.tsv").read_text(encoding="utf-8")
    manifest=json.loads((d/f"{label}.run_manifest.json").read_text(encoding="utf-8"))
    assert [x.split("\t",1)[0] for x in matrix[1:]]==genes
    assert [x.split("\t",1)[0] for x in state[1:]]==genes
    assert all(g in mapping for g in genes)
    assert [x["gene_id"] for x in manifest["gene_identity"]]==genes
final=json.loads((d/"final.finalize_manifest.json").read_text(encoding="utf-8"))
assert final["schema_version"]=="SplitAligner-finalize-manifest-v2"
assert all(g in (d/"final.free.na_classified.txt").read_text(encoding="utf-8") for g in genes)
assert all(t in (d/"species_tree.branch_map.txt").read_text(encoding="utf-8") for t in taxa)
assert all(t in (d/"species_tree.support_b.nwk").read_text(encoding="utf-8") for t in taxa)

def decode_set(key):
    assert key.startswith("HS1:")
    payload=key[4:]
    return [] if not payload else [bytes.fromhex(x).decode("utf-8") for x in payload.split(",")]
seen=set()
for line in (d/"species_tree.primitive_axis.tsv").read_text(encoding="utf-8").splitlines():
    _,key=line.split("\t")
    assert key.startswith("HX1:")
    left,right=key[4:].split("||")
    seen.update(decode_set(left)); seen.update(decode_set(right))
assert seen==taxa,(seen,taxa)
PY
perl -I"$ROOT_DIR/lib" -MSplitAligner::Newick=split_taxon_sets -e \
    'eval { split_taxon_sets("HX1:HS1:ff||HS1:41") }; exit($@ ? 0 : 1)'
pass "RECERT-002-A/B/C Unicode identifiers and structured keys round-trip"

# RECERT-002-D/E: malformed UTF-8 fails closed before publication.
for kind in species gene; do mkdir -p "$TMP/invalid_$kind"; done
python3 - "$TMP/invalid_species" "$TMP/invalid_gene" <<'PY'
from pathlib import Path
import sys
s=Path(sys.argv[1]); g=Path(sys.argv[2])
s.joinpath("species.nwk").write_bytes(b"((A,\xff),(C,D));\n"); s.joinpath("genes.nwk").write_bytes(b"g((A,B),(C,D));\n")
g.joinpath("species.nwk").write_bytes(b"((A,B),(C,D));\n"); g.joinpath("genes.nwk").write_bytes(b"g_\xff((A,B),(C,D));\n")
PY
for kind in species gene; do
    expect_failure "invalid_$kind" bash -c "cd '$TMP/invalid_$kind' && perl '$SA' --mode matrix --species species.nwk --gene genes.nwk --label x"
    [[ ! -e "$TMP/invalid_$kind/x.matrix_with_fuse.txt" ]]
done

for kind in manifest state matrix; do cp -R "$TMP/unicode" "$TMP/invalid_$kind"; done
python3 - "$TMP/invalid_manifest/free.run_manifest.json" <<'PY'
from pathlib import Path
import sys
p=Path(sys.argv[1]); b=p.read_bytes(); p.write_bytes(b[:20]+b"\xff"+b[20:])
PY
python3 - "$TMP/invalid_state" "$TMP/invalid_matrix" <<'PY'
from pathlib import Path
import hashlib,json,sys
for dname,kind in ((sys.argv[1],"state"),(sys.argv[2],"matrix")):
    d=Path(dname); m=d/"free.run_manifest.json"; x=json.loads(m.read_text(encoding="utf-8"))
    if kind=="state":
        p=d/"free.primitive_state.tsv"; b=p.read_bytes(); p.write_bytes(b[:8]+b"\xff"+b[8:]); h=hashlib.sha256(p.read_bytes()).hexdigest()
        x["coordinate_state"]["sha256"]=h; x["outputs"]["primitive_state"]["sha256"]=h
    else:
        p=d/"free.matrix_with_fuse.txt"; b=p.read_bytes(); p.write_bytes(b[:8]+b"\xff"+b[8:]); h=hashlib.sha256(p.read_bytes()).hexdigest()
        x["outputs"]["matrix_with_fuse"]["sha256"]=h; x["coordinate_state"]["matrix_with_fuse_sha256"]=h
    m.write_text(json.dumps(x,sort_keys=True,indent=3,ensure_ascii=False)+"\n",encoding="utf-8")
PY
for kind in manifest state matrix; do
    expect_failure "invalid_consumed_$kind" bash -c "cd '$TMP/invalid_$kind' && perl '$SA' --mode finalize --free free.matrix_with_fuse.txt --fix fix.matrix_with_fuse.txt --final_label invalid"
    [[ ! -e "$TMP/invalid_$kind/invalid.free.na_classified.txt" ]]
done
pass "RECERT-002-D/E malformed UTF-8 fails before publication"

echo "[PASS] RECERT blocker regression suite"
