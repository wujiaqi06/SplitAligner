#!/usr/bin/env bash
set -euo pipefail

TEST_DIR="$(cd "$(dirname "$0")" && pwd)"
ROOT_DIR="$(cd "$TEST_DIR/../.." && pwd)"
SA="$ROOT_DIR/SplitAligner.pl"
TMP="$(mktemp -d "${TMPDIR:-/tmp}/splitaligner_utf8_paths_XXXXXX")"
trap 'rm -rf "$TMP"' EXIT
export LC_ALL=C

fail() { printf '[FAIL] %s\n' "$*" >&2; exit 1; }
pass() { printf '[PASS] %s\n' "$*"; }

expect_failure() {
    local name="$1"
    shift
    if "$@" >"$TMP/$name.stdout" 2>"$TMP/$name.stderr"; then
        fail "$name unexpectedly succeeded"
    fi
}

assert_no_transaction_debris() {
    local dir="$1"
    if find "$dir" -maxdepth 1 -name '.splitaligner-*' -print -quit | grep -q .; then
        fail "transaction debris remains in $dir"
    fi
}

write_ascii_fixture() {
    local dir="$1"
    mkdir -p "$dir"
    cat > "$dir/species.nwk" <<'EOF'
((A,B),(C,D));
EOF
    cat > "$dir/genes.nwk" <<'EOF'
g((A:1,B:2):3,(C:4,D:5):6);
EOF
}

run_matrix() {
    local dir="$1" species="$2" genes="$3" label="$4"
    shift 4
    (
        cd "$dir"
        perl "$SA" --mode matrix --species "$species" --gene "$genes" --label "$label" "$@"
    )
}

# PATH-01: ASCII inputs and label from a Unicode current working directory.
CWD_ASCII="$TMP/日本目录"
write_ascii_fixture "$CWD_ASCII"
run_matrix "$CWD_ASCII" species.nwk genes.nwk x
[[ -f "$CWD_ASCII/x.matrix_with_fuse.txt" ]]
[[ -f "$CWD_ASCII/x.run_manifest.json" ]]
assert_no_transaction_debris "$CWD_ASCII"
pass "PATH-01 matrix mode publishes from a UTF-8 current working directory"

# PATH-01/PATH-06: Unicode identifiers and taxa from another Unicode CWD.
CWD_UNICODE="$TMP/目录_猫"
mkdir -p "$CWD_UNICODE"
python3 - "$CWD_UNICODE" <<'PY'
from pathlib import Path
import sys
d = Path(sys.argv[1])
(d / "species.nwk").write_text("((猫,犬),(牛,馬));\n", encoding="utf-8")
(d / "genes.nwk").write_text("基因_猫_é((猫:1,犬:2):3,(牛:4,馬:5):6);\n", encoding="utf-8")
PY
run_matrix "$CWD_UNICODE" species.nwk genes.nwk x
python3 - "$CWD_UNICODE" <<'PY'
from pathlib import Path
import json, sys
d = Path(sys.argv[1])
gene = "基因_猫_é"
manifest = json.loads((d / "x.run_manifest.json").read_text(encoding="utf-8"))
assert manifest["inputs"]["species_tree"]["path_as_provided"] == "species.nwk"
assert manifest["inputs"]["gene_trees"]["path_as_provided"] == "genes.nwk"
assert manifest["gene_identity"][0]["gene_id"] == gene
assert gene in (d / "x.gene_id_map.tsv").read_text(encoding="utf-8")
assert gene in (d / "x.primitive_state.tsv").read_text(encoding="utf-8")
assert gene in (d / "x.matrix_with_fuse.txt").read_text(encoding="utf-8")
PY
assert_no_transaction_debris "$CWD_UNICODE"
pass "PATH-01/PATH-06 Unicode taxa, gene IDs, and JSON provenance preserve code points"

# PATH-01/PATH-06: Unicode input filenames from an ASCII CWD.
UNICODE_INPUT="$TMP/unicode_input_names"
mkdir -p "$UNICODE_INPUT"
python3 - "$UNICODE_INPUT" <<'PY'
from pathlib import Path
import sys
d = Path(sys.argv[1])
(d / "物种树.nwk").write_text("((A,B),(C,D));\n", encoding="utf-8")
(d / "基因树.nwk").write_text("基因_é((A:1,B:2):3,(C:4,D:5):6);\n", encoding="utf-8")
PY
run_matrix "$UNICODE_INPUT" 物种树.nwk 基因树.nwk x
python3 - "$UNICODE_INPUT" <<'PY'
from pathlib import Path
import json, sys
d = Path(sys.argv[1])
m = json.loads((d / "x.run_manifest.json").read_text(encoding="utf-8"))
assert m["inputs"]["species_tree"]["path_as_provided"] == "物种树.nwk"
assert m["inputs"]["gene_trees"]["path_as_provided"] == "基因树.nwk"
assert m["gene_identity"][0]["gene_id"] == "基因_é"
PY
pass "PATH-01/PATH-06 Unicode input filenames round-trip through provenance"

# PATH-02/03/04/06: Unicode matrix labels compose through all finalize modes.
COMPOSE="$TMP/compose"
write_ascii_fixture "$COMPOSE"
cp "$COMPOSE/genes.nwk" "$COMPOSE/free.nwk"
cp "$COMPOSE/genes.nwk" "$COMPOSE/fix.nwk"
run_matrix "$COMPOSE" species.nwk free.nwk 自由
run_matrix "$COMPOSE" species.nwk fix.nwk 固定 --force
(
    cd "$COMPOSE"
    perl "$SA" --mode finalize \
        --free 自由.matrix_with_fuse.txt \
        --fix 固定.matrix_with_fuse.txt \
        --final_label 结论
)
[[ -f "$COMPOSE/结论.free.na_classified.txt" ]]
[[ -f "$COMPOSE/结论.fix.na_classified.txt" ]]
[[ -f "$COMPOSE/自由.matrix_with_fuse.na_fuse.txt" ]]
[[ -f "$COMPOSE/固定.matrix_with_fuse.na_fuse.txt" ]]
(
    cd "$COMPOSE"
    perl "$SA" --mode finalize_fix \
        --fix 固定.matrix_with_fuse.txt \
        --final_label 固定结论 --force
)
[[ -f "$COMPOSE/固定结论.fix.na_classified.txt" ]]

# Explicit Unicode manifest arguments must resolve identically to auto-discovery.
(
    cd "$COMPOSE"
    perl "$SA" --mode finalize \
        --free 自由.matrix_with_fuse.txt \
        --fix 固定.matrix_with_fuse.txt \
        --free_manifest 自由.run_manifest.json \
        --fix_manifest 固定.run_manifest.json \
        --final_label 显式结论 --force
)
[[ -f "$COMPOSE/显式结论.finalize_manifest.json" ]]

# Support inputs and discovered branch map use matching Unicode basenames.
cp "$COMPOSE/species_tree.forSplit.nwk" "$COMPOSE/支持树.forSplit.nwk"
cp "$COMPOSE/species_tree.branch_map.txt" "$COMPOSE/支持树.branch_map.txt"
(
    cd "$COMPOSE"
    perl "$SA" --mode finalize \
        --free 自由.matrix_with_fuse.txt \
        --fix 固定.matrix_with_fuse.txt \
        --species_tree 支持树.forSplit.nwk \
        --final_label 支持结果 --force
)
[[ -f "$COMPOSE/支持结果.support_b.txt" ]]
[[ -f "$COMPOSE/支持树.support_b.nwk" ]]
python3 - "$COMPOSE" <<'PY'
from pathlib import Path
import json, sys
d = Path(sys.argv[1])
paired = json.loads((d / "显式结论.finalize_manifest.json").read_text(encoding="utf-8"))
fixed = json.loads((d / "固定结论.finalize_manifest.json").read_text(encoding="utf-8"))
support = json.loads((d / "支持结果.finalize_manifest.json").read_text(encoding="utf-8"))
assert paired["final_label"] == "显式结论"
assert fixed["final_label"] == "固定结论"
assert support["final_label"] == "支持结果"
assert "自由.matrix_with_fuse.na_fuse.txt" in paired["outputs"]
assert "固定.matrix_with_fuse.na_fuse.txt" in paired["outputs"]
assert "支持树.support_b.nwk" in support["outputs"]
PY
assert_no_transaction_debris "$COMPOSE"
pass "PATH-02/03/04/06 Unicode basenames compose through paired, fix-only, Support, and manifest paths"

# PATH-05: conflict refusal, successful --force, and failed replacement preserve a clean old run.
TRANSACTION="$TMP/交易目录"
write_ascii_fixture "$TRANSACTION"
run_matrix "$TRANSACTION" species.nwk genes.nwk 交易
python3 - "$TRANSACTION" "$TMP/transaction_initial.json" <<'PY'
from pathlib import Path
import hashlib, json, sys
d = Path(sys.argv[1])
paths = sorted(p for p in d.iterdir() if p.name.startswith("交易") or p.name.startswith("species_tree"))
out = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in paths if p.is_file()}
Path(sys.argv[2]).write_text(json.dumps(out, sort_keys=True), encoding="utf-8")
PY
expect_failure transaction_default_refusal bash -c "cd '$TRANSACTION' && perl '$SA' --mode matrix --species species.nwk --gene genes.nwk --label 交易"
python3 - "$TRANSACTION" "$TMP/transaction_initial.json" <<'PY'
from pathlib import Path
import hashlib, json, sys
d = Path(sys.argv[1]); expected = json.loads(Path(sys.argv[2]).read_text())
actual = {name: hashlib.sha256((d / name).read_bytes()).hexdigest() for name in expected}
assert actual == expected, (actual, expected)
PY
cat > "$TRANSACTION/replacement.nwk" <<'EOF'
g((A:10,B:20):30,(C:40,D:50):60);
EOF
run_matrix "$TRANSACTION" species.nwk replacement.nwk 交易 --force
python3 - "$TRANSACTION" "$TMP/transaction_after_success.json" <<'PY'
from pathlib import Path
import hashlib, json, sys
d = Path(sys.argv[1])
paths = sorted(p for p in d.iterdir() if p.name.startswith("交易") or p.name.startswith("species_tree"))
out = {p.name: hashlib.sha256(p.read_bytes()).hexdigest() for p in paths if p.is_file()}
Path(sys.argv[2]).write_text(json.dumps(out, sort_keys=True), encoding="utf-8")
PY
cat > "$TRANSACTION/duplicate.nwk" <<'EOF'
g((A:1,B:2):3,(C:4,D:5):6);
g((A:7,B:8):9,(C:10,D:11):12);
EOF
expect_failure transaction_failed_force bash -c "cd '$TRANSACTION' && perl '$SA' --mode matrix --species species.nwk --gene duplicate.nwk --label 交易 --force"
python3 - "$TRANSACTION" "$TMP/transaction_after_success.json" <<'PY'
from pathlib import Path
import hashlib, json, sys
d = Path(sys.argv[1]); expected = json.loads(Path(sys.argv[2]).read_text())
actual = {name: hashlib.sha256((d / name).read_bytes()).hexdigest() for name in expected}
assert actual == expected, (actual, expected)
PY
assert_no_transaction_debris "$TRANSACTION"
pass "PATH-05 Unicode-CWD refusal, successful force, failed force preservation, and cleanup"

# PATH-06: malformed UTF-8 scientific input and OS-path octets fail closed.
INVALID="$TMP/invalid_utf8"
mkdir -p "$INVALID"
cp "$CWD_ASCII/species.nwk" "$INVALID/species.nwk"
python3 - "$INVALID/genes.nwk" <<'PY'
from pathlib import Path
import sys
Path(sys.argv[1]).write_bytes(b"g((A,B),(C,D));\n\xff")
PY
expect_failure invalid_scientific_utf8 bash -c "cd '$INVALID' && perl '$SA' --mode matrix --species species.nwk --gene genes.nwk --label x"
[[ ! -e "$INVALID/x.matrix_with_fuse.txt" ]]
expect_failure invalid_filesystem_path_octets perl -I"$ROOT_DIR/lib" -MSplitAligner::TextIO=decode_filesystem_path_utf8 -e 'decode_filesystem_path_utf8(chr(255), "test filesystem path")'
pass "PATH-06 malformed UTF-8 scientific input and exposed path octets fail closed"

printf '[PASS] RECERT-003 UTF-8 filesystem-path regression suite\n'
