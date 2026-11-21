# Compact Perl pipeline: rejection-condition migration map

Date: 2026-11-11. Scope: the frozen 21e397ada04702330f4b90842197461e1108eb88 matrix path versus the v2 development candidate. This is a source audit and coverage map, not an additional test run or a release decision.

The old standalone helpers remain available and byte-identical to the frozen copies. The new matrix path still invokes the species-mode helper, but replaces the old gene-mode split export, split mapping and matrix assembly with CompactPipeline and CompactKernel. A refusal concerning a hand-edited old intermediate file therefore remains at that helper interface; it is not necessarily a valid new end-to-end input case.

Repository links below point to the current files. Line references describe the dated source audit; frozen historical sources and receipts remain in the development archive.

## Evidence and notation

- **T**: old scripts/tree_to_splits.pl; **S**: old scripts/split_branch_label.pl; **G**: old scripts/generate_branch_matrix.pl; **M**: old SplitAligner.pl.
- **P**: candidate lib/SplitAligner/CompactPipeline.pm; **K**: candidate lib/SplitAligner/CompactKernel.pm. Newick and TextIO rules are shared unless stated otherwise.
- **H** means exercised by the independent compact_pipeline harness. Its candidate-004 receipt (retained in the development archive) contains 5,981 checks, 49 bounded commands and no failed check, bound to the source hashes recorded in that receipt. Later source edits are not retrospectively covered by that run.
- **K-source** identifies an existing assertion in tests/compact_kernel/run.pl; this document task did not rerun it. **Static** means source inspection only, not an independently executed negative probe in H.
- Existing failures and earlier receipts remain preserved. No program test, benchmark or implementation edit was added while making this document.

## 1. External tree and gene input

| Condition | Frozen refusal or acceptance | New matrix path / compatibility | Existing coverage |
| --- | --- | --- | --- |
| Missing species, gene or label option | M:148–151. T:36–43, S:59–75 and G:46–59 separately validate their own arguments and inputs. | Candidate main:149–152 retains the external requirement. Standalone helper contracts remain unchanged. | Static; H does not enumerate missing CLI arguments or missing files. |
| Empty label, path separator or control character in label | M:1180–1183; reserved matrix label species_tree at 1186–1195. | Candidate main retains these validation functions. See the distinct label 0 exception below. | Static; not an H negative fixture. |
| Empty species input | Newick:74–78. | Same species parser, before compact calculation. | Static; H covers species BOM, not an empty species file. |
| No nonblank gene record | Newick:100–102,127; M:247 parses the complete file before staging. | P:139 skips blank lines; P:190 refuses a file with no gene records. This can occur after staging, which is cleaned on failure. | H: reject_empty. |
| No opening parenthesis in gene record | Newick:104–106. | P:140–141. | Static; not a distinct H fixture. |
| Empty gene ID before the tree | Newick:108–110 trims surrounding whitespace and rejects empty ID. | P:142–144 retains this rule. | H: reject_empty_gene. |
| Duplicate gene ID | Newick:111–112, using the trimmed ID and reporting the first line. | P:145; IDs are checked before sorting or assigning g-prefixed storage keys. | H: reject_duplicate_gene. |
| Internal TAB in gene ID | T:123–128 rejects it after the old full-file parse. Surrounding whitespace is trimmed first. | P:146 rejects the internal TAB in the same decoded ID. | H: reject_tab_gene. |
| Slashes, underscores, Unicode or numeric gene IDs | Not a blanket rejection: T:130–143 uses separate safe storage keys. Labels and gene IDs have different rules. | Raw gene ID and input_record_id are retained in the ledger; P:192–196 sorts final gene IDs and assigns storage keys. | H positive fixtures include g/a, g_a, numeric and Unicode IDs and reversed input order. Pure gene ID 0 is not separately exercised. |
| Invalid UTF-8 or initial UTF-8 BOM | TextIO:65–72,94–99. | P:49–53 explicitly performs BOM rejection and strict UTF-8 decoding; species input still uses TextIO. | H: reject_invalid_utf8, reject_bom_gene, reject_bom_species. |
| CRLF, blank lines or absent final LF | TextIO:102–108 accepts these. | P:51–52 removes line terminators; blank lines remain ignored. | H: unicode_crlf (including physical input-line provenance). |
| Duplicate terminal taxon | Newick:155–159. | Same parser; K:139 independently protects direct kernel calls and preserves the original taxon in the diagnostic. | H: reject_duplicate_taxon. K-source:158–160. |
| Taxon absent from species reference | M:248 calls Newick:164–173; the standalone T gene command has no reference and cannot make this check. | K:137–138 rejects the original taxon name; P:150–155 adds gene/line context. | H: reject_unknown asserts Z and gene/line context. The earlier candidate-002 failure is retained; candidate-004 contains the repaired result. K-source:157. |
| Unbalanced tree, missing terminator or a second tree on one record | Newick:147–153,421–474. | Same parser; parsed once per accepted record. | H: reject_unbalanced, reject_two_trees. Other syntax branches are not each isolated. |
| Quoted names or bracket comments | Newick:136–137 explicitly refuses them. | Same parser. | H: reject_quoted, reject_comments. |
| Unary node, including unary display root | Newick:428–442 requires at least two children for each internal node. | Same parser. | H: reject_unary. |
| Nonbinary tree / polytomy | Not a general rejection: the parser has a lower bound of two children, not an upper bound of two. | K:128–132 aggregates all children. | H positive numeric_shapes includes star, two-taxon and three-taxon trees. |
| Fewer than two terminal taxa | Newick:155; K also requires at least two retained taxa. | Same parser; K:144 and materialized-ledger K:188 add internal defenses. | Static for a direct one-tip end-to-end refusal. |
| A space embedded in an unquoted terminal name | Tokenization makes the tested input invalid Newick; this is not a new blanket whitespace-name filter. | Same parser. | H: reject_space_taxon; observed diagnostic is the internal-node child-count error. |
| Missing branch length, without a colon | Accepted; Newick and old gene_split_records:337–344 carry NA. It must not be reclassified as an invalid numeric token. | K:158–160 preserves unavailable length; matrix and ledger retain identity. | H: whole_components, priority and numeric_shapes positive missing-length cases. |
| Colon without a value, or a nonnumeric token | Newick:455–470; numeric syntax is the ASCII rule at 40–49. Species-label mode has its explicit B-number allowance. | Same parser. | H: reject_bad_value and reject_unicode_digit; empty value after colon is not a distinct H fixture. |
| Finite-token conversion overflows or nonzero value underflows to zero | Newick:461–466 rejects these. Negative and exact zero values are accepted; smallest representable subnormal is accepted. | Same parser. | H: reject_overflow, reject_underflow; accepted numeric_shapes includes negative, zero and subnormal tokens. |
| Pooling several input edges into one normalized observation overflows | Old gene_split_records:341–344 calls format_pooled_length, which refuses nonfinite sum (Newick:52–60). | K:158–160 uses the same numerical normalization. Both whole and direct need a valid normalized observation. | H: reject_pooled_whole and reject_pooled_direct. K-source:155. |

These entry-point locations can be inspected in [old main][oldmain], [candidate main][newmain], [Newick parser][newick], [TextIO][textio], [P input traversal][pinput] and [K gene computation][kcompute].

## 2. Artificial axes and legacy text interfaces

| Condition | Frozen interface and actual behavior | New representation / migration boundary | Existing coverage |
| --- | --- | --- | --- |
| Missing or duplicate B label in a labeled species tree | T:53–63 rejects a missing B-number edge label or repeated B label. | New main still creates the labeled species tree and calls this species-mode helper. | Static for hand-edited labeled-tree refusal. |
| Duplicate species split caused by suppressing a display root | T:66–117 retains the lowest B representative and records aliases; S:110–126 similarly warns and keeps the lowest B. This is not an unconditional duplicate-split error at these stages. | Species helper and retained-axis policy remain unchanged. K receives the canonical retained axis, not the pre-deduplication edge list. | H scientific and axis-byte comparisons cover ordinary root suppression; not every malformed duplicate input. |
| Bad retained-axis row, B ID, duplicate B, duplicate canonical split, noncanonical split or empty axis | Provenance:44–70 rejects these on the public axis reader. S:77 and G:61 use it. | P:119 still calls that reader. K:33–75 additionally checks the in-memory reference contract. | Static for axis reader corruptions. K-source includes reference and bitset guards; those were not independently rerun for this document. |
| Nonexistent or inconsistent taxon partition in an axis | Old reader validates each split's syntax; S:144–146 additionally compares its species-split B IDs with the primitive ledger. It does not require the complete terminal coordinate set used by K. | K:53–66 requires two nonempty sides, known nonrepeated taxa, exact taxon coverage and canonical split text. K:71–72 requires every terminal coordinate. | K-source:136–138 exercises removing a terminal coordinate. Static for the other individual axis corruptions in H. |
| Deliberately undersized reference axis | Old S permits the artificial B5-only fixture in tests/recert_blocker_regression/run.sh:193–214; a projected internal singleton may remain U with an empty mapped output. | New complete-axis kernel refuses a missing terminal coordinate and treats a singleton internal projection without its terminal partner as an invariant violation (K:114–115). This is an intentional compact internal-contract difference, not a change to legal main-pipeline scientific input. The old helper fixture remains valid at its unchanged interface. | Existing old test is not deleted or silently reinterpreted. H constructs complete axes; K-source checks missing terminal. |
| Empty/bad-header gene map, malformed row, duplicate storage key, duplicate ID or no entries | S:448–470 and G:206–230 reject these. | New matrix creates its gene map from the validated input and sorted index. It does not accept an external hand-edited map. Old helper entry points keep the rejection. | H checks generated map bytes and unsafe gene-ID/storage-key separation; no corrupt legacy-map negative probes. |
| A legacy split file lacks a map entry, duplicates a storage key, or a mapped gene lacks its split file | S:173–177 and369–371 reject these. | New path no longer scans a directory of per-gene split files. Ledger/index count, ID, offset, byte-length and digest checks replace that inventory boundary. | H checks compact indexed blocks and source coverage. It is not a legacy split-directory corruption suite. |
| Malformed species/gene split text or no valid split rows | S:95–107 and185–208 skip short lines; failed canonicalization logs ERROR and continues. Empty gene splits emit a warning and continue; missing final gene-state rows later fail at373–376. Empty valid species axis fails at142. | Such text is not generated or consumed in the new matrix path. Newick is parsed into components; JSONL validation protects the compact boundary. Do not claim that every old logged ERROR was an immediate process rejection. | Existing old helper tests keep their interface meaning. H invalid Newick probes are not equivalent to all corrupt legacy text rows. |
| Duplicate mapped length for one active fusion observation | S:302–303 rejects multiple pooled numeric records where one is expected. | K records one canonical observation per split. K:199 checks strict unique observation-key order; complete source membership is independently checked. | H compares all observations and source nodes; no independent duplicate-observation ledger mutation in H. |
| Overlapping component ownership | G:118–130 rejects a member already owned by another current component. | K:225–228 verifies stored components against the recomputed reference restriction, including membership, key and observation correspondence. | H compares every component to the frozen component records. K-source:127–128 mutates a component. |
| Invalid global fusion-column member | G:196–202 requires B-number syntax and membership in the axis. It does not separately check a repeated member in that target list. | Production targets come from validated reference components. K:276 also rejects repeated or malformed B members for direct construction calls; this is an extra internal defense. Do not claim that a direct K call rejects every unknown but well-formed B number: a missing owner can instead mean deleted. | Static; H targets are valid reference-derived groups. |

See the unchanged [T helper][thelper], [S helper][shelper], [G helper][ghelper], shared [axis reader][axis], [K reference preparation][kref], [K ledger validation][kledger] and [old undersized-axis test][undersized].

Read-only SHA-256 comparison confirmed these candidate helpers equal the frozen originals:

| File | SHA-256 |
| --- | --- |
| scripts/tree_to_splits.pl | 772ad0f4a1afaf110d58063442f60558a75ebca019745f06deb03fe04aed3b4a |
| scripts/split_branch_label.pl | af2462545dca8852fba74d88f93be76684ae99a821975db9e8d925151527d12a |
| scripts/generate_branch_matrix.pl | 20b36fccf4adc43bbe7b17dba70dbc23dc511ce2144ff5f0262c97c2fe27fd04 |

## 3. Missing matrix values are not input rejections

The following outcomes must not be converted into fatal input errors during migration:

| Condition | Frozen G | Compact K | Coverage |
| --- | --- | --- | --- |
| Target includes deleted members or cuts through a larger current component | 139–152: unavailable, with the relevant reason. | 288: unavailable. | H whole_components and priority fixtures. |
| Structurally representable target has an unrecovered component | 158–159: not_recovered. | 295: not_recovered. | H partial and absent-topology cases. |
| Recovered decomposable target requested in direct view | 160–163: direct_only before missing-length or sum checks. | 297: same priority. | H priority_direct covers missing, overflowing, finite and zero lengths. |
| Whole-view target has a component with no numeric length | 164–165: missing_length. | 299: same outcome. | H lengthless and active_missing cases. |
| Sum of several complete observations exceeds finite range | 167–172: nonfinite_sum; no fabricated numeric result. | 304: same outcome. This differs from normalization of one observation, whose overflow is fatal as described above. | H priority_whole and priority_direct distinguish these branches. |

Targets are collected from all reference fusion patterns, including unrecovered patterns (G:66–78; P:160–170), not only from numeric mapped rows. Whole/direct do not change primitive identity, reference states or Support. H compares the scientific outputs byte-for-byte, components, observations and source membership, and runs old finalize/finalize_fix against the new outputs. This does not establish timing or large-panel memory behavior.

## 4. New compact-boundary defenses

| New check | Location | Coverage classification |
| --- | --- | --- |
| Bounded gene count, index bytes, fusion pattern/member counts, columns, input line, ledger block, output row and reference bitset working estimate | P:20–35 declares defaults; P:47,65,101–105,121–122,147,160–184 checks them. These are explicit v2 capacity refusals, not historical text-format compatibility requirements. | H uses ordinary bounded fixtures; it does not probe every configured limit. Separate capacity evidence is owned by the controller. |
| Invalid ledger offset/length, short read, SHA mismatch, multiple JSON records or gene identity/count mismatch | P:71–98. | H verifies valid indexed blocks; H does not corrupt every index field. Separate audit tests are outside this map's executed evidence. |
| Invalid schema/reference identity, malformed bitset, invalid observation length or source-edge ID, incomplete/repeated source coverage | K:181–224. | H checks valid full source coverage. K-source:129–134 exercises source ID, reference identity and underflowing saved values. |
| Stored components, truth, state or primitive values disagree with deterministic reconstruction | K:225–228 and259–262. | H compares valid records with frozen outputs; K-source:125–128 covers altered truth/component records. |
| Bad in-memory tree structure, cyclic or reused node | K:82–94. | Static; not expressible through the accepted Newick grammar in an ordinary pipeline invocation. |

The ledger reader's verification is an integrity check against the stored ledger, not permission to silently replace it with newly recomputed scientific results.

## 5. Intentional behavior and ownership boundaries

1. **Streaming changes refusal order, not error acceptance.** Old main fully parsed all gene records before validating any taxon membership. New P parses, computes and validates a gene before reading the next. A file containing multiple different errors can therefore report a different first error. No test promises identical precedence for multiple corrupt records. H confirms rejected cases leave no output or private staging residue and preserve existing files.
2. **Complete-axis requirement is explicit.** The old artificial helper axis test remains supported by that helper. It cannot be transferred unchanged to a kernel that requires the full retained reference axis. The corresponding scientific guarantee is retained through complete-axis endpoint and fusion tests, not by deleting the old assertion.
3. **Label 0 is an approved compatibility improvement.** Old main accepts it initially, but downstream S/G use Perl truthiness for their required output label and abort. The controller independently measured old exit 25 with no outputs versus new exit 0; see the archived label-zero receipt. New matrix intentionally accepts this valid label. Standalone S/G remain unchanged and retain their historical rejection; this is distinct from numeric gene IDs and does not change science.
4. **Force is not a legacy-directory migration switch.** The new publication set replaces the old per-gene split directory with the compact ledger and its index/dictionary. H verifies refusal to overwrite the frozen legacy layout, valid replacement of a current compact result, removal of stale gene entries from its index, and preservation of a valid result when the replacement input is invalid. Broader transactional failure injection belongs to the separate publication suite; do not infer it solely from H.
5. **Parsing/IO failures remain failures.** The new bounded path adds write/close and frozen-input verification checks. The three unchanged standalone helper commands retain their old external malformed-file behavior. This map does not promise byte-identical diagnostics, error order or exit codes for all OS failures.

## 6. Explicit coverage limits

The independent H campaign has 14 accepted fixtures and 21 rejected fixtures, plus old finalizer consumption, seed variants and force cases. It does not separately enumerate every missing CLI option, filesystem failure, malformed manually supplied axis/map/component file, parser syntax branch, capacity threshold, audit corruption or mixed-error ordering. The source assertions cited above identify implementation locations; they are not newly claimed executed passes. Windows, full 2,275-gene performance, large-panel memory and publication acceptance are outside this document.

For the positive science/provenance migration and the full fixture list, see [compact-test-migration.md][migration] and the [harness source][harness].

[oldmain]: ../SplitAligner.pl
[newmain]: ../SplitAligner.pl
[newick]: ../lib/SplitAligner/Newick.pm
[textio]: ../lib/SplitAligner/TextIO.pm
[pinput]: ../lib/SplitAligner/CompactPipeline.pm
[kcompute]: ../lib/SplitAligner/CompactKernel.pm
[thelper]: ../scripts/tree_to_splits.pl
[shelper]: ../scripts/split_branch_label.pl
[ghelper]: ../scripts/generate_branch_matrix.pl
[axis]: ../lib/SplitAligner/Provenance.pm
[kref]: ../lib/SplitAligner/CompactKernel.pm
[kledger]: ../lib/SplitAligner/CompactKernel.pm
[undersized]: ../tests/recert_blocker_regression/run.sh
[migration]: ../docs/compact-test-migration.md
[harness]: ../tests/compact_pipeline/run.py
