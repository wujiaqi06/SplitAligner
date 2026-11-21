# Compact matrix test migration

This is the explicit migration map for the independent small-fixture harness
`tests/compact_pipeline/run.py`. It does not remove or rewrite the frozen
`21e397` tests. Passing this harness is not a substitute for the complete
publication/failure-injection suite, large-panel validation, or release review.

## Running the harness

Supply explicit source paths. The frozen tree is read-only. Baseline runs occur
inside a new evidence directory and never use a test script's fixed output path.

```sh
python3 tests/compact_pipeline/run.py \
  --frozen /absolute/path/to/frozen/splitaligner \
  --work-dir /absolute/path/to/new/evidence \
  --prepare-only

python3 tests/compact_pipeline/run.py \
  --frozen /absolute/path/to/frozen/splitaligner \
  --candidate /absolute/path/to/candidate/splitaligner \
  --work-dir /absolute/path/to/new/evidence \
  --run-name candidate-001
```

Preparation reads literal fixtures with Python AST parsing; it does not import
or execute the original regression scripts. It runs only small examples through
the frozen CLI, records its Perl source hashes, and pins all expected scientific
files and old source/component ledgers. Candidate execution verifies those hashes
before comparison. Every command, seed configuration, exit code, stdout and
stderr is retained. A failed invocation writes a failure receipt and retains its
partial evidence. A subsequent repaired candidate uses a new `--run-name`.

The 14 scientific files compared per matrix run are the six label-level files
(`matrix_no_fuse`, `matrix_with_fuse`, `fusion_cells`, `primitive_state`,
`gene_id_map`, `primitive_axis`), the three aggregate graph files
(`primitive_truth`, `projected_groups`, `branch_patterns`), and the five existing
`species_tree.*` files. Comparisons use the same label, input bytes and view.
Whole and direct are checked separately. Cross-version manifests are compared
field by field for axis, gene identity, coordinate state and projection graph;
the intentionally changed implementation/storage metadata is not falsely
required to be byte-identical.

## Standalone v2 regression entry points

The five historical entry points below now execute v2 assertions; they do not
silently skip because the old per-gene text layout is absent. The actual frozen
`21e397` checkout and archived fixture tarball remain unchanged.

| Entry point | v2 coverage |
|---|---|
| `tests/direct_view_priority_regression/run.py` | Retains all eight known-answer genes in whole/direct views and fatal single-observation overflow; compares the complete compact JSONL byte-for-byte across views and verifies observation references. Optional Fast/R comparisons remain available. |
| `tests/whole_component_fusion_regression/run.py` | Retains all numerical, reference/recovery, unavailable and nested-group assertions. The pooled target follows `component.observation_id` into the normalized observation; independent LSB-first decoding checks the exact `{A,B} | {E,F,G,H}` split, length `15`, and source nodes `[1,4]`. |
| `tests/phase3_contract_regression/run.sh` | Retains all eight named contract groups. Replacement and formerly colliding IDs are checked in both ledger and index, including unique storage keys, exact block hashes and full byte coverage. Rejected inputs leave the entire output directory empty. |
| `tests/legacy_output_layout_regression/run.py` | The original archived v1 run is still rejected intact. The assertion now names the first missing v2 role, `matrix compact source ledger`, and its exact compact path; namespace hints and all old-file fingerprints remain checked. |
| `tests/input_snapshot_regression/run.py` | Four original transient/persistent matrix and two-finalizer cases use the real compact stage instead of a dead helper hook. The entry point also runs the eight-case `compact_pipeline/input_snapshot.py` suite covering before/after, species/gene and new/force publication boundaries. |

For retained evidence, the two fusion Python scripts and the snapshot script
accept `--work-dir /absolute/path/to/a/new/directory`. The shell suite accepts
`SPLITALIGNER_TEST_TMP`; use a fresh directory for each run. Run with normal
Python assertion checks enabled. Source/output alias and ownership checks stay
in `recert005_io_alias_regression` and `recert006_output_ownership_regression`;
scientific fixture success does not substitute for those publication tests.

## Old assertion to replacement assertion

| Frozen assertion | Purpose | Compact harness replacement |
|---|---|---|
| `whole_component_fusion_regression/run.py:83–87`: component row → split row → source nodes | Complete component, matching empirical split, exact normalized token, and original source edges | `validate_ledger()` decodes **every** old observation and every old component; compares bitset key, normalized length string, source-node list and matched/unmatched reference. The pooled `1,4` sources remain checked as part of the complete equality. |
| `direct_view_priority_regression/run.py:106–113`: component TSV bytes equal across views | View choice must not change reference or empirical component evidence | Full compact JSONL equality between whole/direct, together with unchanged primitive matrices, state, truth and active group records; per-view fusion cells remain byte-identical to the frozen output. |
| `phase3_contract_regression/run.sh:80–87`: one split file after force replacing two genes with one | No stale gene data after replacement | `force_boundary()` checks the new index and the entire JSONL both contain exactly `g1`; a failed replacement preserves the complete current run. |
| `phase3_contract_regression/run.sh:97–109`: two split files for `g/a` and `g_a` | Injective gene/storage identity | Both IDs are present in numeric and Unicode fixtures; exact old gene-map equality, index identity, one ledger block per gene, distinct stable IDs, and complete block coverage replace a file-count requirement. |
| `phase3_contract_regression/run.sh:150–157`, `phase3_narrow_correction_regression/run.sh:78–89`: no old split directory on rejection | Validation must not publish partial results | Complete before/after snapshot equality on rejected cases, including modes/types/inodes and temporary directory absence. This covers real compact outputs, rather than testing a directory v2 never creates. |
| `legacy_output_layout_regression` and old `--force` layout protection | Existing incompatible results must survive unchanged | Candidate `--force` is run over a copy of an actual frozen output; explicit layout/new-namespace diagnostic and complete filesystem snapshot equality are required. Historical fixture archives are left untouched. |
| `phase3_metamorphic_regression`: hash seed and gene ordering | Deterministic identity/output independent of hash traversal or input ordering | Multiple explicit seed/perturb pairs compare scientific files, JSONL, index and dictionary. Reordered input compares scientific files excluding input-line-dependent `gene_id_map`; each ledger independently checks exact input order, input lines and byte positions. |
| `numeric_range_regression`, `direct_view_priority_regression` | Numeric acceptance, normalization overflow versus optional composite overflow | Frozen positive/negative fixtures are exercised in the new harness. Direct/whole priority and length tokens are checked by exact scientific output and all-observation comparisons. |
| `fusion_topology_regression` | Member truth versus active fusion status, paired classification and Support | The untouched frozen controller runs `finalize` with Support and `finalize_fix` directly against candidate matrices/manifests/sidecars. All emitted scientific files and accounting fields are compared to old-on-old results. No adapter changes the candidate output before old consumption. |

## Assertions that remain separate

The low-level `recert_blocker_regression/run.sh:193–214` fixture deliberately
supplies an incomplete artificial axis to the old helper. It is not a normal
complete-reference matrix fixture. Its `! -s` condition allows either an empty
file or no file. Keep that regression for the retained old helper; do not claim
the new core's complete-axis invariant supersedes this interface evidence.

Exact publication inventories, adding/deleting/modifying owned files, hard-link
or symlink replacement, source/output aliases, device guards, publication
rollback, and cleanup failure injection remain the responsibility of their
corresponding ownership/publication suites. This harness supplies an actual
compact-output force boundary check; it does not claim to migrate all those
tests by merely checking that an operation failed.

The obsolete `tree_to_splits.pl` subprocess race hook has been replaced in
`input_snapshot_regression/run.py` by an actual compact-stage wrapper after
input freezing. It explicitly checks that the mutation occurred. The separate
`compact_pipeline/input_snapshot.py` suite captures staged outputs before
publication and exercises both boundaries. The cross-version harness's input
hash/location checks alone are not evidence against concurrent input changes.

## Input and provenance boundaries

- The frozen Perl contract **rejects** a UTF-8 BOM. It accepts valid Unicode,
  CRLF, blank lines, and an unterminated last line. BOM is a rejection fixture,
  not a request to adopt another implementation's acceptance range.
- The dictionary is checked against taxa decoded from the frozen reference axis.
  The independent Python adapter uses the documented LSB-first fixed-width
  bitset and chooses the side excluding the smallest present taxon. It does not
  call the new kernel to produce its own expected values.
- JSONL is input-ordered, observation IDs are one-based local array references,
  and `g%08d` keys are assigned only in the sorted gene map/index. These are
  checked separately, including Unicode byte offsets and complete ledger
  coverage with no gaps, overlaps or extra records.
- Complete empirical observations, including ones unused by any reference
  component, are compared to the old gene split ledger. Every source child
  index occurs exactly once among normalized observations.
- This provenance comparison is independent of the production audit reader.
  Audit-reader tests must additionally prove that reconstructed results are
  compared with saved ledger/scientific records, and that tampered records are
  rejected rather than silently replaced by newly computed ones.

No large-panel timing, memory limit, publication failure, platform portability,
or full release claim follows solely from this small-fixture harness.
