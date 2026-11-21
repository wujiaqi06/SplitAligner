# Materialized execution-spool checks

Run from the candidate source root:

```sh
LC_ALL=C /usr/bin/prove tests/materialized_spool/run.pl
```

This standalone suite uses complete source records and the strong materializer
as its expected-value source, then independently compares the spool roundtrip
and exact UTF-8 scientific rows. It covers raw `1.00`, signed zero, absent
lengths, Unicode and embedded CR in accepted gene IDs, two retained taxa,
recovered fusion groups and unrecovered groups whose truth remains `NA_fuse`.
Whole/direct value construction is checked after owner reconstruction. The
writer must leave source values and JSON scalar types unchanged.

Fault cases cover changed bytes with the original authoritative hash, short
reads, missing final LF, malformed UTF-8, unknown/duplicate/missing rows,
invalid widths/types, duplicate or missing component members, owner/state/value
inconsistency, group/owner inconsistency, offset/length limits, identity and
reference bindings, and reader/writer block capacities.

Some cases deliberately recompute the spool hash after mutation to reach the
format and cross-field checks. Such a rebinding is not available to a corrupted
temporary file in production: the authoritative write-time seal is held in
memory. The reader cannot independently prove that an internally consistent
replacement describes its linked JSON record without reading that JSON or
replaying the input; this suite does not make that claim. The permanent audit
reader remains responsible for source replay and scientific-row comparisons.

Observed result on 2026-11-11: **214 assertions passed**. The bounded test receipt
and before/after source hashes are in the candidate's
`evidence/materialized_spool_round1/RESULTS.json`. An initial fixture incorrectly
used explicit `:NA`, which the existing Newick parser rejects; it was corrected
to an absent branch length. No production fix was needed for this suite.

## Whole-controller integrity probes

`cli_integrity.py` generates a process-local `PERL5OPT` wrapper; production
files and APIs contain no test hooks. Run with a new evidence directory:

```sh
python3 tests/materialized_spool/cli_integrity.py --work-dir /path/to/new/evidence
```

One positive fresh-output call makes `JSON::PP::decode` fatal and counts the
internal materializer: three genes must finish with zero JSON decodes, exactly
three materializations, and no call to `materialize_gene`. This assertion is
limited to a fresh matrix run; ownership checks for `--force` legitimately
decode the previous published manifest.

Seven faults are each exercised in a fresh directory and with existing compact
results under `--force`: spool byte change, wrong spool offset, wrong gene-record
binding, source-ledger mutation at the first spool read, and index/dictionary/
ledger mutation after the compact stage returns. The latter three append one
space, so JSON syntax cannot explain the rejection. They must reach the
write-time artifact-seal guard. The mid-read ledger mutation must be rejected
by the accumulated write hash after all three spool reads have finished.

Each failed call must preserve original inputs and all previous output bytes,
permissions and inodes, and remove private staging, spool and input snapshots.
Every command and per-process probe counter is saved. Suite receipts bind the
production source before and after the complete run.

After the writer's UTF-8 byte-flag correction, the final bounded checks passed:

- `run.pl` plus `utf8_bytes.pl`: **225 assertions** (214 codec + 11 byte-invariant checks),
  with evidence in `evidence/materialized_spool_round2/`.
- `cli_integrity.py`: **15 cases**, with evidence in
  `evidence/spool_cli_integrity_round2/`.

Both final receipts confirm unchanged production-source fingerprints throughout
their suites. Earlier rounds remain separate historical evidence.
