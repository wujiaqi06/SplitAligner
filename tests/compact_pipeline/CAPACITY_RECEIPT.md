# Compact pipeline capacity and ledger-read receipt

Date: 2026-11-11 (Asia/Tokyo). Candidate-only, bounded fixtures; no full-panel
computation, pressure test, benchmark, or frozen-baseline modification.

Command from the candidate source root:

```sh
LC_ALL=C /usr/bin/prove tests/compact_pipeline/capacity.pl
```

Result: **63 assertions passed**, exit 0, on Mac `/usr/bin/perl` v5.34.1.
The test invocation emitted no warnings. Private fixture directories are removed
by `File::Temp` at normal exit.

Files at this check:

| File | SHA-256 |
|---|---|
| `lib/SplitAligner/CompactPipeline.pm` | `42aca475c40f2bd4f3dbe5cc9e8b1d9eadc5a2031b00624fdc6c85bb643f654c` |
| `tests/compact_pipeline/capacity.pl` | `293a3e016a80ad898ce1e6067e29d3f4a0cd68c79cf824e25a602a3847354c25` |

Coverage:

- Every one of the ten capacity categories receives an independent low-limit
  fixture: genes, index bytes, fusion patterns, catalog bytes, fusion members,
  columns, input-line bytes, ledger-block bytes, output-row bytes, and reference
  bitset payload bytes. Each aborts with its named capacity diagnostic and returns
  no result eligible for the controller's publication step.
- The three source inputs and an existing output sentinel outside private
  staging remain unchanged; no public files are added. Partially staged files
  inside the private work directory are permitted. This is **not** a replacement
  for the separate controller/ownership transaction tests.
- Oversized unterminated lines fail before the callback; exact-size EOF input is
  accepted; newline bytes count towards the limit. Unicode, CRLF, blank lines,
  and a UTF-8 character crossing the 65,536-byte read boundary preserve byte
  offsets and lengths. BOM and malformed UTF-8 are rejected.
- Ledger reads accept the exact offset/length range within a larger file and
  reject wrong hashes, short reads, invalid/over-budget boundaries, absent final
  newline, multiple records, invalid JSON, and index-record identity/count
  mismatches.

No implementation files were modified for this test task. These checks establish
the bounded-fixture behavior above, not a total-RSS bound, throughput claim, or
release acceptance.
