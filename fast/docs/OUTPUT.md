# Read the CLI output

`context.json` defines the mapping from matrix column A0, A1, … to canonical unrooted reference splits, their terminal status and their original reference source nodes. Taxa are ordered by UTF-8 bytes. Split hex words have fixed width and are printed most significant word first. Use the context to interpret axis IDs.

`matrix.tsv` stores one row per input gene and one column per primitive reference split. Its values sum the lengths of the original gene nodes identified by a match. `raw_states.tsv` gives a separate topology state for each cell.

| State | Meaning |
| --- | --- |
| D | The projection contains one reference axis and matches a gene edge. |
| U | The projection contains one reference axis and has no matching gene edge. |
| F | Two or more reference axes share this projection. Read the independent recovery flag and value in the composite ledger. |
| S | The reference edge becomes structural when restricted to the observed taxa. |

Numeric availability is independent of topology. A matched edge with missing source lengths retains state D and an unavailable value. Root-node length is ignored. Two edges incident to a binary representation root can contribute to the same unrooted branch. U records absence of an exact match, including in a multifurcation.

| Scope | What to read |
| --- | --- |
| matrix | Reference context, gene index, raw states and numeric matrix. |
| full | Adds coverage, gene mapping, composite and unmatched ledgers. The free role also supplies support. |
| support | Free-tree topology and support, without numeric conversion or matrix. |
| truth | Supply headerless rows of gene ID followed by zero or more tab-separated reference taxon names with `--taxa`. The coverage ledger describes reference restrictions. |

`coverage.jsonl` stores reference projections and *all* original axis members. `mapping.jsonl` stores original node indices, each match direction and both unmatched lists. `composites.jsonl` records recovered and unrecovered multi-axis projections. For free-tree support, the denominator is the count of eligible singleton projections (`n_dec`). The numerator is the count recovered (`n_rec`). A zero denominator yields NA.

By default, full and support scopes write `unmatched.jsonl` and `unmatched.index.tsv`. Use `--unmatched-splits off` to omit those two files. Explicit `on` is available for full and support. Matrix and truth omit the files by default. The manifest records whether the files were written, and the counts of absent reference projections and unmatched gene edges. Truth uses null counts because it has no input gene topology.

A successful output directory contains `manifest.json` with input hashes, build identity and each data file's hash and size. The requested directory must be new. A processing error can leave a sibling staging directory with `FAILED.json`.
