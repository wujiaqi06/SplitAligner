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

`coverage.jsonl` stores reference projections and *all* original axis members. `mapping.jsonl` stores original node indices, each match direction and both unmatched lists. `composites.jsonl` records recovered and unrecovered multi-axis projections. In a free run (`--scope full` or `support`), each composite row includes `group_state`. A shared split displayed on the gene tree is `Mapped`, even if its length is unavailable. A shared split absent from the gene tree is `NA_fuse_topo`. Exact absence also covers unresolved trees; it does not establish a conflicting resolved split. Each member retains raw state F on either outcome. Read the group's `value` once; never allocate it to all members. Fixed-role ledgers report recovery without assigning a free-tree diagnosis. The `matrix` and `truth` scopes do not report group recovery. Turning off unmatched exports does not turn off the composite ledger.

To count five disjoint outcomes per gene, run `python3 examples/count_free_topology.py OUTPUT --unit coordinate` or choose `--unit primitive` (free full/support output). Both count displayed singletons as Mapped, absent singletons as NA_topo and each empty-side reference branch as NA_struct. A displayed or absent group contributes to NA_fuse or NA_fuse_topo respectively. Coordinate mode counts one per group; primitive mode weights by member count. Numeric availability is irrelevant. Group outcomes are not written back to member cells.

For free-tree support, the denominator is the count of eligible singleton projections (`n_dec`). The numerator is the count recovered (`n_rec`). A zero denominator yields NA.

By default, full and support scopes write `unmatched.jsonl` and `unmatched.index.tsv`. Use `--unmatched-splits off` to omit those two files. Explicit `on` is available for full and support. Matrix and truth omit the files by default. The manifest records whether the files were written, and the counts of absent reference projections and unmatched gene edges. Truth uses null counts because it has no input gene topology.

A successful output directory contains `manifest.json` with input hashes, build identity and each data file's hash and size. The requested directory must be new. A processing error can leave a sibling staging directory with `FAILED.json`.
