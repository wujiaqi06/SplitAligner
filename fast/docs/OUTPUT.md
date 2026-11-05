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
| matrix | Reference context, gene index, raw states, numeric matrix, support table and annotated tree. |
| full | Adds coverage, gene mapping, composite and unmatched ledgers. Both roles include the support table and tree. |
| support | Topology and support for fixed or free trees, with numeric branch-length conversion switched off. |
| truth | Supply headerless rows of gene ID followed by zero or more tab-separated reference taxon names with `--taxa`. The coverage ledger describes reference restrictions. |

`coverage.jsonl` stores reference projections and *all* original axis members. `mapping.jsonl` stores original node indices, each match direction and both unmatched lists. `composites.jsonl` records recovered and unrecovered multi-axis projections. In a free run (`--scope full` or `support`), each composite row includes `group_state`. A shared split displayed on the gene tree is `Mapped`, even if its length is unavailable. A shared split absent from the gene tree is `NA_fuse_topo`. Exact absence also covers unresolved trees; it does not establish a conflicting resolved split. Each member retains raw state F on either outcome. Read the group's `value` once; never allocate it to all members. Fixed-role ledgers report recovery without assigning a free-tree diagnosis. The `matrix` and `truth` scopes do not report group recovery. Turning off unmatched exports does not turn off the composite ledger.

To count five disjoint outcomes per gene, run `python3 examples/count_free_topology.py OUTPUT --unit coordinate` or choose `--unit primitive` (free full/support output). Both count displayed singletons as Mapped, absent singletons as NA_topo and each empty-side reference branch as NA_struct. A displayed or absent group contributes to NA_fuse or NA_fuse_topo respectively. Coordinate mode counts one per group; primitive mode weights by member count. Numeric availability is irrelevant. Group outcomes are not written back to member cells.

## Branch support

Every `--genes` run writes `support.tsv` and `support.nwk` in its output directory. This includes `matrix`, `full` (the default) and `support`, for both `--role fixed` and `--role free` (the default). The taxa-only `truth` scope has no gene topologies to score.

Each TSV row is one unrooted reference branch. For an internal branch, a gene is eligible when the projected split belongs to that branch alone and retains at least two taxa on each side. `n_dec` counts eligible genes and `n_rec` counts how many of those trees contain the split. Branch-length availability does not affect either count. For topology-constrained trees that display all their eligible splits, the fraction is 1 and the percentage is 100.

`support.nwk` uses the `support_percent` column as its internal node labels. Zero recoveries among eligible genes gives the label `0`. When `n_dec` is zero, the table records `NA` and the tree leaves that branch unlabeled. The renderer uses the source nodes of the unrooted axis recorded in `context.json`. If a binary display root divides an internal branch, its two halves show the same percentage from a single TSV row. Tips retain their names and the root has no support label. Numeric reference lengths retain their original spelling. The new tree omits missing-length markers, old internal labels and comments. The input annotations remain available in `context.json`.

From the `fast/` directory, try the seven-tip example in two fresh output directories.

```sh
build/saf --reference examples/support/reference.nwk --genes examples/support/free.tsv --out support-free
build/saf --reference examples/support/reference.nwk --genes examples/support/fixed.tsv --role fixed --out support-fixed
cat support-free/support.nwk
```

```text
(((A:1,B:2)0:3,C:4)50:5,((D:6,E:7)50:8,(F:9,G:10):11)50:12);
```

Both free trees can test the A+B split and neither contains it, giving 0/2. One recovers ABC|DEFG, giving 1/2. Both fixed trees recover these splits, giving 2/2. Neither gene contains G, so the F+G branch has no eligible observation in either run.

By default, full and support scopes write `unmatched.jsonl` and `unmatched.index.tsv`. Use `--unmatched-splits off` to omit those two files. Explicit `on` is available for full and support. Matrix and truth omit the files by default. The manifest records whether the files were written, and the counts of absent reference projections and unmatched gene edges. Truth uses null counts because it has no input gene topology.

`genes.tsv` indexes successfully processed rows. A duplicate gene ID, an unknown or duplicate gene-tree taxon, a tree with fewer than two leaves, or a malformed length token such as `0x10` aborts the run. Numeric runs also reject a finite decimal outside the supported range (e.g. `1e400`, or a nonzero input such as `1e-400` that converts to zero), and overflow when summing two root-adjacent lengths. In support mode and with `--no-numeric`, the parser checks token syntax and skips conversion. Consequently every row in a published gene index has status `OK`.

A successful output directory contains `manifest.json` with each data file's SHA-256 and size. For inputs the digest and size cover the exact bytes read, including the optional UTF-8 BOM, empty lines, CRLF and any final newline. Inputs can be regular files, FIFOs or `/dev/stdin`. When a Linux input pathname is not valid UTF-8, the display path uses replacement characters and the optional `path_bytes_hex` field preserves every original byte of its absolute path. JSON strings escape U+2028 and U+2029. Read TSVs using LF as the record delimiter.

`binary_sha256` identifies the executable. `build_commit` defaults to `unknown`. To record the source revision of a clean checkout, configure a fresh build with `-DSAF_BUILD_COMMIT="$(git rev-parse HEAD)"`.

Choose a new output directory under an existing parent. The CLI publishes the complete directory with Linux atomic no-replace rename. A filesystem that cannot do this reports the error and the staging path. On a processing or publication error, the sibling staging directory can contain `FAILED.json`; partial data there are diagnostic. For I/O fault testing, `--inject-io-failure read|write|close` deliberately fails the run.
