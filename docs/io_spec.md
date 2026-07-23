# SplitAligner I/O Specification

## Overview

SplitAligner operates in two major modes:

1. `matrix`: generate branch-mapping matrices from a species tree and one gene-tree file
2. `finalize`: compare fixed-topology and free-topology matrices to classify NA states

This document summarizes the expected inputs and the main outputs of each mode.

### Text encoding contract

- Scientific text inputs are strict UTF-8 without a byte-order mark (BOM).
- Each input is decoded once; gene identifiers and taxon labels are retained as Unicode text.
- Text outputs use UTF-8 with LF line endings and no BOM.
- Malformed UTF-8 is fatal before matrix or finalized outputs are transactionally published.
- SHA-256 values in provenance manifests are calculated over raw file bytes.
- JSON is encoded from decoded character strings and parsed back to decoded character strings.
- Command-line filesystem paths and path values returned by the operating system are strictly decoded once before basename, directory, temporary-workspace, transaction, or manifest-discovery operations.
- SplitAligner does not normalize Unicode code points or filesystem normalization forms.

### I/O namespace and input-immutability contract

- Before computation, each mode registers every declared or discovered input and every file or directory it may publish, replace, or remove.
- Existing inputs are opened directly from the exact decoded `path_as_provided` before canonical resolution. SplitAligner does not lexically remove `.` or `..` before this kernel lookup; therefore a valid directory symlink such as `link/../input` retains normal POSIX resolution semantics.
- The directly opened object must be a regular file. Its raw-byte SHA-256, size, device, and inode are obtained from the held handle; the original pathname is then resolved without lexical preprocessing, and that resolved path must identify the same device/inode. A regular-file component followed by `/..` or `/.`, a trailing slash on a regular file, and directory/FIFO/socket/device inputs fail closed rather than being accepted from canonicalization alone.
- Reading, alias checks, and provenance bind that same opened and resolved input object. Inferred manifests are derived beside the resolved matrix, state sidecars beside the resolved manifest, and Support branch maps beside the resolved Support tree; each discovered file undergoes the same strict lookup before its content is used.
- An input may not equal or alias a destination by canonical path, symlink-resolved path, or existing device/inode identity.
- An input may not be located inside a directory destination. Publication destinations must also be mutually disjoint, including file destinations nested below a directory destination.
- Path containment is component-aware; similarly prefixed sibling paths are not treated as ancestors.
- A destination filename is not proof of SplitAligner ownership. Every successful mode writes a versioned ownership inventory in its run or finalization manifest.
- File inventory entries bind semantic role, relative publication path, type, mode, size, and raw-byte SHA-256. Directory entries bind the complete deterministic recursive tree, including descendant type, mode, file hash or symlink target, and a tree-inventory SHA-256.
- `--force` permits replacement only when the corresponding prior ownership manifest is valid and every existing owned object still matches that inventory. Missing, legacy, malformed, altered, hard-linked, type-swapped, or extra-content outputs fail closed.
- Identical shared `species_tree.*`, `.na_fuse`, and Support-tree outputs may be reused without mutation. Differing shared outputs are immutable and are not replaced by `--force`.
- Every registered input is snapshotted by resolved path, raw-byte SHA-256, size, and device/inode during preflight. Immediately before publication, SplitAligner directly opens the exact original `path_as_provided` again, repeats handle-to-resolved-object identity binding, and verifies that its symlink chain, target identity, size, and content remain unchanged.
- Generated publication destinations use a separate lexical construction path because their final objects may not yet exist; destination construction is never used to resolve an existing input.
- A namespace conflict or changed input exits nonzero before publication and preserves prior successful outputs.
- SplitAligner v1.x supports one publication filesystem device per invocation. The transaction device is obtained from the resolved current working directory; each destination device is obtained from the destination's publication parent, or the nearest existing resolved parent directory when the immediate parent does not yet exist. The existing destination object's inode device is not used as publication-device authority.
- Every destination must match the transaction `st_dev`. In `finalize` and `finalize_fix`, the `.na_fuse` destinations beside the matrix inputs are included, so those matrix directories must be on the same device as the finalization working directory.
- Device mismatch is rejected with `[ERROR][Publication device]` after namespace registration but before ownership preflight, temporary-workspace creation, or helper execution. `--force` cannot bypass the check, and publication-parent paths and devices are re-resolved immediately before publication.
- Read-only inputs, including matrix-mode species and gene trees, explicit manifests, state sidecars, and Support inputs, may reside on another device when no publication destination is derived there beyond the checked `.na_fuse` path.

---

## Input Files

### 1. Species tree

- Format: Newick
- One species tree per run
- Species labels must match those used in the gene trees
- Branch lengths and internal node annotations are allowed
- Internal annotations are ignored during split-based mapping
- Numeric-only terminal labels are valid
- Duplicate terminal labels are rejected
- Quoted labels and bracket comments are currently rejected explicitly

Example:

```text
((A,B),(C,D));
```

Accepted forms also include branch lengths and support values, for example:

```text
((A:0.1,B:0.2):0.2,(C:0.1,D:0.1):0.1):0.1;
```

```text
((A:0.1,B:0.2)100:0.2,(C:0.1,D:0.1)95:0.1)100:0.1;
```

### 2. Gene trees

- Format: Newick, one record per line
- Each line begins with a gene identifier followed immediately by a tree
- Species labels must match the species-tree naming convention
- Branch lengths and node support annotations are allowed
- The current implementation assumes the line-based input format used in the example dataset
- Gene IDs are preserved exactly and must be unique
- Each non-empty line must contain exactly one complete tree
- Duplicate or species-tree-unknown taxa are fatal before outputs are created
- Finite branch lengths use one shared grammar accepting signed integer, decimal, and scientific notation, including `1`, `1.0`, `1.`, `.5`, `-0.1`, and `1e-5`

Example:

```text
GeneA((A:0.1,B:0.2):0.2,(C:0.1,D:0.1):0.1):0.1;
GeneB((A:0.2,B:0.1):0.1,(C:0.1,(D:0.1,E:0.2):0.1):0.1):0.1;
```

### Supported use cases

- free-topology gene trees inferred independently for each gene
- fixed-topology gene trees inferred under a species-tree-constrained framework

---

## `matrix` Mode

Command:

```bash
perl SplitAligner.pl --mode matrix \
  --species input/speciesTree302.nwk \
  --gene input/free_tree.examples.nwk \
  --label free
```

Required arguments:

- `--species`: species tree file
- `--gene`: gene-tree file
- `--label`: output prefix
- `species_tree` is a reserved matrix label because it would collide with the fixed common `species_tree.*` artifact namespace
- `--force`: optional ownership-gated transactional replacement; without it, existing label-owned outputs are an error

Main outputs:

- `species_tree.forSplit.nwk`
- `species_tree.FigTree.tre`
- `species_tree.splits.txt`
- `species_tree.branch_map.txt`
- `species_tree.primitive_axis.tsv`
- `<label>_splits/`
- `<label>_split_branch_label/`
- `<label>.matrix_no_fuse.txt`
- `<label>.matrix_with_fuse.txt`
- `<label>.primitive_axis.tsv`
- `<label>.gene_id_map.tsv`
- `<label>.primitive_state.tsv`
- `<label>.run_manifest.json`

Interpretation:

- `species_tree.splits.txt` lists retained canonical species split keys; `species_tree.primitive_axis.tsv` binds their serialization-local `B` aliases in matrix order
- legacy-safe taxa retain the historical readable split-key form; keys involving a taxon label containing `..` or `|`, or beginning or ending with `.`, use the injective versioned `HX1:` hexadecimal byte encoding, while a single internal period remains legacy-safe
- `<label>_splits/` stores per-gene split representations
- `<label>_split_branch_label/` stores per-gene mappings onto the species-tree branch coordinate system
- `<label>.matrix_no_fuse.txt` contains primitive branches only
- `<label>.matrix_with_fuse.txt` additionally includes fused-branch columns such as `B12|B47`
- `<label>.gene_id_map.tsv` maps injective filesystem storage keys back to exact original gene IDs
- `<label>.primitive_state.tsv` uses schema `SplitAligner-primitive-state-v1`; its header is the exact primitive axis and each cell is `S`, `D`, `F`, or `U`
- `<label>.run_manifest.json` uses schema `SplitAligner-run-manifest-v3` and records the canonical-key schema, ordered coordinate mapping and SHA-256, input hashes, original record count, identity map, matrix hashes, the bound primitive-state provenance, and a `SplitAligner-output-ownership-v1` inventory covering every published file and directory descendant
- primitive matrix columns always come from the complete ledger and are never inferred from gene-observed branch patterns

---

## `finalize` Mode

Command:

```bash
perl SplitAligner.pl --mode finalize \
  --free free.matrix_with_fuse.txt \
  --fix fix.matrix_with_fuse.txt \
  --final_label final \
  --species_tree species_tree.forSplit.nwk
```

Required arguments:

- `--free`: free-topology matrix produced by `matrix` mode
- `--fix`: fixed-topology matrix produced by `matrix` mode
- `--final_label`: output prefix

Optional argument:

- `--species_tree`: `species_tree.forSplit.nwk` for branch-wise `Support` calculation and species-tree annotation
- `--free_manifest` and `--fix_manifest`: explicit run manifests when they cannot be inferred beside the matrix files
- `--force`: transactional replacement of existing finalize outputs

Main outputs:

- `<free_input_stem>.na_fuse.txt`
- `<fix_input_stem>.na_fuse.txt`
- `<final_label>.fix.na_classified.txt`
- `<final_label>.free.na_classified.txt`
- `<final_label>.finalize_manifest.json`
- `<final_label>.support_b.txt` if `--species_tree` is provided
- `<species_prefix>.support_b.nwk` if `--species_tree` is provided

Interpretation:

- `extract_na_fuse.pl` marks primitive-branch `NA` cells that are explained by numeric fused-branch signal as `NA_fuse`
- `confirm_na_structure.pl` compares fixed-topology and free-topology matrices on shared genes using each matrix's explicit primitive-state sidecar
- generic `NA` becomes `NA_struct` only for state `S`; free-side state `U` becomes `NA_topo` only when paired with finite numeric fixed-side primitive evidence in state `D`
- nonnumeric `D`, `F`, and nonqualifying `U` states remain residual generic `NA`
- if `--species_tree` is provided, `finalize` also computes branch-wise `Support` and writes an annotated species tree

Examples:

- `free.matrix_with_fuse.txt -> free.matrix_with_fuse.na_fuse.txt`
- `fix.matrix_with_fuse.txt -> fix.matrix_with_fuse.na_fuse.txt`

The final comparison is defined only for genes shared between the fixed-topology and free-topology inputs. If no shared genes are found, SplitAligner stops with an error.

Before any classified output is created, `finalize` verifies each matrix against its v3 manifest, validates the bound primitive-state hash/schema/axis/genes/cells, and compares the exact ordered canonical coordinate mappings. A missing or altered state sidecar, earlier v2 manifest, matrix-hash mismatch, ledger mismatch, or Support-tree mismatch is fatal. Matching textual `B1...Bn` headers alone are not sufficient, and no unsafe header-only fallback is enabled.

### Branch-wise `Support`

When `--species_tree` is provided, SplitAligner computes `Support(b)` for each branch on the species-tree backbone.

- denominator: number of numeric branch-length entries for branch `b` in the fixed-topology matrix on shared genes
- numerator: number of numeric branch-length entries for branch `b` in the free-topology matrix on shared genes

The resulting outputs are:

- `<final_label>.support_b.txt`: branch-wise summary table
  - columns: `branch_id`, `branch_type`, `n_shared_genes`, `n_fix_non_na`, `n_free_non_na`, `support_percent`, `discordance_percent`
  - `n_fix_non_na` and `n_free_non_na` are compatibility column names; in the current implementation they count numeric branch evidence rather than arbitrary non-NA strings
- `<species_prefix>.support_b.nwk`: standard Newick tree with internal-node `Support` values written in the bootstrap position

---

## NA States

### `NA`

Generic missing value before final classification.

### `NA_fuse`

The branch is absent as a primitive branch but represented through a fused branch after taxon pruning.

### `NA_struct`

The projected branch is structurally absent after projection only when a projected side disappears entirely, so the branch has no projected identity and is not evaluable for that gene.

For internal branches, a projected `1|k` split is not emitted as an independently observable primitive internal branch, but it is retained for fused-path bookkeeping. If a numeric fused coordinate explains the primitive absence, the primitive cell is classified as `NA_fuse`, not `NA_struct`.

### `NA_topo`

The branch has numeric fixed-side primitive evidence but is absent from the free-topology gene tree, consistent with topology-induced discordance.

---

## Notes

- Canonical split representation plus the ordered ledger is authoritative branch identity; textual `B` IDs are serialization-local aliases bound to that ledger.
- Missing taxa can collapse multiple species-tree branches into the same projected split.
- Fused branches are explicitly tracked rather than ignored.
- The `finalize` step is intended to separate coverage-driven missingness from discordance-driven missingness.
- Matrix construction is transactional. Without `--force`, existing label-owned outputs are not overwritten. With `--force`, only intact outputs covered by the matching prior ownership inventory are eligible for replacement, and replacement occurs only after a complete new run succeeds. Legacy manifests remain readable where otherwise supported but cannot authorize destructive replacement.
- Cross-device publication is outside the v1.x transaction contract. Run finalization on the matrix filesystem or copy the complete matrix run into the transaction filesystem before finalizing.
