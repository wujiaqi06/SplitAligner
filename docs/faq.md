# SplitAligner FAQ

## 1. What problem does SplitAligner solve?

SplitAligner provides a split-based framework for reconciling gene trees and a species tree when branch identity becomes unstable due to missing taxa, fused branches, and topological discordance.

---

## 2. Why not compare branches directly?

Direct branch-to-branch comparison can fail when taxa are missing. After taxon pruning, multiple species-tree branches may collapse into the same projected split, so branch identity is no longer preserved in a naive one-to-one sense.

---

## 3. What is a fused branch?

A fused branch is a projected branch pattern produced when multiple original species-tree branches become indistinguishable after removing taxa absent from a gene tree.

---

## 4. What is the difference between `NA_struct`, `NA_fuse`, and `NA_topo`?

- `NA_struct`: structural absence after projection
- `NA_fuse`: signal exists only through a fused branch
- `NA_topo`: the branch has numeric fixed-side primitive evidence but is absent in the free-topology gene tree

---

## 5. Why does finalize mode stop with an error when there are no shared genes?

The final comparison is defined only on genes shared between the fixed-topology and free-topology matrices. If no shared genes exist, the comparison is not meaningful and the program stops.

---

## 6. Do species names need to match exactly?

Yes. Every gene-tree taxon must occur exactly in the species tree. Unknown or duplicate taxa are rejected during whole-file preflight before matrix outputs are created.

---

## 7. Is SplitAligner intended for fixed-topology trees, free-topology trees, or both?

Both. The matrix-generation step can be applied to either kind of gene tree input, and the final classification step is designed to compare the two.

---

## 8. What should a gene-tree input line look like?

Each non-empty line should contain a gene identifier followed immediately by a Newick tree.

Example:

```text
GeneA((A:0.1,B:0.2):0.2,(C:0.1,D:0.1):0.1):0.1;
```

This line-based format allows SplitAligner to generate one split file per gene and to preserve gene identities throughout matrix construction.

Gene identifiers are preserved exactly and must be unique. SplitAligner uses injective internal storage keys, so distinct IDs such as `g/a` and `g_a` do not collide.

---

## 9. What happens when an internal branch projects to a `1|k` split?

It is not emitted as an independently observable primitive internal branch. However, it is retained for fused-path bookkeeping. If another branch projects to the same reduced split and the fused coordinate carries numeric signal, the primitive branch is classified as `NA_fuse`.

---

## 10. Why can residual generic `NA` remain after `finalize`?

Residual `NA` can remain whenever the explicit state is nonstructural but the evidence required for a more specific label is unavailable. Examples include a directly mapped branch with no branch length (`D`), a fused mapping with no numeric composite (`F`), or a projected-unmapped free cell (`U`) without qualifying finite numeric fixed-side evidence in state `D`.

---

## 11. Are textual `B` branch IDs stable across equivalent Newick serializations?

Not necessarily. `B` IDs are serialization-local aliases assigned by the historical traversal convention. The authoritative identity is the canonical unrooted split bound to that alias in the ordered `<label>.primitive_axis.tsv` ledger and run manifest.

---

## 12. Why does `finalize` require run manifests?

Two matrices can have identical-looking `B1...Bn` headers while those aliases refer to different canonical splits. `finalize` verifies matrix and primitive-state hashes, requires the complete ordered FREE and FIX ledgers to match, and validates state rows and columns against each matrix. It also verifies an optional Support tree against the same ledger. Missing, state-less v2, or mismatched provenance is a fatal error before classified outputs are created.

---

## 13. What happens if an output label already exists?

SplitAligner exits nonzero by default. `--force` is not a filename-based overwrite switch: it may replace only an intact prior run covered by a valid SplitAligner ownership inventory. An unowned file, legacy or malformed manifest, changed owned file, extra directory descendant, symlink/hardlink substitution, or differing shared artifact is rejected and preserved. For an eligible prior run, the replacement is built in a fresh temporary workspace and published transactionally, so stale per-gene files cannot survive and a failed publication restores the prior output.

---

## 14. Which Newick forms are accepted?

One complete species tree and one complete tree per non-empty gene record are required. Finite branch lengths may use signed integer, decimal, or scientific notation, including `1`, `1.`, `.5`, `-0.1`, and `1e-5`. Numeric terminal labels are valid. Quoted labels and bracket comments are currently rejected explicitly rather than silently parsed.

---

## 15. Can taxon labels contain the old split delimiters?

Yes. Labels containing `..`, `||`, or a single `|`, and labels beginning or ending with `.`, are represented with a versioned `HX1:` hexadecimal byte key so that distinct taxa and splits cannot collide. A single internal period, such as the one in `NC_045512.2`, retains the historical readable split form. Colon, ASCII comma, and ASCII whitespace remain structural Newick delimiters and are not accepted as unquoted taxon content.

---

## 16. Why is `species_tree` rejected as a matrix label?

Matrix mode always writes common artifacts named `species_tree.*`. Using `--label species_tree` would make a label-owned primitive ledger collide with that namespace, so the controller rejects the reserved label before parsing inputs or creating a temporary workspace.

---

## 17. Can input paths, working directories, and output labels contain Unicode?

Yes, on platforms where filesystem paths are exposed as valid UTF-8. SplitAligner strictly decodes both command-line paths and paths returned by the operating system before composing basenames, temporary workspaces, manifests, Support outputs, or transactional destinations. It does not apply Unicode normalization. Malformed UTF-8 path bytes are rejected when the platform exposes them, while file-content hashes remain defined over raw bytes.

---

## 18. Can finalization publish outputs across mounted filesystems?

Not in SplitAligner v1.x. The transaction workspace, rollback workspace, and every destination publication parent must resolve to one filesystem device because publication uses atomic `rename()` operations. The existing destination object's inode device is not the publication-device authority. In particular, `.na_fuse` is written beside each input matrix, so run `finalize` or `finalize_fix` from the matrix filesystem or first copy the complete matrix run into the transaction filesystem. SplitAligner rejects a parent-device mismatch before creating a workdir or running helpers, and `--force` cannot override it. Read-only inputs that do not acquire a derived publication destination may still reside on another device.

---

## 19. How are existing input paths validated?

SplitAligner directly opens the exact supplied pathname before canonicalizing it, requires a regular file, and binds the resolved path to the opened file's device and inode. Hash and size are calculated from that opened object. Invalid forms such as a regular-file component followed by `/..`, `/.`, or a trailing slash therefore fail closed, while a valid directory symlink followed by `/..` retains operating-system POSIX semantics. The same strict lookup is repeated immediately before publication and applies to declared and discovered matrices, manifests, state sidecars, Support trees, and branch maps.
