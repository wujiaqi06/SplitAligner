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
- `NA_topo`: the branch is present in the fixed-tree matrix but absent in the free-tree matrix

---

## 5. Why does finalize mode stop with an error when there are no shared genes?

The final comparison is defined only on genes shared between the fixed-topology and free-topology matrices. If no shared genes exist, the comparison is not meaningful and the program stops.

---

## 6. Do species names need to match exactly?

Yes. Species labels in the gene trees must be consistent with the species labels in the species tree.

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
