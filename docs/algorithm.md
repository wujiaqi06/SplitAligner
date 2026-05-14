# SplitAligner algorithm overview

## Overview

SplitAligner is a split-based gene tree–species tree reconciliation framework designed for robust branch mapping under missing taxa, fused branches, and gene-tree/species-tree discordance.

Its core idea is to define **branch identity on the species-tree backbone** using canonicalized unrooted edge splits, then project that split space onto each gene tree according to the taxa actually present in that gene.

This allows branch mapping to remain well-defined even when naive branch-to-branch comparison breaks down after taxon pruning.

---

## Core principle

In a complete species tree, each branch defines a bipartition of taxa.
SplitAligner represents each branch by this bipartition, written in canonicalized split form.

For each gene tree:

1. taxa absent from that gene are removed from the species-tree split space,
2. the projected split representation is recalculated,
3. projected splits are compared with the observed splits in the gene tree,
4. exact and fused branch correspondences are recorded.

Thus, branch identity is not defined by node order or graphical position, but by its split representation after projection.

---

## Why branch comparison becomes difficult

When taxa are missing, the projected species-tree structure may change in nontrivial ways.

Two major consequences follow:

### 1. Degenerate projected branches

A species-tree branch may lose projected identity after taxon pruning and no longer define an independently observable primitive branch on the reduced taxon set.

The current implementation distinguishes three cases:

- empty projected side
  - the branch has no projected identity and is structurally missing
- internal `>=2|>=2`
  - the projected split remains an independently observable primitive internal branch
- internal `1|k`
  - the projected split is not independently observable as a primitive internal branch, but it is retained for fused-path bookkeeping because adjacent or terminal branches may project to the same reduced split and carry numeric fused signal

This corresponds to:

- `NA_struct` for empty-side loss of projected identity

### 2. Fused projected branches

Multiple distinct species-tree branches may collapse into the same projected split after pruning absent taxa.

This means that the signal is represented through a fused branch rather than through a unique one-to-one correspondence.

If a numeric fused coordinate explains the primitive absence, the primitive branch is classified as `NA_fuse`.

This corresponds to:

- fused-branch labels in the branch map
- `NA_fuse` during downstream matrix interpretation

---

## Workflow

The SplitAligner workflow consists of the following steps.

### Step 1. Label the species tree

The species tree is relabeled so that each branch can be tracked explicitly in downstream processing.

### Step 2. Convert species-tree branches to splits

Each species-tree branch is converted into a canonicalized unrooted edge split.

This defines the reference branch coordinate system.

### Step 3. Convert gene trees to splits

Each gene tree is also converted into split form so that comparison is performed in the same representation.

### Step 4. Project the species-tree split space

For each gene tree, taxa absent from that gene are removed from the species-tree split space before branch comparison.

This projection step is central to the method.

### Step 5. Detect exact and fused correspondences

Projected species-tree splits are compared with observed gene-tree splits.

Possible outcomes include:

- exact mapping,
- fused mapping,
- structural absence,
- topological discordance.

### Step 6. Generate branch matrices

The reconciled mappings are summarized as gene-by-branch matrices for downstream analysis.

These matrices can optionally include fused-branch completion.

### Step 7. Compare fixed-tree and free-tree matrices

In the final stage, fixed-topology and free-topology outputs are compared to classify missing values into biologically meaningful categories.

### Step 8. Compute branch-wise concordance

Optionally, SplitAligner computes a branch-wise concordance score, `Support`, on the species-tree backbone after final classification.

In the current implementation, for each branch `b`, `Support(b)` is computed on genes shared between the fixed-topology and free-topology inputs as:

`100 × n_free_numeric(b) / n_fix_numeric(b)`

where both counts include only numeric branch-length evidence. `NA`, `NA_fuse`, `NA_struct`, `NA_topo`, `NaN`, `Inf`, and empty strings are not counted.

The resulting summary can also be written back onto the species tree as a branch-annotated output for visualization and downstream interpretation.

---

## Interpretation of NA states

SplitAligner distinguishes several forms of missingness.

### `NA`
Generic missing value before final classification.

### `NA_struct`
The projected species-tree branch becomes structurally absent only when a projected side disappears entirely after taxon pruning, so the branch has no projected identity for that gene.

Internal `1|k` projections are treated separately: they are not emitted as independently observable primitive internal branches, but they remain fuse-eligible for downstream bookkeeping.

### `NA_fuse`
The branch is absent as a primitive branch but represented through a fused branch after taxon pruning.

### `NA_topo`
The branch has numeric fixed-side primitive evidence but is absent from the free-topology gene tree, consistent with topology-induced discordance.

---

## Key conceptual advantage

The key advantage of SplitAligner is that it separates three fundamentally different reasons why a branch may appear to be "missing":

1. the branch is structurally undefined after projection,
2. the signal is represented through a fused branch,
3. the signal is absent due to topological discordance.

A naive branch comparison often conflates these cases.
SplitAligner keeps them explicit.

---

## Relation to Other Branch Summaries

SplitAligner is related to other split-based or concordance-style summaries in that it evaluates gene-tree information relative to a fixed species-tree backbone. Its emphasis, however, is different.

Naive branch comparison typically asks whether a gene tree contains a branch that appears to correspond to a species-tree branch. This becomes unreliable when missing taxa alter the projected branch structure, because branch identity itself may no longer be preserved in a one-to-one way after pruning.

By contrast, SplitAligner first asks whether the species-tree branch remains well-defined on the gene-specific taxon set. Only after this projection step does it evaluate whether the corresponding split is mapped, fused, structurally absent, or topologically absent.

Similarly, concordance-oriented summaries are often designed to count supporting versus conflicting signal around species-tree branches. SplitAligner can contribute to that broader goal, but its primary role is to provide a standardized branch coordinate system and an explicit missingness decomposition at the gene-by-branch level.

In short, SplitAligner is not only a support-summary method. It is a projection-aware branch reconciliation framework for constructing comparable gene-by-branch matrices under heterogeneous taxon coverage and gene-tree/species-tree discordance.

---

## Branch-Wise Support

SplitAligner can optionally summarize branch-wise concordance on the species-tree backbone using `Support(b)`, defined in the current implementation as the ratio of free-topology numeric-evidence counts to fixed-topology numeric-evidence counts for the same branch on shared genes.

This makes `Support` a branch-resolved cross-gene concordance summary on the projected branch coordinate system. It is therefore distinct from bootstrap values or posterior support measures, even when written to the species tree in the bootstrap position for compatibility with standard Newick viewers.

---

## Output interpretation

The final branch matrices provide a standardized branch coordinate system across genes, while preserving biologically meaningful distinctions among different forms of missingness.

This makes SplitAligner particularly suitable for large-scale comparative analyses in which branch-level signals need to be aggregated across many genes under heterogeneous taxon coverage.

---

## Summary

SplitAligner is based on one central rule:

> branch identity must be defined in projected split space, not assumed to survive taxon pruning unchanged.

By making projected split structure explicit, SplitAligner provides a robust framework for branch mapping in phylogenomics.
