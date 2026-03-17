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

A species-tree branch may lose its discriminating power after taxon pruning and no longer define a meaningful split in the reduced taxon set.

This corresponds to:

- `NA_struct`

### 2. Fused projected branches

Multiple distinct species-tree branches may collapse into the same projected split after pruning absent taxa.

This creates a fused branch pattern rather than a unique one-to-one correspondence.

This corresponds to:

- fused branch labels in the branch map
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

---

## Interpretation of NA states

SplitAligner distinguishes several forms of missingness.

### `NA`
Generic missing value before final classification.

### `NA_struct`
The projected species-tree branch becomes structurally absent after taxon pruning.

### `NA_fuse`
The corresponding signal is not present as a primitive branch, but is captured by a fused projected branch.

### `NA_topo`
The branch is present in the fixed-topology matrix but absent in the free-topology matrix, consistent with topology-induced discordance.

---

## Key conceptual advantage

The key advantage of SplitAligner is that it separates three fundamentally different reasons why a branch may appear to be "missing":

1. the branch is structurally undefined after projection,
2. the signal is absorbed into a fused branch,
3. the signal is absent due to topological discordance.

A naive branch comparison often conflates these cases.
SplitAligner keeps them explicit.

---

## Output interpretation

The final branch matrices provide a standardized branch coordinate system across genes, while preserving biologically meaningful distinctions among different forms of missingness.

This makes SplitAligner particularly suitable for large-scale comparative analyses in which branch-level signals need to be aggregated across many genes under heterogeneous taxon coverage.

---

## Summary

SplitAligner is based on one central rule:

> branch identity must be defined in projected split space, not assumed to survive taxon pruning unchanged.

By making projected split structure explicit, SplitAligner provides a robust framework for branch mapping in phylogenomics.
