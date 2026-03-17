# SplitAligner

**SplitAligner** is a split-based gene tree–species tree reconciliation framework for robust branch mapping under missing taxa, fused branches, and gene-tree/species-tree discordance.

It defines **branch identity on a fixed species-tree backbone** using canonicalized unrooted edge splits, projects the species-tree split space onto each gene tree according to its observed taxon set, and generates standardized **gene-by-branch mapping / branch-length matrices** for downstream comparative analyses.

SplitAligner explicitly distinguishes several biologically meaningful forms of missingness:

- `NA_struct`: structural missingness caused by degenerate projected splits after taxon pruning
- `NA_fuse`: fusion-row missingness, where signal is represented on a composite fused branch
- `NA_topo`: topology-induced missingness, where a decisive projected split is absent from a free-topology gene tree

This repository contains the SplitAligner source code, example datasets, and documentation for reproducing the core branch-mapping workflow.

---

## Introduction

On a fixed species-tree spine,

we ask one thing: does branch *b* still hold?

Project the split:

If it collapses — `NA_struct`.

If branches fuse — `Bs1|Bs3`, `NA_fuse`.

If topology turns away — `NA_topo`.

No ghosts, no leaks:

**Total = Mapped + NA_struct + NA_fuse + NA_topo.**

A quiet ledger—where every absence has a name.

---

## Why SplitAligner?

In phylogenomics, branch identity is often treated as if it were stable across all gene trees. In practice, missing taxa can collapse multiple species-tree branches into the same projected split, making naive branch-to-branch comparison unreliable.

SplitAligner addresses this problem by:

- defining branch identity using canonicalized unrooted edge splits,
- projecting species-tree splits onto the taxon set observed in each gene tree,
- distinguishing exact and fused branch correspondences,
- generating branch matrices for large-scale comparative analyses,
- separating structural, fusion-related, and topology-induced missingness.

---

## Preprint

**SplitAligner: A Gene–Species Tree Reconciliation Framework Using Split-Based Branch Mapping**  
bioRxiv (2026), Jiaqi Wu.  
https://doi.org/10.64898/2026.02.24.707838
