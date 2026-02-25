# SplitAligner

**SplitAligner** is a split-based gene–species tree reconciliation framework that defines **branch identity on a fixed species-tree backbone** under pervasive **missing taxa** and **gene-tree/species-tree discordance**. It produces standardized **branch-by-gene mapping / branch-length tables** and supports explicit decomposition of missingness into:

- `NA_struct`: structural missingness (degenerate projected split due to taxon coverage)
- `NA_fuse`: fusion-row missingness (signal represented on a composite fused branch)
- `NA_topo`: topology-induced missingness (decisive projected split absent from a free-topology gene tree)

This repository hosts the source code, example commands, and scripts to reproduce the figures in the preprint.

---

## Preprint

**SplitAligner: A Gene–Species Tree Reconciliation Framework Using Split-Based Branch Mapping**  
bioRxiv (2026), Wu.

---

## Key ideas

Given a **species tree** `S` and a **gene tree** `G_g` with a gene-specific taxon set `T_g`:

1. Decompose trees into **splits (bipartitions)**.
2. **Project** each species-tree split onto `T_g` to evaluate whether a branch is informative (decisive) or degenerates (`NA_struct`).
3. Under missing taxa, multiple species-tree branches can become indistinguishable on `T_g`, forming a **fusion group**. SplitAligner reports a composite fused-branch identity (e.g., `Bs1|Bs3`) and marks member branches as `NA_fuse`.
4. Under free-topology gene trees, SplitAligner additionally detects **topology-induced missingness** (`NA_topo`).

---

## Inputs / Outputs

### Inputs
- Species tree in Newick (required for the main pipeline)
- Gene trees in Newick (fixed-topology and/or free-topology)
- Optional: gene lists, branch label files, mapping options (see `options/`)

### Outputs
- **Branch-by-gene mapping / branch-length tables** (TSV)
- Optional: tables including **composite fused branches**
- Branch-wise summary tables, including **Support (%)** and missingness decomposition
- Annotated species trees (e.g., with branch-wise Support)

---

## Quick start (example)

> **NOTE**: The repository includes small example datasets under `examples/` for a quick sanity check.

```bash
# 1) Prepare / normalize gene trees (example; adjust to your pipeline)
Rscript scripts/Change_Tree_format.R

# 2) Run SplitAligner mapping (example command; update filenames)
perl SplitAligner.pl \
  --species examples/species_tree.tre \
  --genes   examples/gene_trees.nwk \
  --matrix  results/gene_branch_matrix.tsv \
  --fusion  none

