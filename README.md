# SplitAligner

**SplitAligner** is a split-based gene–species tree reconciliation framework that defines **branch identity** on a fixed species-tree backbone under two pervasive realities of phylogenomics: **missing taxa** and **gene-tree/species-tree discordance**.

The core idea is to represent each internal branch as a **split (bipartition)**, project species-tree splits onto gene-specific taxon sets, and then perform branch-wise split alignment. This yields standardized **gene×branch mapping tables** and an explicit decomposition of missingness into:

- `NA_struct`: structural missingness (degenerate projected split due to taxon coverage)
- `NA_fuse`: fusion-row missingness (signal represented on composite fused-branch rows under fixed topology)
- `NA_topo`: topology-induced missingness (decisive projected split absent from the free-topology gene tree)

SplitAligner also provides a branch-wise concordance score (**Support**) and accounting identities for consistent bookkeeping.

---

## Preprint

**SplitAligner: A Gene–Species Tree Reconciliation Framework Using Split-Based Branch Mapping**  
bioRxiv (2026) — DOI/link will be added once available.

---

## Repository status

This repository is being populated. The full source code, example commands, and data products used in the preprint will be released here in standard plain-text formats (e.g., Newick, TSV), together with scripts to reproduce figures and key analyses.

Planned contents:
- `src/` SplitAligner core scripts (Perl)
- `examples/` minimal runnable example (species tree, a few gene trees, expected outputs)
- `scripts/` figure reproduction scripts (R)
- `data/` (selected) derived data products used in figures (TSV/Newick)
- `docs/` usage notes and format specifications

---

## Quick overview of outputs

SplitAligner produces:
1. **Fixed-topology mapping table** (branch assignment/branch length on the species-tree backbone; fusion is explicitly represented)
2. **Free-topology mapping table** (captures discordance-driven absence as `NA_topo`)
3. Optional **hybrid (mix-fixed) table** where `NA_topo` entries are imputed using fixed-topology estimates and explicitly flagged
4. **Branch-wise Support annotation** on the species tree

---

## How to cite

If you use SplitAligner, please cite the preprint:

> Wu, J. (2026). SplitAligner: A Gene–Species Tree Reconciliation Framework Using Split-Based Branch Mapping. bioRxiv.

(BibTeX will be added once the bioRxiv DOI is assigned.)

---

## Contact

Jiaqi Wu  
Graduate School of Integrated Sciences for Life, Hiroshima University, Japan  
Email: wujiaqi06@gmail.com, wujiaqi@hiroshima-u.ac.jp

---

## License

Code and documentation: to be finalized (recommended: MIT for code).  
Preprint PDF: CC BY 4.0 (as posted on bioRxiv).
