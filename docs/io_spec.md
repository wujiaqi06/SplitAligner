# SplitAligner I/O Specification

## Overview

SplitAligner operates in two major modes:

1. `matrix`: generate branch-mapping matrices from a species tree and one gene-tree file
2. `finalize`: compare fixed-topology and free-topology matrices to classify NA states

This document summarizes the expected inputs and the main outputs of each mode.

---

## Input Files

### 1. Species tree

- Format: Newick
- One species tree per run
- Species labels must match those used in the gene trees
- Branch lengths and internal node annotations are allowed
- Internal annotations are ignored during split-based mapping

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

Main outputs:

- `species_tree.forSplit.nwk`
- `species_tree.FigTree.tre`
- `species_tree.splits.txt`
- `species_tree.branch_map.txt`
- `<label>_splits/`
- `<label>_split_branch_label/`
- `<label>.matrix_no_fuse.txt`
- `<label>.matrix_with_fuse.txt`

Interpretation:

- `species_tree.splits.txt` defines the canonical branch coordinate system on the species-tree backbone
- `<label>_splits/` stores per-gene split representations
- `<label>_split_branch_label/` stores per-gene mappings onto the species-tree branch coordinate system
- `<label>.matrix_no_fuse.txt` contains primitive branches only
- `<label>.matrix_with_fuse.txt` additionally includes fused-branch columns such as `B12|B47`

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

Main outputs:

- `<free>.na_fuse.txt`
- `<fix>.na_fuse.txt`
- `<final_label>.fix.na_classified.txt`
- `<final_label>.free.na_classified.txt`
- `<final_label>.support_b.txt` if `--species_tree` is provided
- `<species_prefix>.support_b.nwk` if `--species_tree` is provided

Interpretation:

- `extract_na_fuse.pl` marks primitive-branch `NA` cells that are explained by non-NA fused-branch signal as `NA_fuse`
- `confirm_na_structure.pl` compares fixed-topology and free-topology matrices on shared genes
- the final classified outputs distinguish structural missingness from topology-induced missingness
- if `--species_tree` is provided, `finalize` also computes branch-wise `Support` and writes an annotated species tree

The final comparison is defined only for genes shared between the fixed-topology and free-topology inputs. If no shared genes are found, SplitAligner stops with an error.

### Branch-wise `Support`

When `--species_tree` is provided, SplitAligner computes `Support(b)` for each branch on the species-tree backbone.

- denominator: number of non-NA entries for branch `b` in the fixed-topology matrix on shared genes
- numerator: number of non-NA entries for branch `b` in the free-topology matrix on shared genes

The resulting outputs are:

- `<final_label>.support_b.txt`: branch-wise summary table
- `<species_prefix>.support_b.nwk`: standard Newick tree with internal-node `Support` values written in the bootstrap position

---

## NA States

### `NA`

Generic missing value before final classification.

### `NA_fuse`

The branch is absent as a primitive branch but represented through a fused branch after taxon pruning.

### `NA_struct`

The projected branch is structurally absent after projection and is not evaluable for that gene.

### `NA_topo`

The projected branch is decisive under the fixed-topology comparison but absent from the free-topology gene tree, consistent with topology-induced discordance.

---

## Notes

- Branch identity is defined by split representation rather than by node order or drawing position.
- Missing taxa can collapse multiple species-tree branches into the same projected split.
- Fused branches are explicitly tracked rather than ignored.
- The `finalize` step is intended to separate coverage-driven missingness from discordance-driven missingness.
