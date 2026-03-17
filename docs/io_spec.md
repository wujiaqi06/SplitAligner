# SplitAligner I/O specification

## Overview

SplitAligner operates in two major modes:

1. `matrix` mode: generate branch-mapping matrices from a species tree and a set of gene trees
2. `finalize` mode: compare fixed-topology and free-topology matrices to classify NA states

---

## Input files

### 1. Species tree

- Format: Newick
- One species tree per run
- Species labels must be consistent with those used in the gene trees

Example:
```text
((A,B),(C,D));
```
### 2. Gene trees
	•	Format: Newick
	•	Species labels must match the species-tree naming convention
	•	The current implementation assumes the project-specific input format used in the provided example dataset

### Two major use cases are supported:
	•	free-topology gene trees: inferred independently for each gene
	•	fixed-topology gene trees: inferred under a fixed-topology framework

⸻

### Matrix mode

Command
```
SplitAligner.pl --mode matrix \
  --species input/speciesTree302.nwk \
  --gene input/free_tree.examples.nwk \
  --label free
```
Required arguments
	•	--species : species tree file
	•	--gene : gene tree file
	•	--label : output prefix

Main outputs
	•	species_tree.forSplit.nwk
	•	speciesTree.FigTree.tre or equivalent visualization tree output
	•	species_tree.splits.txt
	•	species_tree.branch_map.txt
	•	<label>_splits/
	•	<label>_split_branch_label/
	•	<label>.matrix_no_fuse.txt
	•	<label>.matrix_with_fuse.txt

⸻

### Finalize mode

Command
```
SplitAligner.pl --mode finalize \
  --free free.matrix_with_fuse.txt \
  --fix fix.matrix_with_fuse.txt \
  --final_label final
```
Required arguments
	•	--free : free-topology matrix
	•	--fix : fixed-topology matrix
	•	--final_label : output prefix

Main outputs
	•	<free>.na_fuse.txt
	•	<fix>.na_fuse.txt
	•	<final_label>.fix.na_classified.txt
	•	<final_label>.free.na_classified.txt

The final NA-classified matrices are generated for genes shared between the fixed-tree and free-tree inputs.
If no shared genes are found, the program stops with an error.

⸻

### NA states

NA

Generic missing value before final classification.

NA_fuse

The branch is absent as a primitive branch but represented through a fused branch after taxon pruning.

NA_struct

The branch is missing in both fixed-tree and free-tree matrices, consistent with structural absence after projection.

NA_topo

The branch is present in the fixed-tree matrix but absent in the free-tree matrix, consistent with topology-induced discordance.

⸻

### Notes
	•	Branch identity is defined by split representation, not by node order.
	•	Missing taxa can collapse multiple species-tree branches into the same projected split.
	•	Fused branches are explicitly tracked rather than ignored.
---
