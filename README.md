# SplitAligner

**SplitAligner** is a split-based gene tree–species tree reconciliation framework for robust branch mapping under missing taxa, fused branches, and gene-tree/species-tree discordance.

It defines branch identity on a fixed species-tree backbone using canonicalized unrooted edge splits, projects the species-tree split space onto each gene tree according to its observed taxon set, and generates standardized gene-by-branch mapping / branch-length matrices for downstream comparative analyses.

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
- generating branch matrices for comparative analyses,
- separating structural, fusion-related, and topology-induced missingness.

---

## Repository structure

```text
SplitAligner/
  SplitAligner.pl
  README.md
  LICENSE
  CITATION.cff
  .gitignore

  scripts/
    label_species_tree.pl
    tree_to_splits.pl
    split_branch_label.pl
    generate_branch_matrix.pl
    extract_na_fuse.pl
    confirm_na_structure.pl

  examples/
    302mammal/
      input/
        speciesTree302.nwk
        free_tree.examples.nwk
        fix_tree.examples.nwk
      run.sh
      expected/
        free.matrix_with_fuse.txt
        fix.matrix_with_fuse.txt
        final.fix.na_classified.txt
        final.free.na_classified.txt
        species_tree.branch_map.txt

  docs/
    io_spec.md
    algorithm.md
    faq.md
```
## Installation

SplitAligner is implemented in Perl.

### Requirements

- Perl 5
- Standard Perl modules:
  - `Getopt::Long`
  - `File::Basename`
  - `File::Path`

No external R scripts are required for the main workflow.

### Setup

Clone the repository, make the main script executable, and add the project root to your `PATH`:

```bash
git clone https://github.com/wujiaqi06/SplitAligner.git
cd SplitAligner
chmod +x SplitAligner.pl
export PATH="$PWD:$PATH"
```

### Quick start

The main workflow consists of two stages:
	1.	generate branch matrices from a species tree and a set of gene trees,
	2.	finalize NA classification by comparing the fixed-tree and free-tree matrices.

Matrix mode

Generate the branch matrix for free-topology gene trees:
```
SplitAligner.pl --mode matrix \
  --species input/speciesTree302.nwk \
  --gene input/free_tree.examples.nwk \
  --label free
```
Generate the branch matrix for fixed-topology gene trees:
```
SplitAligner.pl --mode matrix \
  --species input/speciesTree302.nwk \
  --gene input/fix_tree.examples.nwk \
  --label fix
```
Finalize mode
Compare the fixed-tree and free-tree outputs and classify NA states:
```
SplitAligner.pl --mode finalize \
  --free free.matrix_with_fuse.txt \
  --fix fix.matrix_with_fuse.txt \
  --final_label final
```
### Example
A runnable example is provided in:
```
examples/302mammal/
```
Run the example as follows:
```
cd examples/302mammal
bash run.sh
```
Key expected outputs are provided in:
```
examples/302mammal/expected/
```
## Command-line interface
--mode matrix

Required arguments:
	•	--species : species tree in Newick format
	•	--gene : gene trees in Newick format
	•	--label : output label prefix

Example:
```
SplitAligner.pl --mode matrix \
  --species input/speciesTree302.nwk \
  --gene input/free_tree.examples.nwk \
  --label free
```
--mode finalize

Required arguments:
	•	--free : matrix generated from free-topology gene trees
	•	--fix : matrix generated from fixed-topology gene trees
	•	--final_label : output prefix for final NA-classified matrices
Example:
```
SplitAligner.pl --mode finalize \
  --free free.matrix_with_fuse.txt \
  --fix fix.matrix_with_fuse.txt \
  --final_label final
```
### Input files

Species tree
	•	Format: Newick
	•	One species tree per run
	•	Species names must be consistent with the names used in the gene trees

Gene trees
	•	Format: Newick
	•	Species names must match the species-tree naming convention
	•	The current workflow assumes gene trees are supplied in the project-specific input format used by the example dataset

Fixed-tree vs free-tree inputs

SplitAligner can be applied to both:
	•	free-topology gene trees, typically inferred independently for each gene
	•	fixed-topology gene trees, typically constrained to a species-tree topology or another fixed-topology framework

The final NA classification step compares these two outputs.

### Output files

Matrix mode outputs

Running --mode matrix generates the following main files:
	•	species_tree.forSplit.nwk
species tree relabeled for split processing
	•	speciesTree.FigTree.tre or equivalent visualization tree output
species tree annotated for visualization
	•	species_tree.splits.txt
species-tree split definitions
	•	species_tree.branch_map.txt
mapping between branch labels and split definitions
	•	<label>_splits/
split representations for gene trees
	•	<label>_split_branch_label/
branch-mapped gene-tree split outputs
	•	<label>.matrix_no_fuse.txt
branch matrix without fused-branch completion
	•	<label>.matrix_with_fuse.txt
branch matrix with fused-branch completion

Finalize mode outputs

Running --mode finalize generates:
	•	<free>.na_fuse.txt
	•	<fix>.na_fuse.txt
	•	<final_label>.fix.na_classified.txt
	•	<final_label>.free.na_classified.txt

The final NA-classified matrices are generated for genes shared between the fixed-tree and free-tree inputs.
If no shared genes are found, the program stops with an error.

⸻

### Interpretation of NA states

SplitAligner distinguishes several biologically meaningful NA classes.
	•	NA
generic missing value before final classification
	•	NA_fuse
the branch is absent as a primitive branch but is represented through a fused branch after taxon pruning
	•	NA_struct
the branch is missing in both fixed-tree and free-tree matrices, consistent with structural absence after projection
	•	NA_topo
the branch is present in the fixed-tree matrix but absent in the free-tree matrix, consistent with topology-induced discordance

⸻

### Workflow overview

The SplitAligner workflow can be summarized as follows:
	1.	Label the species tree.
	2.	Convert species-tree branches into canonicalized unrooted edge splits.
	3.	Convert each gene tree into split form.
	4.	Project the species-tree split space after pruning taxa absent from each gene tree.
	5.	Identify exact and fused branch correspondences.
	6.	Generate gene-by-branch matrices.
	7.	Compare fixed-tree and free-tree matrices to classify NA states.

⸻

### Notes
	•	Branch identity is defined by split representation, not by node order.
	•	Missing taxa can collapse multiple species-tree branches into the same projected split.
	•	Fused branches are explicitly tracked rather than ignored.
	•	Final NA classification is based only on genes shared between the fixed-tree and free-tree inputs.

⸻

### Preprint

SplitAligner: A Gene–Species Tree Reconciliation Framework Using Split-Based Branch Mapping
bioRxiv (2026), Jiaqi Wu.
https://doi.org/10.64898/2026.02.24.707838

⸻

Documentation

Additional documentation is available in:
	•	docs/io_spec.md
	•	docs/algorithm.md
	•	docs/faq.md

⸻

### Citation

If you use SplitAligner in your work, please cite the associated preprint/manuscript.

A CITATION.cff file is provided for GitHub citation support.

⸻

### Contact

Jiaqi Wu
Graduate School of Integrated Sciences for Life, Hiroshima University
Email: wujiaqi@hiroshima-u.ac.jp, wujiaqi06@gmail.com

