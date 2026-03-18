# SplitAligner

**SplitAligner** is a split-based gene tree-species tree reconciliation framework for robust branch mapping under missing taxa, fused branches, and gene-tree/species-tree discordance.

It defines branch identity on a fixed species-tree backbone using canonicalized unrooted edge splits, projects that split space onto each gene tree according to the taxa observed in that gene, and generates standardized gene-by-branch matrices for downstream comparative analyses.

SplitAligner explicitly distinguishes biologically meaningful forms of missingness:

- `NA_struct`: structural absence after taxon pruning causes a projected branch to become undefined
- `NA_fuse`: signal is represented on a fused branch rather than on a primitive branch
- `NA_topo`: the projected branch is supported in the fixed-topology analysis but absent in the free-topology gene tree

This repository contains the SplitAligner source code, example datasets, and documentation needed to reproduce the core branch-mapping workflow.

---

## Introduction

On a fixed species-tree spine,

we ask one thing: does branch `b` still hold?

Project the split:

If it collapses - `NA_struct`.

If branches fuse - `Bs1|Bs3`, `NA_fuse`.

If topology turns away - `NA_topo`.

No ghosts, no leaks:

`Total = Mapped + NA_struct + NA_fuse + NA_topo.`

A quiet ledger, where every absence has a name.

---

## Why SplitAligner?

In phylogenomics, branch identity is often treated as though it remains stable across all gene trees. In practice, missing taxa can collapse distinct species-tree branches into the same projected split, and free-topology gene trees can lose decisive projected branches entirely. Under these conditions, naive branch-to-branch comparison becomes unreliable.

SplitAligner addresses this problem by:

- defining branch identity in split space rather than by graphical position or node order
- projecting the species-tree split space onto the observed taxon set of each gene
- recording both exact and fused branch correspondences
- generating standardized branch matrices for downstream comparative analyses
- separating structural, fusion-related, and topology-induced missingness

The central rule is simple:

> branch identity should be defined in projected split space, not assumed to survive taxon pruning unchanged.

---

## Conceptual Positioning

SplitAligner reframes branch identity as a projection problem in split space.

Instead of asking whether a gene tree simply "supports" a species-tree branch, SplitAligner asks whether that branch remains well-defined after taxon pruning, whether it becomes structurally degenerate, whether its signal is absorbed into a fused branch, or whether it is absent because of topological discordance.

This shift turns branch reconciliation from a naive tree-comparison problem into a standardized gene-by-branch matrix construction framework under controlled missingness.

In that sense, SplitAligner is not only a branch-mapping tool. It is also a branch-coordinate infrastructure for downstream comparative analyses, where each gene-branch cell can be interpreted within an explicit missingness model rather than as an undifferentiated absence.

---

## Main Features

- Split-based branch mapping on a fixed species-tree backbone
- Explicit handling of missing taxa during per-gene projection
- Recognition of fused branch patterns after taxon pruning
- Matrix generation for fixed-topology and free-topology gene-tree sets
- Final classification of `NA`, `NA_fuse`, `NA_struct`, and `NA_topo`
- Reproducible example workflow included in `examples/302mammal/`

---

## Repository Structure

```text
SplitAligner/
  SplitAligner.pl          main controller
  README.md
  LICENSE
  CITATION.cff

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
      expected/
      run.sh

  docs/
    algorithm.md
    io_spec.md
    faq.md
```

---

## Requirements

- Perl 5
- Core Perl modules:
  - `Getopt::Long`
  - `Getopt::Std`
  - `File::Basename`
  - `File::Path`
  - `File::Spec`
  - `FindBin`
  - `Cwd`

No external R or non-core Perl dependencies are required for the main workflow.

---

## Installation

Clone the repository:

```bash
git clone https://github.com/wujiaqi06/SplitAligner.git
cd SplitAligner
```

Optionally make the main controller executable:

```bash
chmod +x SplitAligner.pl
```

You can run SplitAligner directly from the repository root:

```bash
perl SplitAligner.pl --help
```

Or add the repository root to your `PATH`:

```bash
export PATH="$PWD:$PATH"
SplitAligner.pl --help
```

To make that persistent:

macOS (`zsh`)

```bash
echo 'export PATH="'"$PWD"'":$PATH' >> ~/.zshrc
source ~/.zshrc
```

Linux (`bash`)

```bash
echo 'export PATH="'"$PWD"'":$PATH' >> ~/.bashrc
source ~/.bashrc
```

---

## Quick Start

A runnable example is provided in `examples/302mammal/`.

From the repository root:

```bash
cd examples/302mammal
bash run.sh
```

This example performs:

1. matrix generation for free-topology gene trees
2. matrix generation for fixed-topology gene trees
3. final NA classification by comparing the two matrix sets

Expected reference outputs are provided in `examples/302mammal/expected/`.

---

## Workflow Overview

SplitAligner runs in two major stages.

### Stage 1: `matrix`

1. Label the species tree with stable branch identifiers
2. Convert species-tree branches into canonicalized unrooted edge splits
3. Convert each gene tree into split form
4. Project the species-tree split space after pruning taxa absent from each gene
5. Detect exact and fused branch correspondences
6. Generate gene-by-branch matrices

### Stage 2: `finalize`

1. Mark primitive-branch `NA` cells that are explained by fused-branch signal as `NA_fuse`
2. Compare fixed-topology and free-topology matrices on shared genes
3. Classify remaining missing values as `NA_struct` or `NA_topo`

---

## Command-Line Interface

### `--mode matrix`

Generate branch matrices from a species tree and one gene-tree file.

Required arguments:

- `--species`: species tree in Newick format
- `--gene`: gene-tree file in SplitAligner line-based format
- `--label`: output label or prefix, for example `free` or `fix`

Example: free-topology gene trees

```bash
perl SplitAligner.pl --mode matrix \
  --species input/speciesTree302.nwk \
  --gene input/free_tree.examples.nwk \
  --label free
```

Example: fixed-topology gene trees

```bash
perl SplitAligner.pl --mode matrix \
  --species input/speciesTree302.nwk \
  --gene input/fix_tree.examples.nwk \
  --label fix
```

Main outputs:

- `species_tree.forSplit.nwk`
- `species_tree.FigTree.tre`
- `species_tree.splits.txt`
- `species_tree.branch_map.txt`
- `<label>_splits/`
- `<label>_split_branch_label/`
- `<label>.matrix_no_fuse.txt`
- `<label>.matrix_with_fuse.txt`

### `--mode finalize`

Finalize NA classification from two `matrix_with_fuse` outputs, typically one from fixed-topology gene trees and one from free-topology gene trees.

Required arguments:

- `--free`: free-topology `matrix_with_fuse` file
- `--fix`: fixed-topology `matrix_with_fuse` file
- `--final_label`: output prefix for the classified matrices

Example:

```bash
perl SplitAligner.pl --mode finalize \
  --free free.matrix_with_fuse.txt \
  --fix fix.matrix_with_fuse.txt \
  --final_label final
```

Main outputs:

- `<free>.na_fuse.txt`
- `<fix>.na_fuse.txt`
- `<final_label>.fix.na_classified.txt`
- `<final_label>.free.na_classified.txt`

The final classification step is defined only for genes shared between the fixed-topology and free-topology inputs. If no shared genes are found, SplitAligner stops with an error.

---

## Input Formats

### Species tree

- Format: Newick
- One species tree per run
- Species labels must be consistent with those used in the gene trees
- Branch lengths and internal node annotations are allowed
- Internal annotations are ignored during split-based branch mapping

Accepted examples:

```text
((A,B),(C,D));
```

```text
((A:0.1,B:0.2):0.2,(C:0.1,D:0.1):0.1):0.1;
```

```text
((A:0.1,B:0.2)100:0.2,(C:0.1,D:0.1)95:0.1)100:0.1;
```

### Gene trees

- Format: Newick, one record per line
- Each line begins with a gene identifier followed immediately by a tree
- Species labels must match the species-tree naming convention
- Branch lengths and node support annotations are allowed
- Internal annotations are ignored during split-based mapping

Example:

```text
GeneA((A:0.1,B:0.2):0.2,(C:0.1,D:0.1):0.1):0.1;
GeneB((A:0.2,B:0.1):0.1,(C:0.1,(D:0.1,E:0.2):0.1):0.1):0.1;
GeneC((A:0.1,B:0.2):0.2,(C:0.1,D:0.1):0.4):0.1;
```

---

## Output Files

### Outputs from `matrix`

- `species_tree.forSplit.nwk`
  - species tree relabeled for downstream split processing
- `species_tree.FigTree.tre`
  - species tree annotated for visualization
- `species_tree.splits.txt`
  - canonical species-tree split definitions
- `species_tree.branch_map.txt`
  - mapping between branch identifiers and species-tree subtrees
- `<label>_splits/`
  - per-gene split representations
- `<label>_split_branch_label/`
  - per-gene mapped branch patterns after projection to the species-tree axis
- `<label>.matrix_no_fuse.txt`
  - primitive-branch matrix only
- `<label>.matrix_with_fuse.txt`
  - primitive branches plus fused-branch columns

### Outputs from `finalize`

- `<free>.na_fuse.txt`
  - primitive-branch matrix in which fused-supported `NA` cells are relabeled as `NA_fuse`
- `<fix>.na_fuse.txt`
  - same transformation for the fixed-topology matrix
- `<final_label>.fix.na_classified.txt`
  - fixed-topology matrix after final NA classification
- `<final_label>.free.na_classified.txt`
  - free-topology matrix after final NA classification

---

## Interpretation of NA States

- `NA`
  - generic missing value before final classification
- `NA_fuse`
  - the branch is absent as a primitive branch but represented through a fused branch after taxon pruning
- `NA_struct`
  - the projected branch is structurally absent in both fixed-topology and free-topology comparisons
- `NA_topo`
  - the branch is present in the fixed-topology matrix but absent in the free-topology matrix, consistent with topology-induced discordance

These categories are intended to prevent biologically distinct sources of missingness from being conflated in downstream analyses.

---

## Documentation

Additional documentation is available in:

- `docs/algorithm.md`
- `docs/io_spec.md`
- `docs/faq.md`

---

## Citation

If you use SplitAligner in your work, please cite the software repository and the associated preprint.

Preferred citation:

> Wu J. 2026. SplitAligner: A Gene-Species Tree Reconciliation Framework Using Split-Based Branch Mapping. bioRxiv. https://doi.org/10.64898/2026.02.24.707838

GitHub citation metadata is provided in `CITATION.cff`.

---

## Contact

Jiaqi Wu  
Graduate School of Integrated Sciences for Life, Hiroshima University  
Email: `wujiaqi@hiroshima-u.ac.jp`  
Email: `wujiaqi06@gmail.com`
