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
  SplitAligner.pl          # main controller
  README.md
  LICENSE
  CITATION.cff
  .gitignore

  scripts/                 # internal pipeline steps
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

⸻

## Installation

SplitAligner is implemented in Perl.

## Requirements
	•	Perl 5
	•	Standard Perl modules (core):
	•	Getopt::Long
	•	File::Basename
	•	File::Path
	•	FindBin
	•	File::Spec

## Setup

Clone the repository and run via perl:

git clone https://github.com/wujiaqi06/SplitAligner.git
cd SplitAligner

(Optionally) make the controller executable:

chmod +x SplitAligner.pl


⸻

## Quick start

A runnable example is provided in:

examples/302mammal/

Run the example as follows:

cd examples/302mammal
bash run.sh

Key expected outputs are provided in:

examples/302mammal/expected/


⸻

## Command-line interface

### --mode matrix

Generate branch matrices from a species tree and a single gene-tree set.

#### Required arguments
	•	--species: species tree in Newick format
	•	--gene: gene trees in the SplitAligner line-based input format (one record per line)

#### Optional arguments
	•	--label: output label prefix (default: inferred from the gene-tree filename)
	•	--axis: species-axis label prefix (default: species_tree; useful for parallel runs)

#### Example (free-topology gene trees)
```text
perl ../../SplitAligner.pl --mode matrix \
  --species input/speciesTree302.nwk \
  --gene input/free_tree.examples.nwk \
  --label free
```
#### Example (fixed-topology gene trees)
```text
perl ../../SplitAligner.pl --mode matrix \
  --species input/speciesTree302.nwk \
  --gene input/fix_tree.examples.nwk \
  --label fix
```

⸻

### --mode finalize

Finalize NA classification by comparing two with-fuse matrices (typically from fixed-tree vs free-tree runs).

```text
[Finalize mode: input requirements]
Please provide TWO "with-fuse" branch matrices (the large matrices):
  --free <FREE matrix file>   (header should include fused branches, e.g., B101|B102; often named *matrix_with_fuse.txt)
  --fix  <FIX  matrix file>   (header should include fused branches, e.g., B101|B102; often named *matrix_with_fuse.txt)
  --final_label <output prefix>

Finalize will run:
  1) extract_na_fuse.pl -i <FREE matrix file>   (produces <FREE>.na_fuse.txt)
  2) extract_na_fuse.pl -i <FIX  matrix file>   (produces <FIX>.na_fuse.txt)
  3) confirm_na_structure.pl -fix <FIX>.na_fuse.txt -free <FREE>.na_fuse.txt -o <final_label>

All executed commands are recorded in: commands.txt
```
#### Example

perl ../../SplitAligner.pl --mode finalize \
  --free free.matrix_with_fuse.txt \
  --fix fix.matrix_with_fuse.txt \
  --final_label final


⸻

### Input files

#### Species tree
	•	Format: Newick
	•	One species tree per run
	•	Species names must be consistent with the names used in the gene trees
	•	The following species-tree formats are accepted:
```text
((A,B),(C,D));
```
```text
((A:0.1,B:0.2):0.2,(C:0.1,D:0.1):0.1):0.1;
```
```text
((A:0.1,B:0.2)100:0.2,(C:0.1,D:0.1)95:0.1)100:0.1;
```
These correspond to:
	1.	topology only,
	2.	topology with branch lengths,
	3.	topology with branch lengths and internal node annotations (e.g., support values).

Internal tree annotations are removed during preprocessing and are not used in downstream calculations.

Gene trees
	•	Format: Newick (line-based records; one record per line)
	•	Species names must match the species-tree naming convention
	•	The current workflow assumes the project-specific input format used in the example dataset:
	•	each line begins with a gene ID followed immediately by a Newick tree
  •	Gene trees may include node support values (e.g., bootstrap). SplitAligner does not use these values for branch mapping or matrix construction.

Examples:
```text
GeneA((A:0.1,B:0.2):0.2,(C:0.1,D:0.1):0.1):0.1;
GeneB((A:0.2,B:0.1):0.1,(C:0.1,(D:0.1,E:0.2):0.1):0.1):0.1;
GeneC((A:0.1,B:0.2):0.2,(C:0.1,D:0.1):0.4):0.1;
```
Internal node annotations are allowed but are not used in split-based branch mapping.

⸻

### Output files

#### Matrix mode outputs

Running --mode matrix generates:
	•	speciesTree.forSplit.nwk
Species tree relabeled for split processing.
	•	speciesTree.FigTree.tre
Species tree annotated for visualization (FigTree-ready).
	•	<axis>.splits.txt
Species-tree split definitions (default <axis>=species_tree).
	•	<axis>.branch_map.txt
Mapping between branch labels and split definitions.
	•	<label>_splits/
Split representations for gene trees.
	•	<label>_split_branch_label/
Branch-mapped gene-tree split outputs.
	•	<label>.matrix_no_fuse.txt
Branch matrix without fused-branch augmentation.
	•	<label>.matrix_with_fuse.txt
Branch matrix including fused-branch columns (when present).

In addition, SplitAligner writes:
	•	commands.txt (exact executed commands, for provenance/debug)
	•	run.log (timestamps and runtimes)

#### Finalize mode outputs

Running --mode finalize generates:
	•	<FREE_MATRIX>.na_fuse.txt
e.g., free.matrix_with_fuse.txt.na_fuse.txt
	•	<FIX_MATRIX>.na_fuse.txt
e.g., fix.matrix_with_fuse.txt.na_fuse.txt
	•	<final_label>.fix.na_classified.txt
	•	<final_label>.free.na_classified.txt

The final NA-classified matrices are generated for genes shared between the fixed-tree and free-tree inputs.
If no shared genes are found, the program stops with an error.

⸻

### Interpretation of NA states

SplitAligner distinguishes several biologically meaningful NA classes:
	•	NA
Generic missing value before final classification.
	•	NA_fuse
The branch is absent as a primitive branch but is represented through a fused branch after taxon pruning.
	•	NA_struct
The branch is missing in both fixed-tree and free-tree matrices, consistent with structural absence after projection.
	•	NA_topo
The branch is present in the fixed-tree matrix but absent in the free-tree matrix, consistent with topology-induced discordance.

⸻

### Workflow overview
	1.	Label the species tree.
	2.	Convert species-tree branches into canonicalized unrooted edge splits.
	3.	Convert each gene tree into split form.
	4.	Project the species-tree split space after pruning taxa absent from each gene tree.
	5.	Identify exact and fused branch correspondences.
	6.	Generate gene-by-branch matrices.
	7.	Compare fixed-tree and free-tree matrices to classify NA states.

⸻

### Preprint

SplitAligner: A Gene–Species Tree Reconciliation Framework Using Split-Based Branch Mapping
bioRxiv (2026), Jiaqi Wu.
DOI: 10.64898/2026.02.24.707838

⸻

### Documentation

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
Email: wujiaqi@hiroshima-u.ac.jp and wujiaqi06@gmail.com
