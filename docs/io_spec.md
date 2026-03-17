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
