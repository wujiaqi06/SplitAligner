<p align="center">
  <img src="assets/SplitAligner_logo.png" alt="SplitAligner logo" width="720">
</p>

# SplitAligner

Missing species and differences between gene trees can make it hard to tell which branches correspond. SplitAligner compares gene trees with a reference species tree and gives you a table of branch lengths, one row per gene. When several reference branches share a single gene-tree length, it records that length once for the group.

[Run the Perl example](#quick-start) · [Run Fast on Linux](#splitalignerfast-on-linux) · [Use R](https://github.com/wujiaqi06/SplitAlignerR#installation)

## Introduction

**SplitAligner**

On a fixed species-tree spine,

we ask one thing: does branch *b* still hold?

Project the split:

If it collapses — **NA_struct**.

If branches fuse — **Bs1|Bs3**, **NA_fuse**.

If topology turns away — **NA_topo**.

No ghosts, no leaks:

**Total = Mapped + NA_struct + NA_fuse + NA_topo + NA_fuse_topo**.

A quiet ledger—where every absence has a name.

## Quick Start

To try the five-gene example on macOS or Linux, use Perl 5 and Bash:

```bash
git clone https://github.com/wujiaqi06/SplitAligner.git
cd SplitAligner
bash examples/run.sh toy
head -n 3 examples/302mammal/run/final.free.na_classified.txt | cut -f1-4
```

The preview shows two genes and three reference branches. For gene `A1BG`, branch `B1` has length `0.1809315089`. For your own data, see the [input formats and commands](docs/workflow-reference.md#input-formats).

## SplitAlignerFast on Linux

To run the Fast command-line program, start in the repository root. You need CMake 3.16 or later, C and C++17 compilers, and the OpenSSL development library. Choose a new output directory for each run.

```bash
cmake -S fast -B ../splitalignerfast-build -DCMAKE_BUILD_TYPE=Release
cmake --build ../splitalignerfast-build -j2
fast/examples/quickstart.sh ../splitalignerfast-build/saf ../quartet-output
cat ../quartet-output/matrix.tsv
```

Look at column `A2` in this example. The gene called `same` has a matching branch of length `9`. For `different`, that branch is absent and the entry is `NA`. [Run Fast on your own trees](fast/README.md)

In our [single-thread Linux benchmark](docs/fast-timing.md), Fast read 2,275 gene trees from the 302-mammal example and wrote the branch matrix in about **4 seconds**, or the full output in **8–9 seconds**.

On macOS, the Fast core also runs through [SplitAlignerR](https://github.com/wujiaqi06/SplitAlignerR#installation) or the [Python binding](fast/docs/LIBRARY.md).

[Project overview](https://wujiaqi06.github.io/#software) · [Tutorial](https://wujiaqi06.github.io/splitaligner/SplitAligner_tutorial.html) · [Animation](https://wujiaqi06.github.io/splitaligner/SplitAligner_explainer.html) · [Method paper (preprint)](https://doi.org/10.64898/2026.02.24.707838) · [Software citation](CITATION.cff)

## Acknowledgments

I thank Chao Zhang (张超) for his technical help and for sharing his `fastSplitAligner.py` implementation.
