# SplitAligner

<p align="center">
  <img src="assets/SplitAligner_logo.png" alt="SplitAligner logo" width="720">
</p>

SplitAligner compares each gene tree with a reference species tree and writes a **gene × branch table** of branch lengths. When a gene lacks species, two reference branches can correspond to the same gene-tree branch. SplitAligner records that length once for the group.

The identity in the poem counts each original reference branch once per gene. Each fusion member counts as `NA_fuse`. The group table records `NA_fuse_topo` if the shared split is absent from the free tree.

**Current prerelease:** [v2.0.0-rc.1](https://github.com/wujiaqi06/SplitAligner/releases/tag/v2.0.0-rc.1).

[Run Fast](#splitalignerfast) · [Perl quick start](#quick-start) · [Use SplitAlignerR](https://github.com/wujiaqi06/SplitAlignerR#installation)

## Introduction

> On a fixed species-tree spine,\
> we ask one thing: does branch *b* still hold?
>
> Project the split:\
> If it collapses — **NA_struct**.\
> If branches fuse — **Bs1|Bs3**, **NA_fuse**.\
> If topology turns away — **NA_topo**.
>
> No ghosts, no leaks:\
> **Total = Mapped + NA_struct + NA_fuse + NA_topo**.\
> A quiet ledger—where every absence has a name.

## SplitAlignerFast

An Apple Silicon executable for macOS 14.4 or later is available in the [release downloads](https://github.com/wujiaqi06/SplitAligner/releases/tag/v2.0.0-rc.1). Other platforms can build from source.

Fast processes one set of gene trees at a time. The following source build works on Linux and macOS and requires CMake 3.16 or later, C/C++17 compilers, OpenSSL development files, and Python 3 for the example. For native Windows, see the [Windows build instructions](fast/docs/WINDOWS.md). Choose a new output directory for each run.

```bash
git clone https://github.com/wujiaqi06/SplitAligner.git
cd SplitAligner
cmake -S fast -B ../splitalignerfast-build -DCMAKE_BUILD_TYPE=Release
cmake --build ../splitalignerfast-build -j2
fast/examples/quickstart.sh ../splitalignerfast-build/saf ../quartet-output
cat ../quartet-output/matrix.tsv
```

```text
gene_id    A0  A1  A2  A3  A4
same        1   2   9   4   5
different   1   4  NA   2   5
```

Column `A2` is the split between {A, B} and {C, D}. The gene called `same` has a matching branch of length `9`. For `different`, that branch is absent and the entry is `NA`. Across both genes the split has 50% support. Each gene-tree run also writes [the support table and annotated reference tree](fast/docs/OUTPUT.md#branch-support). [Run Fast on your own trees](fast/README.md) · [Use the core from C, C++ or Python](fast/docs/LIBRARY.md)

The shared fusion matrix puts each group in one column. In a gene where the group forms, its value is the measured edge length. In other genes the value can be a sum when complete, recovered edges cover exactly those members and each edge has a length. Use `--fusion-view direct` for one normalized observation per value, or `--fusion-view off` to omit the shared fusion matrix. [Fusion values and their source edges](fast/docs/OUTPUT.md)

Historical measurements and their exact workloads are retained in [the benchmark notes](docs/fast-timing.md). They refer to the recorded source versions and output settings.

The same Fast core is available through [SplitAlignerR](https://github.com/wujiaqi06/SplitAlignerR#installation) and the [Python binding](fast/docs/LIBRARY.md).

## Quick Start

The Perl implementation follows the paper’s explicit branch-mapping steps. Version 2 uses compact internal split keys and writes scientific tables plus an indexed JSON audit ledger. Use the [gene audit command](docs/perl-v2-compact.md) to trace a matrix entry back to its input tree. The individual helper scripts remain available. Start v2 in a new output directory; `--force` does not convert a v1 result directory.

To try the small example on macOS or Linux, use Perl 5 and Bash:

```bash
git clone https://github.com/wujiaqi06/SplitAligner.git
cd SplitAligner
bash examples/run.sh toy
head -n 3 examples/302mammal/run/final.free.na_classified.txt | cut -f1-4
```

The preview shows two genes and three reference branches. For gene `A1BG`, branch `B1` has length `0.1809315089`. [Read the output in a six-species example](docs/six-species.md). For your own data, see the [input formats and commands](docs/workflow-reference.md#input-formats).

For native Windows Perl requirements and its test command, see [the Windows guide](docs/windows-candidate.md).

[Project overview](https://wujiaqi06.github.io/#software) · [Tutorial](https://wujiaqi06.github.io/splitaligner/SplitAligner_tutorial.html) · [Animation](https://wujiaqi06.github.io/splitaligner/SplitAligner_explainer.html) · [Method paper (preprint)](https://doi.org/10.64898/2026.02.24.707838) · [Software citation](CITATION.cff)

## Acknowledgments

I thank Chao Zhang (张超; [@chaoszhang](https://github.com/chaoszhang)) of the School of Life Sciences, Peking University, for his technical help and for sharing his `fastSplitAligner.py` implementation.
