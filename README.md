<p align="center">
  <img src="assets/SplitAligner_logo.png" alt="SplitAligner logo" width="720">
</p>

# SplitAligner

Which branch in a gene tree corresponds to the branch you want to study in the species tree? Give SplitAligner a reference species tree and a set of gene trees. It returns a table with one row per gene and one column per reference branch. A cell gives the matching branch length when available. A separate table holds lengths that belong to several reference branches together.

[Watch the animation](https://wujiaqi06.github.io/splitaligner/SplitAligner_explainer.html) · [Follow the tutorial](https://wujiaqi06.github.io/splitaligner/SplitAligner_tutorial.html) · [Run the Perl example](#quick-start) · [Use Fast on Linux](#splitalignerfast-on-linux) · [Install the R package](https://github.com/wujiaqi06/SplitAlignerR#installation)

## Introduction

**SplitAligner**

On a fixed species-tree spine,

we ask one thing: does branch *b* still hold?

Project the split: **σ_g(b)**.

If it collapses — **NA_struct**.

If branches fuse — **Bs1|Bs3**, **NA_fuse**.

If topology turns away — **NA_topo**.

No ghosts, no leaks:

**Total = Mapped + NA_struct + NA_fuse + NA_topo + NA_fuse_topo**.

A quiet ledger—where every absence has a name.

## What happens when a species is missing?

Cut branch *b* in the reference tree and read the species on either side. Those two lists are its *split*. For gene *g*, cross out every species missing from that gene. The shortened lists are the projected split, **σ_g(b)**.

Suppose branch *x* separates **A B | C D E F**, and branch *y* separates **A B C | D E F**. When C is missing, both read **A B | D E F**. Look for that split on the gene tree:

| What the gene tree shows | Cell for x | Cell for y | Result for the group {x, y} |
| --- | --- | --- | --- |
| A matching branch of length 0.31 | `NA_fuse` | `NA_fuse` | 0.31, recorded once |
| No branch with that split | `NA_fuse` | `NA_fuse` | `NA_fuse_topo` |

The tip branch for C has an empty side after C is removed. Its cell reads `NA_struct`. To interpret a branch that survives on its own, the Perl workflow compares two trees for the gene. One follows the reference topology (the *fixed* tree). The other is inferred freely (the *free* tree). If the fixed tree gives the branch a finite length and the free tree lacks its split, the branch reads `NA_topo`.

The poem’s **Total** counts outcomes across the original reference branches. For the recovery rate of an individual branch, a member of an `NA_fuse_topo` group contributes zero to both numerator and denominator. [Counting units and recovery rules](docs/workflow-reference.md#definition-summary)

## Quick Start

The Perl example uses five mammal genes and requires Perl 5 with its core modules. [Download or clone the repository](docs/workflow-reference.md#installation), then run this from its top directory:

```bash
bash examples/run.sh toy
```

To see the first three reference branches for the first two genes:

```bash
head -n 3 examples/302mammal/run/final.free.na_classified.txt | cut -f1-4
```

For gene **A1BG**, B1 has length **0.1809315089** and B2 reads **`NA_struct`**. Open `final.free.projected_groups.tsv` in the same directory for the fused groups and their lengths. Use the branch map from the same run to look up B1. The labels can change if you reorder the tips in the reference tree. [Why B labels are serialization-local aliases](docs/workflow-reference.md#determinism-and-coordinate-identity)

[Use your own trees](docs/workflow-reference.md#input-formats) · [Run the full 2,275-gene example](docs/workflow-reference.md#running-the-examples)

## SplitAlignerFast on Linux

Use the `saf` command to map one collection of gene trees. Build it from the repository root with CMake 3.16 or later, C and C++17 compilers, and the OpenSSL development library. Python 3.8 or later enables the Python examples and tests.

```bash
cmake -S fast -B ../splitalignerfast-build -DCMAKE_BUILD_TYPE=Release
cmake --build ../splitalignerfast-build -j2
(cd ../splitalignerfast-build && ctest --output-on-failure)
```

The bundled example has two four-species trees, `same` and `different`. From the repository root, choose an output directory that does not yet exist:

```bash
../splitalignerfast-build/saf \
  --reference fast/examples/reference.nwk \
  --genes fast/examples/genes.tsv \
  --scope full --role free --out ../quartet-output

cat ../quartet-output/matrix.tsv
cat ../quartet-output/raw_states.tsv
```

In `raw_states.tsv`, find the column that changes from **D** to **U**. This is the internal branch. The `same` tree has this split and its length in `matrix.tsv` is **9**. The `different` tree lacks the split and its length entry is **NA**. With your own data, replace the reference file and supply a headerless file with one gene ID, a tab and a Newick tree per line.

Full output includes `composites.jsonl`, with one record per fused group. Add `--unmatched-splits off` if you want to omit the two unmatched-split files. If you downloaded the standalone Fast source archive, start with the README in the extracted directory. [Fast installation](fast/README.md) · [Output guide](fast/docs/OUTPUT.md) · [Reuse a mapping in C or Python](fast/docs/LIBRARY.md)

## Documentation

- [Perl workflow reference](docs/workflow-reference.md) — fixed and free trees, command options, file formats, and interpretation of the output.
- [Choosing between Perl, Fast and R](docs/implementations.md) — including how to reuse a mapping across branch-length draws.
- [Algorithm](docs/algorithm.md) · [Benchmark rules](docs/benchmark_rules.md) · [FAQ](docs/faq.md)

## Citation

Jiaqi Wu (2026). *SplitAligner: A Gene-Species Tree Reconciliation Framework Using Split-Based Branch Mapping.* [bioRxiv preprint](https://doi.org/10.64898/2026.02.24.707838) · [Software citation](CITATION.cff)

## Acknowledgments

We thank Chao Zhang (张超) for his technical help and for sharing his `fastSplitAligner.py` implementation.

## Contact

Jiaqi Wu, Graduate School of Integrated Sciences for Life, Hiroshima University. [wujiaqi@hiroshima-u.ac.jp](mailto:wujiaqi@hiroshima-u.ac.jp) · [wujiaqi06@gmail.com](mailto:wujiaqi06@gmail.com)
