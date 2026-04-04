# Benchmark Final

This folder packages the current Benchmark V1 working bundle for SplitAligner.

## Core Principle

This benchmark is only meaningful if the R-side oracle stays independent of SplitAligner.

Its job is to infer `NA_struct` and `NA_fuse` explicitly from tip deletion, degree-2 contraction, and node/edge membership tracking on the frozen full-tree axis. It must not re-implement SplitAligner in R.

Because `NA_struct` and `NA_fuse` are defined here as consequences of pruning and contraction on the original full-tree primitive-edge axis, they are recoverable from graph-only node/edge incidence and contraction history alone. Reintroducing split logic would no longer validate SplitAligner independently, but would merely restate it in another encoding.

- `NA_struct` and `NA_fuse` are inferred from explicit deletion-to-edge mapping and contraction history.
- The R benchmark must never use split-based missingness logic.
- The R benchmark must never compute projected splits, split keys, or split-space equivalence classes.
- Split logic belongs only to the Perl SplitAligner workflow.

## Development Note

Why the benchmark oracle must remain structurally independent:

1. Initial implementations repeatedly drifted back toward split-based logic, which would have compromised oracle independence from SplitAligner.
2. Once a split-free R-side oracle was enforced, it became clear that tip deletion induces structural collapse and degree-2 contraction, causing `NA_fuse` spillover beyond the immediately deleted local region.
3. Purely local, recursive, or offspring-node-based heuristics were not sufficient:
   - recursive propagation overcalled entire clades,
   - local offspring tracking missed higher-level spillover,
   - global deletion without explicit edge tracing still missed fused-edge cases.
4. The final benchmark therefore uses only tree surgery, degree-2 contraction, and explicit node/edge tracking to infer `NA_struct` and `NA_fuse`, without relying on split-based representations.

## Purpose

This benchmark evaluates `NA_struct` and `NA_fuse` under fixed-topology pruning on a frozen full-tree primitive-edge axis.

In plainer terms, all benchmark labels are defined relative to the branch identities of the original full tree, before any pruning-induced collapse.

It is an oracle-validated benchmark bundle, not a generic phylogeny-accuracy benchmark.

## Hard Rule

The R-side oracle must remain structurally independent of SplitAligner.

- R must not use split-based logic to classify missingness.
- R must not compute projected splits or split-space equivalence classes.
- Split-based branch mapping remains exclusive to the Perl SplitAligner workflow.
- The R oracle uses only graph-based tree surgery, degree-2 contraction, and explicit node/edge tracking.

## Directory Layout

- `inputs/`
  - example tree input used for the packaged benchmark run
- `scripts/`
  - `benchmark.R`: graph-only benchmark driver
  - `oracle_utils.R`: graph-only oracle utilities
  - `plot_fulltree_collapse.R`: cumulative collapse plot on the full tree
- `outputs/`
  - `unrooted/`
  - `rooted/`
  - benchmark tables, tree files, and plotted PDFs separated by tree semantics
- `docs/`
  - benchmark specification and tasklist

## Main Inputs and Outputs

Primary generated files are written into one of:

- `outputs/unrooted/`
- `outputs/rooted/`

Each semantics-specific output directory contains:

- `benchmark.species_tree.nwk`
- `benchmark.gene_trees.nwk`
- `benchmark.table.txt`
- `benchmark.oracle_events.tsv`
- `benchmark.trajectory.tsv`
- `benchmark.fulltree_collapse.pdf`

## Run

Edit `scripts/benchmark.R` to choose `TREE_SEMANTICS <- "unrooted"` or `TREE_SEMANTICS <- "rooted"`. Benchmark V1 defaults to `unrooted` to match SplitAligner. By default, the packaged run therefore writes into `outputs/unrooted/`.

From this folder:

```bash
Rscript scripts/benchmark.R
```

This writes the benchmark tables into `outputs/unrooted/` or `outputs/rooted/` and automatically generates a matching collapse PDF in the same directory.

If needed, the plotting step can also be rerun by itself:

```bash
Rscript scripts/plot_fulltree_collapse.R
```

## Notes

- The packaged default is `MAX_DELETE_FRACTION = 0.7` because it produces a clearer collapse trajectory for visualization and a more expressive benchmark figure.
- `species_tree.nwk` and `gene_trees.nwk` are formatted to help align the benchmark with the Perl SplitAligner workflow.
- The benchmark tables are intended for direct sanity checking: any numeric cell that has already absorbed contraction-induced branch-length summation should instead have been marked `NA_fuse`.
- The rooted and unrooted modes are intentionally different. In `unrooted` mode, the benchmark now performs explicit pseudo-root normalization and unrooted degree-2 contraction, so branches that become indistinguishable only under the unrooted interpretation are merged exactly as in SplitAligner. In `rooted` mode, the benchmark retains rooted branch distinctions and therefore serves as a separate, biologically interpretable companion analysis rather than as the primary SplitAligner comparison target.
- Historical experiments and failed drafts are intentionally excluded from this packaged folder.
