# Benchmark Rules

This document records hard constraints for future benchmark implementations associated with SplitAligner.

These rules are methodological boundaries, not optional coding preferences.

## Core Principle

The split-based definition of branch identity and missingness classification belongs to the Perl implementation of SplitAligner.

Any future benchmark implementation in R must not reproduce, approximate, or partially reuse split-based logic to classify missingness states.

## Hard Rules

1. R benchmark code must not use splits to classify missingness.
2. R benchmark code must not infer `NA_struct`, `NA_fuse`, or related missingness states through split-space projection.
3. All split-based branch mapping and split-based missingness classification must remain exclusive to the Perl SplitAligner workflow.
4. If an R benchmark is implemented, it must use an explicit non-split strategy based on `ape` node correspondence or equivalent node-based tree matching.
5. Under that R benchmark design, `NA_fuse` and `NA_struct` must be computed explicitly from node-matching logic rather than by any split-derived shortcut.
6. R code must not use a hybrid or disguised split-based procedure that effectively reproduces SplitAligner while presenting itself as a separate benchmark.

## Practical Interpretation

For benchmark comparisons:

- The Perl implementation is the only allowed source of split-based branch mapping.
- An R implementation may be developed for comparison, but it must remain a genuinely node-based alternative.
- The benchmark must not blur the distinction between SplitAligner's projected split-space framework and an external node-based comparator.

## Audit Target

Packaged benchmark audits compare `splitaligner_perl` strictly against `benchmark_unrooted`.

`benchmark_rooted` outputs are retained as reference outputs for rooted-label conventions and are expected to differ in root-adjacent cases. Such rooted/unrooted differences are not benchmark failures.

A benchmark scenario passes only when `splitaligner_perl` matches `benchmark_unrooted` with zero unexpected mismatches.

## Why This Boundary Exists

The goal of the benchmark is to compare distinct methodological strategies fairly.

If the R benchmark uses split-based logic, even partially, then it is no longer an independent comparator. It becomes a partial reimplementation of SplitAligner, which would invalidate the benchmark design and weaken the interpretability of the comparison.

## Scope

These rules apply to future benchmark code, scripts, helper functions, and downstream wrappers associated with this repository.

They do not restrict future R utilities for reading SplitAligner outputs, plotting results, or summarizing finalized matrices after the Perl workflow has already produced them.
