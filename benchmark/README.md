# Benchmark

This directory contains the SplitAligner benchmark bundle and audit scaffold.

## Scope

The current packaged benchmark bundle contains two deterministic 10-tip toy scenarios:

- `outputs/t10_global_deletion/`
  Global deletion series derived from the 10-tip benchmark tree.
- `outputs/t8_to_t3_local_deletion/`
  Local deletion series focused on the mid-tree collapse path from `t8` to `t3`.

Each scenario contains three benchmark-coordinate output locations:

- `benchmark_rooted/`
  Rooted oracle outputs retained as a biologically interpretable companion analysis.
- `benchmark_unrooted/`
  Unrooted oracle outputs used as the formal benchmark target for SplitAligner Perl.
- `splitaligner_perl/`
  Perl SplitAligner outputs generated from the same benchmark inputs.

## Why Both Rooted and Unrooted Outputs Exist

Both rooted and unrooted benchmark outputs are retained because they answer different questions.

- `benchmark_rooted/` preserves rooted branch distinctions and remains useful as a companion benchmark view.
- `benchmark_unrooted/` performs pseudo-root normalization and unrooted contraction so that branch identity is defined in the same way as the Perl SplitAligner workflow.

For this reason, `splitaligner_perl/` must be audited strictly against `benchmark_unrooted/`, not against `benchmark_rooted/`.

## Audit Rule

Pass/fail is defined only with respect to `benchmark_unrooted/`.

- Expected rooted-vs-unrooted differences are not treated as SplitAligner failures.
- A passing audit has zero unexpected mismatches between `splitaligner_perl/` and `benchmark_unrooted/`.

## Audit Outputs

Each scenario produces both oracle outputs and audit outputs.

Scenario oracle outputs:

- `benchmark_unrooted/benchmark.table.txt`
  Primitive-coordinate oracle matrix used as the formal audit target.
- `benchmark_unrooted/oracle_cell_status_long.tsv`
  Long-form primitive-cell oracle state table.
- `benchmark_unrooted/oracle_fusion_groups.tsv`
  Oracle fused-coordinate expectations for grouped `NA_fuse` paths.
- `benchmark_unrooted/benchmark.oracle_events.tsv`
  Event-level contraction motif record.
- `benchmark_unrooted/benchmark.trajectory.tsv`
  Trajectory summary table.
- `benchmark_unrooted/benchmark.schedule.tsv`
  Frozen deletion schedule used for the packaged scenario.

Scenario audit outputs:

- `audit/<scenario>/expected_vs_splitaligner.full.tsv`
  Full primitive-cell comparison table for all mapped cells.
- `audit/<scenario>/expected_vs_splitaligner.diff.tsv`
  Diff-only subset of the primitive-cell comparison table.
- `audit/<scenario>/unexpected_mismatches.tsv`
  Primitive-cell mismatches retained for quick inspection.
- `audit/<scenario>/fusion_groups.full.tsv`
  Full fused-coordinate audit table.
- `audit/<scenario>/fusion_groups.diff.tsv`
  Diff-only fused-coordinate audit table.
- `audit/<scenario>/pass_fail_summary.txt`
  Scenario-level pass/fail summary.
- `audit/<scenario>/endpoint_collapse_summary.txt`
  Tiny endpoint-collapse invariant summary for `1|k` internal-branch regression checks.
- `audit/<scenario>/endpoint_collapse_cases.tsv`
  Case-level endpoint-collapse checks, including singleton and terminal-fused cases.
- `audit/<scenario>/endpoint_collapse_internal_only_groups.tsv`
  Internal-only fused-path bookkeeping checks.

The current packaged audit validates:

- primitive-cell state agreement between `splitaligner_perl/` and `benchmark_unrooted/`
- active fused-coordinate agreement plus synthetic composite sum consistency
- endpoint-collapse invariants for `1|k` internal-branch cases

It does not currently package broader trajectory statistics such as ARI/Jaccard-style partition summaries.

The endpoint-collapse invariant test is implemented by:

- `scripts/test_endpoint_collapse_invariants.py`

This is a scenario-specific regression check for the packaged endpoint-collapse / local-deletion benchmark bundle. It uses fixed expected branch and gene identifiers from the current packaged scenarios and is not intended as a generic arbitrary-scenario validator.

This small regression script checks that:

- an endpoint-collapsed internal singleton is not emitted as a primitive numeric branch
- an endpoint-collapsed internal branch fused with a terminal branch produces a numeric fused coordinate and `NA_fuse` primitive members
- internal-only fused bookkeeping groups remain explicitly represented when they occur

In the current packaged audit, finite numeric Perl composite columns are separated into two cases. If all primitive component branches are finite numeric for a gene, the composite column is treated as a synthetic compatibility value and must equal the sum of its primitive components. If one or more primitive components are not finite numeric, the composite column is treated as an active fused-coordinate signal and must correspond to an oracle-expected merge group for the same gene. Composite columns with `NA`, `NA_fuse`, `NA_struct`, `NA_topo`, `NaN`, `Inf`, `-Inf`, or empty values are not treated as numeric fused-coordinate evidence.

## Directory Layout

- `docs/`
  Benchmark specification and implementation notes.
- `inputs/`
  Shared benchmark input trees and related source material.
- `scripts/`
  R-side oracle driver and utilities.
- `outputs/rooted_species_tree_branch_labels.pdf`
  Shared rooted 10-tip branch-label diagram used by both toy scenarios.
- `outputs/unrooted_species_tree_branch_labels.pdf`
  Shared unrooted 10-tip branch-label diagram used by both toy scenarios.
- `outputs/t10_global_deletion/`
  Scenario-specific output folders plus `branch_label_map.tsv`.
- `outputs/t8_to_t3_local_deletion/`
  Scenario-specific output folders plus `branch_label_map.tsv`.
- `audit/t10_global_deletion/`
  Audit tables comparing Perl outputs against `benchmark_unrooted/`.
- `audit/t8_to_t3_local_deletion/`
  Audit tables comparing Perl outputs against `benchmark_unrooted/`.

## Branch Label Maps

Each scenario includes a `branch_label_map.tsv` file with these columns:

- `splitaligner_branch`
- `benchmark_unrooted_branch`
- `benchmark_rooted_branch`
- `split`
- `note`

These maps are used to align branch coordinates across the Perl and benchmark outputs without altering the expected biological or mathematical results.

## Current Packaging State

The benchmark directory is scaffolded around the two packaged deterministic toy scenarios described above. Scenario folders may contain committed rooted, unrooted, Perl, and audit outputs for the current benchmark release.

When the benchmark workflow is rerun, generated rooted, unrooted, and Perl outputs should be written directly into the matching scenario subfolders above, and audit artifacts should be written into the corresponding `audit/` directory.

Random deletion remains a valid Benchmark V1 extension, but it is not part of the current packaged audit bundle.
