# Benchmark V1 Tasklist

## Goal

Build an oracle-validated fixed-topology pruning benchmark for `NA_struct` and `NA_fuse` on the semantics-normalized full-tree primitive-edge axis.

## Freeze

- full-tree primitive-edge axis includes terminal and internal edges
- external edge IDs are standardized to `Bxxx`
- V1 covers `NA_struct` and `NA_fuse` only
- ordered deletion is required for the packaged benchmark bundle
- ordered deletion uses a frozen plotting convention
- pruning stops at `max_delete_fraction = 0.7` and must satisfy `min_tips_remaining = 3`
- random deletion remains an optional extension for broader Benchmark V1 analyses

## Deliverables

- `Benchmark_V1_Spec.md`
- `README.md`
- `benchmark.table.txt`
- `oracle_cell_status_long.tsv`
- `oracle_fusion_groups.tsv`
- `benchmark.oracle_events.tsv`
- `benchmark.trajectory.tsv`
- `benchmark.schedule.tsv`
- ordered-pruning input gene-tree series
- SplitAligner comparison outputs
- audit outputs:
  - `expected_vs_splitaligner.full.tsv`
  - `expected_vs_splitaligner.diff.tsv`
  - `unexpected_mismatches.tsv`
  - `fusion_groups.full.tsv`
  - `fusion_groups.diff.tsv`
  - `pass_fail_summary.txt`
- stepwise plots
- trajectory plots
- random-pruning replicate series only if the optional random-deletion extension is enabled

## Work Items

### 1. Axis and mapping

- define the semantics-normalized full-tree primitive-edge axis from `T0`
- read or generate the authoritative `Bxxx` branch mapping
- ensure oracle-side edge IDs and SplitAligner-side edge IDs use the same axis

### 2. Schedule generation

- implement ordered deletion schedule
- freeze ladderize + cladewise + top-to-bottom plotting convention
- implement packaged deterministic toy scenarios
- implement random deletion schedule as a random tip permutation only if the optional extension is enabled
- enforce stopping rules: `max_delete_fraction = 0.7`, `min_tips_remaining = 3`
- store schedule metadata and seeds

### 3. Pruning-series generation

- generate the ordered pruning series `T0..Tk`
- generate random pruning replicate series only if the optional random-deletion extension is enabled
- store retained tip sets `S_i`
- export tree series in a format directly usable by SplitAligner

### 4. Event-level oracle

- for each deletion step, record `deleted_tip`
- record `new_na_struct_edge`
- record `local_e_term`
- record `local_e_sib`
- record `local_e_up`
- write `benchmark.oracle_events.tsv`

### 5. State-level oracle

- track primitive-edge survival on the semantics-normalized full-tree axis through pruning and contraction
- classify `NA_struct` by explicit node/edge state tracking
- derive fusion groups from contraction-induced merge classes
- write `benchmark.table.txt`
- write `oracle_cell_status_long.tsv`
- write `oracle_fusion_groups.tsv`

### 6. Closure checks

- verify all primitive edges are partitioned into `NA_struct` and non-`NA`
- verify all non-`NA` edges are fully partitioned by merge-class membership
- verify no overlaps and no omissions

### 7. Trajectory summaries

- compute `na_struct_count`
- compute `fusion_burden = Σ(|group| - 1)`
- compute `max_group_size`
- compute `num_groups_ge2`
- compute `n_tips`
- write `benchmark.trajectory.tsv`
- write `benchmark.schedule.tsv`

### 8. SplitAligner runs

- run SplitAligner on the ordered pruning series
- run SplitAligner on random replicate series only if the optional random-deletion extension is enabled
- collect matrix outputs on the same benchmark inputs

### 9. Packaged audit

- compare `benchmark_unrooted` primitive-cell states against `splitaligner_perl`
- audit active fused-coordinate agreement plus synthetic composite sum consistency against `benchmark.matrix_with_fuse.txt`
- write full audit tables and diff-only audit tables
- record scenario-level pass/fail summary

### 10. Optional broader Benchmark V1 comparisons

- compute `NA_struct exact match`
- compute `primitive-edge state exact match`
- compute `fusion partition exact match`
- compute `ARI`
- compute `Jaccard`

### 11. Visualization

- generate stepwise local-motif plots
- generate stepwise oracle-state plots
- generate ordered trajectory plots
- generate random `mean ± sd` trajectory plots only if the optional random-deletion extension is enabled
- generate ordered-vs-random comparison plots only if the optional random-deletion extension is enabled

### 12. Narrative alignment

- keep wording consistent with `Benchmark_V1_Spec.md`
- describe fusion as arising from subtree contraction
- describe SplitAligner as reporting or recovering fusion structure
- do not expand V1 to `NA_topo`

## Suggested Implementation Order

1. axis and `Bxxx` mapping
2. packaged deterministic schedules
3. pruning-series export
4. event-level oracle
5. state-level oracle
6. closure checks
7. trajectory summaries
8. SplitAligner runs
9. packaged primitive-cell and fused-coordinate audit
10. optional broader metrics
11. plots

## Definition of Done

- oracle outputs are generated on the semantics-normalized full-tree primitive-edge axis
- ordered packaged scenarios are reproducible, and random trajectories are reproducible if the optional random-deletion extension is enabled
- oracle closure passes at every step
- SplitAligner outputs are comparable on the same `Bxxx` axis
- packaged primitive-cell and fused-coordinate audits pass end-to-end
- benchmark outputs are ready to support manuscript and README updates
