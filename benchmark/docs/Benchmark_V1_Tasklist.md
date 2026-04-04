# Benchmark V1 Tasklist

## Goal

Build an oracle-validated fixed-topology pruning benchmark for `NA_struct` and `NA_fuse` on the full-tree primitive-edge axis.

## Freeze

- full-tree primitive-edge axis includes terminal and internal edges
- external edge IDs are standardized to `Bxxx`
- V1 covers `NA_struct` and `NA_fuse` only
- ordered deletion and random deletion are both required
- ordered deletion uses a frozen plotting convention
- random deletion uses permutation-based trajectories with fixed seed and `nrep = 100`
- pruning stops at `max_delete_fraction = 0.7` and must satisfy `min_tips_remaining = 3`

## Deliverables

- `Benchmark_V1_Spec.md`
- `oracle_states.tsv`
- `oracle_events.tsv`
- `trajectory.tsv`
- ordered-pruning input gene-tree series
- random-pruning replicate series
- SplitAligner comparison outputs
- stepwise plots
- trajectory plots

## Work Items

### 1. Axis and mapping

- define the full-tree primitive-edge axis from `T0`
- read or generate the authoritative `Bxxx` branch mapping
- ensure oracle-side edge IDs and SplitAligner-side edge IDs use the same axis

### 2. Schedule generation

- implement ordered deletion schedule
- freeze ladderize + cladewise + top-to-bottom plotting convention
- implement random deletion schedule as a random tip permutation
- enforce stopping rules: `max_delete_fraction = 0.7`, `min_tips_remaining = 3`
- store schedule metadata and seeds

### 3. Pruning-series generation

- generate the ordered pruning series `T0..Tk`
- generate random pruning replicate series
- store retained tip sets `S_i`
- export tree series in a format directly usable by SplitAligner

### 4. Event-level oracle

- for each deletion step, record `deleted_tip`
- record `new_na_struct_edge`
- record `local_e_term`
- record `local_e_sib`
- record `local_e_up`
- write `oracle_events.tsv`

### 5. State-level oracle

- track primitive-edge survival on the full-tree axis through pruning and contraction
- classify `NA_struct` by explicit node/edge state tracking
- derive fusion partitions from contraction-induced merge classes
- write `oracle_states.tsv`

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
- write `trajectory.tsv`

### 8. SplitAligner runs

- run SplitAligner on the ordered pruning series
- run SplitAligner on random replicate series as needed
- collect matrix outputs on the same benchmark inputs

### 9. SplitAligner-side fusion reconstruction

- reconstruct primitive-edge equivalence classes from fused composite labels
- incorporate primitive-row `NA_fuse` assignments
- derive the SplitAligner-side fusion partition per step

### 10. Comparison metrics

- compute `NA_struct exact match`
- compute `primitive-edge state exact match`
- compute `fusion partition exact match`
- compute `ARI`
- compute `Jaccard`

### 11. Visualization

- generate stepwise local-motif plots
- generate stepwise oracle-state plots
- generate ordered trajectory plots
- generate random `mean ± sd` trajectory plots
- generate ordered-vs-random comparison plots

### 12. Narrative alignment

- keep wording consistent with `Benchmark_V1_Spec.md`
- describe fusion as arising from subtree contraction
- describe SplitAligner as reporting or recovering fusion structure
- do not expand V1 to `NA_topo`

## Suggested Implementation Order

1. axis and `Bxxx` mapping
2. ordered schedule
3. pruning-series export
4. event-level oracle
5. state-level oracle
6. closure checks
7. trajectory summaries
8. SplitAligner runs
9. fusion-partition reconstruction
10. metrics
11. plots

## Definition of Done

- oracle outputs are generated on the full-tree primitive-edge axis
- ordered and random trajectories are reproducible
- oracle closure passes at every step
- SplitAligner outputs are comparable on the same `Bxxx` axis
- `NA_struct` and fusion recovery metrics can be computed end-to-end
- benchmark outputs are ready to support manuscript and README updates
