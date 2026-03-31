# Benchmark V1 Spec

## 1. Scope and Non-goals

Benchmark V1 evaluates only `NA_struct` and `NA_fuse` under fixed-topology pruning. It does not evaluate `NA_topo`, species-tree inference accuracy, or generic phylogenetic accuracy.

The purpose of this benchmark is to validate SplitAligner's core object:

- a reproducible species-tree branch coordinate under progressive taxon loss
- structurally explicit missingness under progressive taxon loss
- fusion as contraction-induced edge merging rather than a software-invented artifact

This benchmark is therefore an oracle-validated missingness benchmark, not a tree-reconstruction benchmark.

Fusion groups are not heuristic objects introduced by the software, but contraction-induced merged-edge classes on the fixed full-tree primitive-edge axis.

## 1.1 Oracle Independence Constraint

Benchmark V1 requires the oracle implementation to remain structurally independent of SplitAligner.

This is a hard rule.

The R-side benchmark oracle must not reimplement SplitAligner's split-based missingness logic, must not use projected splits to infer `NA_struct` or `NA_fuse`, and must not act as a second R version of SplitAligner.

Instead, the R-side oracle is restricted to APE-based tree surgery and explicit node-edge relations, including:

- `drop.tip`-based pruning
- degree-2 contraction
- parent-child adjacency
- explicit local and cumulative edge-state tracking on the full-tree primitive-edge axis

In Benchmark V1:

- split-based branch mapping is reserved for the Perl SplitAligner implementation
- R is used only to construct an independent structural oracle for `NA_struct` and `NA_fuse`

If the R-side benchmark begins to use split-based logic to classify missingness, then the oracle ceases to be independent and the benchmark loses its methodological value.

## 2. Core Object

The unified coordinate axis is the set of all primitive edges on the full species tree `T0`, including both terminal and internal edges.

All oracle outputs, SplitAligner outputs, and evaluation metrics must be mapped back to this full-tree primitive-edge axis. External branch identifiers should be standardized to SplitAligner branch labels (`Bxxx`) using the species-tree branch map as the authoritative mapping source.

## 3. Definitions

### 3.1 Full tree and pruning series

Let `T0` be the full species tree with tip set `S0`.

Benchmark data are generated as a nested pruning trajectory:

- `T0 = full tree`
- `T1 = drop one tip from T0 and contract degree-2 nodes`
- `T2 = drop one additional tip from T1 and contract degree-2 nodes`
- ...
- `Tk`

Each `Ti` is treated as a fixed-topology gene tree: topology is always the induced subtree of `T0`, and only the retained taxon set changes.

### 3.2 Primitive-edge axis

For each primitive edge `e` in `T0`, assign a stable edge identity on the full-tree axis. This axis includes both terminal and internal edges and remains the only benchmark coordinate system used for oracle outputs, SplitAligner outputs, and comparison metrics.

### 3.3 Structural state under pruning

At pruning step `i`, tree `Ti` is produced from `T0` by repeated tip deletion and degree-2 contraction.

For each primitive edge `e` on the full-tree axis, the node-based oracle assigns one of two structural outcomes:

- `survives`: edge identity is still represented through a surviving contracted edge in `Ti`
- `NA_struct`: edge identity is no longer structurally represented after pruning and contraction

In Benchmark V1, terminal primitive edges are retained on the same full-tree axis as internal edges. Once the incident tip of a terminal primitive edge is absent from the retained taxon set, that edge becomes structurally undefined and must be labeled `NA_struct`.

### 3.4 Merge classes and fusion truth

At each step `i`, the contracted tree `Ti` contains a set of surviving edges. Under the node-based oracle, each surviving edge in `Ti` is associated with the set of primitive edges from `T0` whose identity has collapsed onto that surviving edge.

This set is called a merge class.

Interpretation:

- merge class size `= 1`: uniquely surviving primitive edge
- merge class size `>= 2`: fusion group

Fusion truth in Benchmark V1 is therefore defined by contraction-induced merge classes derived from explicit node/edge tracking, not by projected-split equivalence.

### 3.5 Local contraction motif

Under the degree-2 contraction rule used by `drop.tip` in the R `ape` package, each tip deletion yields a local contraction motif at the parent node of the deleted tip.

For a deleted tip `x`, define:

- `e_term(x)`: the incident terminal primitive edge of `x`
- `e_sib(x)`: the neighboring edge descending from the same parent
- `e_up(x)`: the immediate upstream edge from that parent

This local motif is recorded as an event-level object for explanation, sanity checking, and visualization. It is not itself the final global fusion truth.

It records the newly induced local contraction event at that deletion step, whereas final branch states are determined only by the state-level node/edge oracle on the full primitive-edge axis.

## 4. Benchmark Data Generation

### 4.1 Canonical ordered tip list

Benchmark V1 first defines a single deterministic ordered tip list on the full tree. This list is the basis for all ordered schedules.

The canonical ordered tip list should be constructed by:

- ladderizing the tree
- reordering cladewise
- extracting the final tip plotting order from top to bottom

This ordered tip list is benchmark-defined and visualization-dependent. Once frozen, the plotting convention and resulting tip order must remain fixed across all downstream analyses.

### 4.2 Canonical ordered deletion

Canonical ordered deletion starts from the beginning of the frozen ordered tip list and deletes tips one by one in that order.

This remains the primary benchmark trajectory because it provides a single canonical pruning path on which contraction propagation and fusion growth can be tracked consistently.

### 4.3 Anchored local ordered deletion

In addition to canonical ordered deletion, Benchmark V1 permits anchored local ordered deletion on the same frozen ordered tip list.

An anchored local ordered schedule is defined by:

- `start_index`: the position of the starting tip in the frozen ordered tip list
- `direction`: `forward` or `backward`
- `delete_count` or `delete_fraction`

The deletion path is then the contiguous window on the frozen ordered tip list induced by that start position and direction.

Examples:

- `forward`: `t_s, t_{s+1}, t_{s+2}, ...`
- `backward`: `t_s, t_{s-1}, t_{s-2}, ...`

This schedule family is intended for local mechanism analysis, sensitivity analysis, and targeted fusion benchmarks. It does not replace canonical ordered deletion as the default primary trajectory.

### 4.4 Maximum deletion length and stopping rules

Benchmark V1 does not delete tips until tree collapse.

The stopping constraints are:

- `max_delete_fraction = 0.7`
- `min_tips_remaining = 3`

For any schedule, the maximum allowed deletion length is:

- `k_max = min(direction_available_length, floor(n_tips * max_delete_fraction), n_tips - min_tips_remaining)`

where `direction_available_length` is the number of deletable tips available from the chosen `start_index` in the chosen direction, including the start position itself.

Using 1-based indexing on the frozen ordered tip list:

- `forward`: `direction_available_length = n_tips - start_index + 1`
- `backward`: `direction_available_length = start_index`

If `delete_fraction` is used, it must be converted into a deletion count by:

- `delete_count = floor(n_tips * delete_fraction)`

and then truncated by `k_max`.

If `delete_count` is used directly, it must satisfy `delete_count <= k_max`.

### 4.5 Random deletion

Random deletion is generated by first drawing a random permutation of the full tip set, then deleting tips in that order one by one.

Default settings:

- `nrep = 100`
- fixed random seed
- `max_delete_fraction = 0.7`
- `min_tips_remaining = 3`

Random trajectories are summarized across replicates using `mean ± sd`.

The same stopping rules apply to canonical ordered deletion, anchored local ordered deletion, and random deletion.

Schedule roles within Benchmark V1 are therefore:

- canonical ordered deletion: primary benchmark trajectory
- anchored local ordered deletion: local mechanism, sensitivity, and targeted-fusion trajectories
- random deletion: statistical control trajectories summarized across replicates

## 5. Oracle Structure

The oracle has two levels.

### 5.1 Level 1: event-level local oracle

For each deletion step:

- identify the deleted tip
- record the incident terminal edge
- record the local contraction motif `{e_term, e_sib, e_up}`

Under the benchmark contraction rule, each tip deletion deterministically renders the incident terminal branch structurally undefined on the retained taxon set, yielding at least one `NA_struct` event.

This level is used for:

- mechanism explanation
- stepwise sanity checking
- local visualization

### 5.2 Level 2: state-level node/edge oracle

For each step `i`, the oracle operates only through explicit tree surgery and node/edge tracking:

- prune the specified tip from the current tree
- contract degree-2 nodes
- track which primitive edges from `T0` are still represented by each surviving edge in `Ti`
- mark primitive edges that are no longer represented by any surviving edge as `NA_struct`
- define fusion groups as surviving merge classes with size `>= 2`

This is the benchmark ground truth.

## 6. Closure and Sanity Conditions

For each step:

- all primitive edges must be partitioned into `NA_struct` or exactly one surviving merge class
- all surviving primitive edges must be completely partitioned by merge-class membership
- there must be no overlaps and no omissions

This provides the oracle-side closure check for the benchmark.

## 7. SplitAligner Outputs and Mapping

SplitAligner outputs must be mapped to the same full-tree primitive-edge axis as the oracle.

The authoritative external branch ID is `Bxxx`, derived from the species-tree branch map generated by SplitAligner, typically via `species_tree.branch_map.txt`.

Oracle-side edge IDs and SplitAligner outputs must not be compared using raw edge indices from the R `ape` package alone.

## 8. Evaluation Targets

Benchmark V1 compares SplitAligner against the oracle on three levels.

### 8.1 Structural missingness recovery

For each step, compare the SplitAligner `NA_struct` set against the oracle `NA_struct` set on the primitive-edge axis.

Primary goal:

- exact edge-level match

### 8.2 Fusion recovery

For each step, compare the oracle fusion partition against SplitAligner's fused representation.

SplitAligner-side fusion recovery should be evaluated by reconstructing, for each step, the primitive-edge merge classes implied by its fused composite labels and primitive-row `NA_fuse` assignments, and comparing that partition against the oracle fusion partition.

Primary goal:

- structure-level fusion partition recovery

Secondary goal:

- consistency of primitive-row `NA_fuse` labeling with the recovered fusion structure

### 8.3 Trajectory consistency

Compare ordered and random pruning trajectories in terms of structural and fusion burden.

## 9. Metrics

Primary metrics:

- `NA_struct exact match`
- `primitive-edge state exact match`
- `fusion partition exact match`
- `ARI`
- `Jaccard`

Trajectory statistics:

- `na_struct_count(i)`
- `fusion_burden(i) = Σ(|group| - 1)`
- `max_group_size(i)`
- `num_groups_ge2(i)`
- `n_tips(i)`

Ordered deletion is summarized as a single deterministic trajectory. Random deletion is summarized as `mean ± sd` across replicates.

For `primitive-edge state exact match`, each primitive edge at each step is compared in terms of its benchmark-visible state on the full-tree axis, such as `mapped`, `NA_struct`, or `NA_fuse`.

## 10. Required Output Tables

### 10.1 `oracle_states.tsv`

Recommended fields:

- `step_id`
- `deleted_tip`
- `n_tips`
- `edge_id`
- `edge_type`
- `state_label`
- `is_struct_undef`
- `merge_group_id`

This is the main long-format oracle truth table.

### 10.2 `oracle_events.tsv`

Recommended fields:

- `step_id`
- `deleted_tip`
- `new_na_struct_edge`
- `local_e_term`
- `local_e_sib`
- `local_e_up`

This stores the event-level local motif record for each deletion step.

### 10.3 `trajectory.tsv`

Recommended fields:

- `step_id`
- `n_tips`
- `na_struct_count`
- `fusion_burden`
- `max_group_size`
- `num_groups_ge2`
- `schedule_type`
- `replicate_id`

This stores trajectory-level summaries for ordered and random schedules.

## 11. Benchmark Questions

Benchmark V1 must answer the following questions:

1. Can SplitAligner exactly recover oracle-defined `NA_struct` under fixed-topology pruning?
2. Can SplitAligner recover oracle-defined fusion groups on the same primitive-edge axis?
3. Does ordered pruning produce a distinct fusion trajectory relative to random pruning?
4. How does the local contraction motif relate to the global merge-class growth pattern across steps?
5. Does oracle closure hold at every pruning step?

## 12. Implementation Order

Recommended implementation order:

1. define and freeze the full-tree branch axis
2. generate ordered pruning series
3. generate random pruning replicates
4. implement event-level local oracle
5. implement state-level node/edge oracle
6. write `oracle_states.tsv`, `oracle_events.tsv`, and `trajectory.tsv`
7. run SplitAligner on the same pruning series
8. map SplitAligner outputs to the `Bxxx` axis
9. compare oracle and SplitAligner outputs
10. generate trajectory plots and stepwise diagrams

## 13. Narrative Constraints

Documentation, README text, and manuscript wording should remain consistent with this benchmark design.

Recommended wording principles:

- do not describe fusion as something invented by the software
- describe fusion as arising from subtree contraction under progressive taxon loss
- describe SplitAligner as reporting, representing, or recovering fusion structure
- do not frame this benchmark as a species-tree inference benchmark
- do not expand V1 to `NA_topo`

Suggested narrative summary:

> Benchmark V1 is an oracle-validated fixed-topology pruning benchmark designed to evaluate branch-state persistence, structural missingness (`NA_struct`), and contraction-induced merge classes (`NA_fuse`) on a unified full-tree primitive-edge axis.
