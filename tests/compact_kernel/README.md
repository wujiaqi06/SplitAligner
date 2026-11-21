# Compact kernel fixture checks

Run from the source root:

```sh
perl tests/compact_kernel/run.pl
```

For an independent comparison with a frozen source tree, set
`SPLITALIGNER_BASELINE_DIR` to that tree before running. The test reads its
`lib/SplitAligner/Newick.pm` into a separate in-memory namespace; it does not
modify or install the frozen source. The 2026-11-11 candidate was checked against
the Mac timing snapshot whose Perl files match commit
`21e397ada04702330f4b90842197461e1108eb88`.

For acceptance, use the explicit frozen-tree form and retain the command and
the frozen module hash with the log:

```sh
SPLITALIGNER_BASELINE_DIR=/absolute/path/to/frozen/21e397 \
  perl tests/compact_kernel/run.pl
```

At dev2 this executes **358 assertions**, including the assertion that the
external module loaded. Together with capacity (66), materialized spool (214)
and UTF-8 byte tests (11), the total is **649**. Without that environment
variable, this kernel suite executes **357**, and the four-suite total is
**648**. That local smoke mode compares with the candidate's own legacy
Newick helpers; it is not an independent frozen-source comparison. Record the
actual mode and counts rather than labelling either invocation simply "649".

These are small scientific and representation fixtures, not a full-panel timing
run. The checks cover byte widths and masks, numeric-scalar contamination,
canonical split identity, input-node provenance, binary display-root pooling,
Unicode/structured labels, primitive truth/state/values, complete component
construction, missing lengths, zeros, overflow, record corruption, and tree
object release. A subsequent full pipeline comparison remains necessary.

## Internal API

`prepare_reference($axis, $taxa)` takes the existing ordered axis records
(`branch_id`, `split`) and the taxon-label array. It returns a context containing
`axis`, `axis_ids`, `taxa`, `taxon_ids`, `reference_id`, `n`, and `byte_length`,
plus private byte-vector structures. It rejects incomplete terminal coverage.

`compute_gene($root, $ctx)` receives one already parsed tree. It does not parse
input again or modify the root. Its JSON-encodable record contains:

```text
schema = splitaligner-compact-gene-v1
reference_id
tg                  # fixed-width lowercase hexadecimal byte vector
source_node_count   # every non-root Newick preorder node/edge
observations[]      # {key, length, source_nodes[]}
components[]        # {members, key, observation_id}
truth[]
state[]
primitive[]
```

Observations are sorted by hexadecimal key. Their 1-based array position is the
observation ID; this is not a legacy split-text row. All input observations are
retained, even when no reference component matches. Each source node occurs
exactly once across the normalized observation units. Lengths and primitive
values are strings so original single-edge tokens survive JSON serialization.

Components are sorted by numeric branch-member pattern. An unmatched component
has a null `observation_id`. `truth`, `state`, and `primitive` follow `axis_ids`.
Reference-fused members remain `NA_fuse` in truth even when their compatibility
state is `U` because their component was not recovered.

`materialize_gene($record, $ctx)` validates the reference binding, bit vectors,
observation and component identities, source-node coverage, and stored derived
arrays before returning `{truth, state, primitive, groups, owners, patterns}`.
Groups use the existing canonical text split and are sorted by that text.
Only group output decodes the long taxon labels. This internal validation does
not replace replay against the original input for an audit of source-edge
lengths or topology.

`construct_cell($target_members, $owners, $view)` returns the existing eight
fusion-cell fields, in their established order, for `whole` or `direct`.
`compare_patterns($a, $b)` compares canonical `B...|B...` patterns numerically.

The bitwise module uses Core Perl 5.10-compatible syntax. It deliberately does
not enable the newer numeric `bitwise` feature and forces string operands in
its byte operations. Width, tail, Unicode-flag and ID checks protect the byte
contract; callers should not implement parallel raw `&`/`^` expressions.
