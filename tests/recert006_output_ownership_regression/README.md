# Output ownership coverage in the compact layout

Run `LC_ALL=C LANG=C TMPDIR=/private/tmp bash tests/recert006_output_ownership_regression/run.sh`
on macOS, or use a literal existing POSIX temporary directory. All 16 regression
groups must complete. The final coverage line reports executed CLI commands and
rejection cases; these counts are not counts of individual assertions.

The suite inventories the real v2 destinations: one aggregate directory and the
compact source ledger, index and taxon dictionary. Each compact file is checked
as an unowned regular file, an unowned directory containing user data, a modified
owned file, a missing owned file, a directory substituted for an owned file, a
mode-changed file, a symlink and a hardlink. All rejected attempts must preserve
the prior filesystem snapshots and leave no transaction workspace.

The retained aggregate directory still receives recursive add/modify/delete,
type/mode and nested-link substitutions. Later original groups remain active:
manifest integrity, finalizer ownership, safe force replacement, common-axis
reuse/conflict, Unicode identity, injected late rollback, exact inventories and
byte-deterministic manifests. All three compact files are included in whole-run
snapshot and rollback assertions.

The old `_splits` pathname is no longer a matrix publication target. An unrelated
directory using that spelling is explicitly preserved during successful new and
force runs; this does not authorize migration of an actual v1 result directory.
The frozen-v1 result/manifest refusal is covered by the compact pipeline suite.
