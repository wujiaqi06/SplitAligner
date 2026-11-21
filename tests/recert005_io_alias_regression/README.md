# I/O alias coverage in the compact layout

Run `LC_ALL=C LANG=C TMPDIR=/private/tmp bash tests/recert005_io_alias_regression/run.sh`
on macOS, or set `TMPDIR` to a literal existing temporary directory on POSIX.
The script uses a new private workspace and exits nonzero on any failed check.

The v2 matrix layout owns one aggregate directory, `<label>_split_branch_label`,
and three compact files: `<label>.compact.jsonl`, `<label>.compact.index.json`,
and `<label>.taxon_dictionary.json`. The old `_splits` publication-directory
fixture is replaced by exact-path, symlink-resolved and same-inode input aliases
against **each** compact file. The actual aggregate-directory nested-input
rejection remains tested, including Unicode paths and labels.

The original subsequent checks still execute: species aliases with/without
force, scientific matrix aliases, explicit/inferred manifests, coordinate
sidecars, Support tree/map, cross-role finalizers, duplicate destinations,
safe replacement, rollback, and input mutation after preflight. The low-level
I/O property test also retains arbitrary nested-directory and overlapping
destination checks independently of the matrix layout.

A bare, unowned `<label>_splits` directory is not a v2 destination. Its user
files must survive a successful new run and force rerun. This is distinct from
an actual frozen v1 run with its manifest, whose refusal is tested by
`tests/compact_pipeline/run.py` and the layout regression.
