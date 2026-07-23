# Bundled example workflows

The bundled example assets and generated outputs use separate namespaces.

```text
302mammal/
  input/       immutable toy inputs
  expected/    immutable historical reference outputs
  run/         generated ownership-managed outputs

preprint_302mammal/
  input/       immutable 2,275-gene inputs
  run/         generated ownership-managed outputs
```

The `run/` directories are created when needed and are intentionally absent
from a fresh source archive.

The runner accepts only a real `run/` directory at the exact physical child
path shown above. It rejects a symbolic link, regular file, or special object
at `run/` before SplitAligner starts, so generated outputs cannot be redirected
into immutable `input/` or `expected/` assets or an external directory. A
repository checkout reached through a symbolic-link ancestor remains supported.

From the repository root:

```bash
bash examples/run.sh toy
bash examples/run.sh preprint
```

Each command writes every matrix, intermediate directory, manifest, classified
output, Support table, and annotated tree inside its corresponding `run/`
directory. Repeating the command is supported: `--force` may replace only the
outputs authorized by the manifests from the preceding owned run.

To start with a fresh workspace, preserve the current workspace as a uniquely
named sibling archive:

```bash
bash examples/archive_run.sh toy
bash examples/archive_run.sh preprint
```

`archive_run.sh` validates that `run/` is the exact physical direct child of the
selected example root, then atomically renames the entire directory to a name
such as `run.saved.20260721T120000Z.12345`. It never recursively deletes files.
An absent `run/` is a safe no-op; a symbolic link, non-directory, mounted or
cross-device workspace, identity change, collision, or failed rename is handled
without intentionally deleting the original contents. The reported archive
must be inspected before any manual deletion because it may contain unrelated
files that the ownership-aware workflow correctly preserved. A later
`examples/run.sh` invocation creates a fresh `run/` directory.

Do not delete or modify `input/` or `expected/`. Historical unowned generated
files are not adopted, and no ownership manifest is fabricated for them.
