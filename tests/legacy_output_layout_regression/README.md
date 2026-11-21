# Historical output-layout fixture

`baseline_25ab8e0.tar.gz` contains genuine Perl CLI output produced from Git
commit `25ab8e0`, using its unmodified `SplitAligner.pl`, `lib/` and `scripts/`.
It was generated during the 2026-11-06 review, before testing the current driver.
`fixture_provenance.json` records the full commit, commands, archive checksum,
and every input/output file checksum. The ownership manifests inside the archive
also record the exact implementation source hashes.

The small five-taxon fixture has two matrix runs (`free`, `fix`) and a paired
finalize run (`old`) with Support outputs. These genuine historical manifests
are not hand-edited or resealed to make the tests pass. The test runs without
Git or network access.

The current matrix adds `fusion_cells.tsv`, so the old matrix inventory does not
cover the current publication set. Omitting Support on a rerun of the archived
finalize creates the reverse case: unexpected old output paths. Both operations
must fail without changing any archived file and explain the fresh-label route.
Ordinary current-version force reruns must continue to work.
