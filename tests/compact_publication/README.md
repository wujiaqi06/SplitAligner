# Compact publication regression

Run `LC_ALL=C perl tests/compact_publication/run.pl [new-evidence-directory]`.
Without a directory, fixtures are temporary. An explicit evidence directory must
be new or empty; failed-transaction backups and diagnostics remain there for
inspection. No large input or benchmark is run.

The driver loads the actual controller `publish_transaction`,
`publication_test_fail_after`, and `remove_path` bodies without executing matrix
construction. It uses the real ownership, provenance, namespace and device
modules. Failures are injected by locally wrapping calls at the relevant boundary,
without adding hooks to production code. This is a transaction test, not an
end-to-end scientific-output or release acceptance test.

Coverage follows the historical publication fault probes:

- P1: generated file/directory changed after manifest, before plan (clean and
  replacement); authority remains the manifest inventory.
- P2/P4: generated object changed after plan in clean/replacement runs.
- P3/P8: early/late replacement destination changed after plan, rejected before
  any output moves.
- P6/P7: owned/shared reuse destinations changed after plan, rejected before moves.
- P5: generated file/directory changed after all pre-publication checks. The
  changed installed object is preserved; the prior object remains in the reported
  recovery directory. This is preserved recoverability, **not automatic complete
  restoration**. Untouched outputs roll back normally.
- E10: an exception after moving a valid old object into backup restores the exact
  prior output. Actual backup tampering is instead retained and reported as
  incomplete recovery, never deleted or falsely claimed as restored.
- Manifest serialization seal: bytes changed before/after plan fail even if the
  altered JSON is semantically equivalent. The seal is produced from the intended
  serialization bytes, not by rereading mutable output.
- Unexpected inventory entries and ordinary late-failure rollback.

`source_binding.json` binds the loaded controller, two changed modules and test
driver. Each failing scenario records `error.txt`; P5/E10 recovery directories
remain in explicit evidence runs. The fixture records are synthetic and small.
The earlier `round1` receipt used swapped P6/P7 names and a nonhistorical P3 name;
it is retained as history. `round2` and later use the corrected probe meanings.
