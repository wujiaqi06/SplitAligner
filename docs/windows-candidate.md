# Windows candidate: validation pending

This isolated checkout contains native Windows adaptations and test entry
points. It is **not a Windows release**. Do not infer Windows support from the
Mac regression results. The previous frozen Mac/Linux sources and timing
measurements are unchanged.

## Perl: bounded first target

The initial native target is 64-bit Perl 5.40 or later, with `Win32API::File`
and `Win32`, on a local fixed NTFS volume. A full
[Strawberry Perl installation](https://strawberryperl.com/) is the intended
starting environment; the test must confirm its actual version and modules.
File identity checks remain enabled. Older Win32 `stat` implementations are
not accepted as substitutes for the required volume/file identity.

Use short ASCII paths for the source, work directory and command-line
arguments. Spaces are supported. UTF-8 taxon names and gene IDs are supported
**inside input files**; use a gene-ID file when auditing a non-ASCII ID. This
first port deliberately rejects UNC/network paths, SUBST drives, reparse
points (including junctions and symlinks), short-name aliases, alternate data
streams, reserved device names, trailing dots/spaces, non-ASCII filesystem
paths, and absolute paths of 240 or more characters. Publication paths must
be on the working-directory drive. These are explicit scope limits, not
silent path rewriting or reduced safety checks.

From the checkout root in native PowerShell, with Perl and Python available:

```powershell
./tests/windows_platform/run-native.ps1 -WorkDir C:/SA-check-perl-001
```

Choose a new results directory. The runner sets only its process-local
temporary-directory environment, restoring it afterwards. It records
platform/kernel/spool tests and a small end-to-end smoke test, including
whole/direct output, CRLF, UTF-8 contents, no-overwrite behavior and auditing.
Read both `RESULTS.json` and the individual logs. Failures and skipped tests
must remain visible. No large-panel timing claim is made by this check.


## Integration scope

This candidate combines the reviewed Perl v2 diagnostic and test-layout repairs
with the earlier native Windows port and accepted PowerShell runner correction.
Native Windows acceptance of this combined source remains pending. Fast and R
are maintained separately; this checkout is staged for Perl validation only.
