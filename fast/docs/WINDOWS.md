# Native Windows source candidate

The preceding Windows source was tested with Rtools GCC/MinGW on native
Windows x64. The current candidate incorporates that port and its returned
test fixes. CI tests the assembled source with UCRT64 GCC/MinGW; consult the
run results for this revision rather than inheriting a previous PASS. MSVC is
not claimed as validated. WSL is a separate Linux environment.

## Build and test

The initial target is 64-bit Windows, a C++17 compiler, CMake/CTest 3.26 or later
for this verification script, OpenSSL development headers/libraries of the same
architecture/toolchain, and 64-bit Python 3.8 or later.
The verification script accepts an explicit CMake generator. For the GCC/MinGW
target use `-Generator "MinGW Makefiles"` with matching GCC, Make and OpenSSL
on PATH. Visual Studio builds are an additional, unqualified target.
The script does not install dependencies or change system settings. An existing
OpenSSL installation can be supplied as `-OpenSSLRoot`; a configured dependency
toolchain can instead be supplied as `-ToolchainFile`.

From the repository root, in PowerShell:

```powershell
./fast/tests/verify_windows.ps1 -Results C:/SplitAlignerCheck/fast-001 -Generator "MinGW Makefiles" -OpenSSLRoot C:/path/to/openssl
```

Use a new local work directory. `C:/path/to/openssl` is a placeholder for the
actual installation, not a bundled dependency. A first run on local NTFS avoids
confounding network-share and removable-drive behavior with native support.
This is a source/build test, not a precompiled Windows distribution.

The script records versions, source hashes, compiler/build logs, the C ABI and
Python-binding checks, and CLI results. It installs into its own results folder
and tests the installed CLI, DLL and Python adapter there. Review `ctest.log` for skipped probes:
POSIX FIFO/executable-replacement tests do not run on Windows; symlink tests
skip only if Windows reports the missing symlink privilege. These skipped
capabilities are not counted as validated.

## What changed

- Wide command-line arguments and UTF-8 filesystem conversion preserve Unicode
  paths; outputs remain binary-mode UTF-8/LF, and input hashes include CRLF.
- Decimal parsing uses a private C numeric locale on Windows.
- DLL symbols are exported and the static library has a distinct Windows name.
- Portable population counts avoid requiring a CPU-specific instruction.
- The fusion spool is closed and removed before publishing the output folder.
- Windows directory publication uses `MoveFileExW` without replacement or copy
  flags. Existing targets must remain untouched; OS errors are retained.

Use `--genes -` or `--taxa -` for standard input. The reference remains a file.
The manifest records `"-"` as the stream name and hashes the actual bytes read.
The Python adapter accepts an explicit `SAF_LIBRARY` pointing to the built DLL.
Default discovery also checks the build directory and its Release/Debug folders.

The two-quartet command-line example has a native Python entry point:
`python examples/quickstart.py PATH_TO_SAF NEW_OUTPUT_DIRECTORY`. It uses the
same inputs and CLI options as `quickstart.sh`, which requires a POSIX shell.

## Acceptance boundary

Native testing must cover build/install, Unicode and space paths, input hashes,
whole/direct/off views, missing/zero lengths, source provenance, temporary-file
cleanup and no-overwrite publication. Scientific outputs are compared with
fixed expectations, not merely checked for a zero exit status. Long paths,
network shares, unusual filesystem configurations, Windows ARM and additional
compiler families require separate evidence; this candidate does not certify
them by inference.

For development on POSIX, `-DSAF_TEST_PORTABILITY_BRANCHES=ON` exercises the
portable population count and close-before-delete spool. This is useful
regression coverage and still **does not execute the Windows APIs**.

Implementation references: [Microsoft MoveFileExW](https://learn.microsoft.com/en-us/windows/win32/api/winbase/nf-winbase-movefileexw)
and [CMake Windows symbol exports](https://cmake.org/cmake/help/latest/prop_tgt/WINDOWS_EXPORT_ALL_SYMBOLS.html).
