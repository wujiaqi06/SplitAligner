# SplitAlignerFast

SplitAlignerFast maps unrooted gene-tree edges onto a reference tree using exact taxon bitsets. It preserves the original source-node indices of a match and records how several reference edges can project onto one edge after pruning. The streaming `saf` command reads Newick trees. Applications can reuse reference, coverage and mapping plans through the C ABI or Python.

Version: **0.2.0-rc.3**, C ABI 2.

## Precompiled Apple Silicon binary

The macOS arm64 archive on the [release page](https://github.com/wujiaqi06/SplitAligner/releases)
requires macOS 14.4 or later. Extract it and run `bin/saf --version`. It
includes a static OpenSSL dependency, so the CLI needs neither Homebrew nor
Python. The archive also contains the C ABI libraries and optional Python
adapter. Read `README-BINARY.txt` and `BUILD_INFO.json` in the archive for
build and signing details. Intel Macs require a source build.

## Build on Linux or macOS

Use CMake 3.16 or later, a C++17 compiler, a C compiler, POSIX threads and the OpenSSL 1.1 or later development library. Python 3.8 or later runs the binding and optional tests. On Debian or Ubuntu the build packages are `build-essential`, `cmake` and `libssl-dev`.

On macOS, use the Xcode Command Line Tools and a native OpenSSL installation. If CMake cannot locate a Homebrew OpenSSL installation, add `-DOPENSSL_ROOT_DIR="$(brew --prefix openssl@3)"` to the configure command. This source supports native Apple Silicon builds; availability of precompiled distributions is separate.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j2
(cd build && ctest --output-on-failure)
```

The executable is `build/saf`. The C ABI library is `build/libsplitaligner_fast.so` on Linux and `build/libsplitaligner_fast.dylib` on macOS. To install both under a local prefix:

```sh
cmake --install build --prefix "$PWD/install"
```

For native Windows GCC/MinGW, use the [Windows build guide](docs/WINDOWS.md).
MSVC and Windows ARM are not qualified by the existing tests.

## Run a quartet example

```sh
./examples/quickstart.sh ./build/saf quartet-output
cat quartet-output/raw_states.tsv
cat quartet-output/support.nwk
```

The two input trees share four taxa with the reference. The first has the same five unrooted splits. The second recovers all four terminal edges and misses the reference's internal split. Its own internal split is the unmatched gene edge. In `matrix.tsv`, the matched internal edge in the first tree has length 9, the sum of the two root-adjacent lengths 3 and 6. The second tree has `NA` on that reference axis. The support tree is `((A,B)50,(C,D)50);`. The reference's internal split is present in one of two eligible genes. Its two labels show where the display root divides the same branch. Use a fresh output directory for each run.

For your own data, supply a reference Newick file and a headerless TSV with one gene ID, one tab and one Newick tree on each line.

```sh
build/saf --reference species.nwk --genes genes.tsv --out aligned --scope full --role free
```

For piped genes, use `--genes /dev/stdin` (for example, `gzip -dc genes.tsv.gz | build/saf --reference species.nwk --genes /dev/stdin --out aligned`). The CLI reads each input once and hashes its original bytes. An optional UTF-8 BOM at the start of the reference, gene, or taxa file is accepted. If your existing file places the ID immediately before the opening `(` of each tree, select `--gene-format prefix-newick`. Every character before that parenthesis belongs to the ID, including spaces. The default format is tab-separated.

Every run with gene trees writes `support.tsv` and the annotated reference tree `support.nwk`, including `--role fixed` and `--scope matrix`. [Compare fixed and free support in a seven-tip example](docs/OUTPUT.md#branch-support).

`--scope matrix` writes original-branch values and states, plus the shared fused-group matrix and its cell provenance. In a numeric run, a group column may contain one recovered edge length or the sum of an exact set of recovered whole components. Choose `--fusion-view direct` for single-edge observations, or `--fusion-view off` for original-branch values alone. Use `full` to retain coverage, mapping and composite records, `support` for topology-only support, and `truth` with `--taxa` to inspect the reference restriction for each taxon set. `--unmatched-splits off` omits the two dedicated unmatched exports. Use `--no-numeric` to skip branch-length conversion and omit the numeric matrices. The scope defaults and the two event counts are described in [the output guide](docs/OUTPUT.md).

For a free tree, read `composites.jsonl` in full mode to see whether each fused group was recovered. A recovered group has `group_state: "Mapped"`; an absent group has `"NA_fuse_topo"`. The member axes stay F in both cases. A missing branch length alone does not change the group state.

## Embed the core

The C header is `include/splitaligner/fast_c_api.h`. It exposes ABI 2 and caller-owned output buffers. The Python binding uses the same shared library through `ctypes`.

```sh
SAF_LIBRARY="$PWD/build/libsplitaligner_fast.so" PYTHONPATH="$PWD/adapters" python3 examples/reuse_mapping.py
```

On macOS, use `libsplitaligner_fast.dylib` in that command. The Python binding also chooses the native suffix automatically when loading from the default `build` directory.

A reference owns the complete taxon universe. A coverage plan restricts that reference to the taxa observed for one gene. A mapping binds one gene topology and its *original node order* to that coverage. Reuse the mapping for new lengths on the same nodes. When topology, root representation or node order changes, construct a new mapping. See [the library guide](docs/LIBRARY.md).

## License and citation

The project's software is distributed under the [MIT License](LICENSE). The CLI links the system OpenSSL crypto library for SHA-256 hashing. OpenSSL retains its own license. Software citation metadata is in [CITATION.cff](CITATION.cff).

## Acknowledgments

We are especially grateful to Chao Zhang (张超;
[@chaoszhang](https://github.com/chaoszhang)) of the School of Life Sciences,
Peking University, for his generous technical help and for sharing his
`fastSplitAligner.py` implementation. His
implementation provided an important reference for efficient split-based
branch mapping and a comparison baseline during the development and
validation of SplitAlignerFast.
