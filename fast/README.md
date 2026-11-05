# SplitAlignerFast

SplitAlignerFast maps unrooted gene-tree edges onto a reference tree using exact taxon bitsets. It preserves the original source-node indices of a match and records how several reference edges can project onto one edge after pruning. The streaming `saf` command reads Newick trees. Applications can reuse reference, coverage and mapping plans through the C ABI or Python.

## Build on Linux

Use CMake 3.16 or later, a C++17 compiler, a C compiler, POSIX threads and the OpenSSL 1.1 or later development library. Python 3.8 or later runs the binding and optional tests. On Debian or Ubuntu the build packages are `build-essential`, `cmake` and `libssl-dev`.

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j2
(cd build && ctest --output-on-failure)
```

The executable is `build/saf` and the C ABI library is `build/libsplitaligner_fast.so`. To install both under a local prefix:

```sh
cmake --install build --prefix "$PWD/install"
```

## Run a quartet example

```sh
./examples/quickstart.sh ./build/saf quartet-output
cat quartet-output/raw_states.tsv
```

The two input trees share four taxa with the reference. The first has the same five unrooted splits. The second recovers all four terminal edges and misses the reference's internal split. Its own internal split is the unmatched gene edge. In `matrix.tsv`, the matched internal edge in the first tree has length 9, the sum of the two root-adjacent lengths 3 and 6. The second tree has `NA` on that reference axis. Use a fresh output directory for each run.

For your own data, supply a reference Newick file and a headerless TSV with one gene ID, one tab and one Newick tree on each line.

```sh
build/saf --reference species.nwk --genes genes.tsv --out aligned --scope full --role free
```

For piped genes, use `--genes /dev/stdin` (for example, `gzip -dc genes.tsv.gz | build/saf --reference species.nwk --genes /dev/stdin --out aligned`). The CLI reads each input once and hashes its original bytes. An optional UTF-8 BOM at the start of the reference, gene, or taxa file is accepted. If your existing file places the ID immediately before the opening `(` of each tree, select `--gene-format prefix-newick`. Every character before that parenthesis belongs to the ID, including spaces. The default format is tab-separated.

`--scope matrix` writes primitive values and states. Use `full` to retain coverage, mapping and composite records, `support` for topology-only free-tree support, and `truth` with `--taxa` to inspect the reference restriction for each taxon set. `--unmatched-splits off` omits the two dedicated unmatched exports. Use `--no-numeric` to skip branch-length conversion and omit the numeric matrix. The scope defaults and the two event counts are described in [the output guide](docs/OUTPUT.md).

For a free tree, read `composites.jsonl` in full mode to see whether each fused group was recovered. A recovered group has `group_state: "Mapped"`; an absent group has `"NA_fuse_topo"`. The member axes stay F in both cases. A missing branch length alone does not change the group state.

## Embed the core

The C header is `include/splitaligner/fast_c_api.h`. It exposes ABI 2 and caller-owned output buffers. The Python binding uses the same shared library through `ctypes`.

```sh
SAF_LIBRARY="$PWD/build/libsplitaligner_fast.so" PYTHONPATH="$PWD/adapters" python3 examples/reuse_mapping.py
```

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
