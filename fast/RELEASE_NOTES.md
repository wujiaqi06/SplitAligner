# SplitAlignerFast 0.2.0-rc.2

This Linux source candidate packages the C++17 exact-bitset core, C ABI 2, streaming Newick CLI and Python binding. The CLI now reports an unrecovered shared projection as `NA_fuse_topo` once per group in free full/support runs. Fused members keep raw state F. Production source: `91860d941c6b591b078fc99b02bdf7933fdba6c6`; the native core and ABI are unchanged. The value in the CLI manifest identifies that implementation source. Source archive checksums identify this distribution.

- Perform exact matching with all taxon-bitset words and retain complete projection membership and original source nodes.
- Reuse immutable plans across branch-length draws with the same topology and node numbering.
- Choose matrix, full mapping, topology-only support or coverage truth output. Dedicated unmatched exports are optional in full and support modes.
- Run the included C ABI, forced-collision and native regressions with CTest. Enabling Python 3.8 or later at configure time adds binding and CLI example tests.
- Install the executable, both library forms, headers, Python module and notice with CMake.

The supported target of this candidate is Linux. The CLI uses procfs to hash its executing binary. Each run handles one collection of gene trees and emits raw primitive states and optional detail tables. The paired fixed/free finalization is provided by the separate SplitAligner tools. See the [output guide](docs/OUTPUT.md) for state and support definitions.
