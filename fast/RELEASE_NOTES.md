# SplitAlignerFast 0.2.0-rc.1

This Linux source candidate packages the C++17 exact-bitset core, C ABI 2, streaming Newick CLI and Python binding. Production C++ and Python code comes from source revision `57308470c1ebba9194a00638573bdf1efcf47e3c`. The value in the CLI manifest identifies that implementation source. Source archive checksums identify this distribution.

- Perform exact matching with all taxon-bitset words and retain complete projection membership and original source nodes.
- Reuse immutable plans across branch-length draws with the same topology and node numbering.
- Choose matrix, full mapping, topology-only support or coverage truth output. Dedicated unmatched exports are optional in full and support modes.
- Run the included C ABI, forced-collision and native regressions with CTest. Enabling Python 3.8 or later at configure time adds binding and CLI example tests.
- Install the executable, both library forms, headers, Python module and notice with CMake.

The supported target of this candidate is Linux. The CLI uses procfs to hash its executing binary. Each run handles one collection of gene trees and emits raw primitive states and optional detail tables. The paired fixed/free finalization is provided by the separate SplitAligner tools. See the [output guide](docs/OUTPUT.md) for state and support definitions.
