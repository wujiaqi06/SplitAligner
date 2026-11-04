# Use the native library

The shared library exposes the ABI in `include/splitaligner/fast_c_api.h`. Set `SAF_LIBRARY` to its full path before importing the Python adapter, or pass the path to `Library(path)`.

Pass one parent index and one taxon ID per node. Use parent -1 for the root and taxon -1 for internal nodes. The caller supplies the universe's taxon dictionary and must use it consistently in every plan. A reference contains all declared taxa. Each gene contains a duplicate-free subset equal to its coverage. Node numbering is the caller's original order.

Reference, coverage and mapping plans are immutable. A child keeps its parent alive. The Python binding supports context managers. The C caller destroys each handle with its matching destroy function. An existing mapping can be queried in parallel using separate output buffers. Mutable support updates, merges and handle destruction require caller synchronization.

`examples/reuse_mapping.py` illustrates two calls to `apply` on one mapping. The result holds separate primitive, projection and gene-edge families. Each family contains values and a parallel validity vector. Use the validity vector for unavailable values. The API accepts finite zero and negative values and checks for overflow when combining finite source lengths. C callers can query the required capacities and reuse their own buffers. On failure, inspect the returned status and copy the current thread's `saf_last_error()` text.

Inspect `ref.axes()` for primitive splits, `cov.truth()` for projection members, and `mapping.ledger()` for both match directions and all original gene-edge sources. Every plan has a local axis or projection order. Store the dictionary together with any exported results.
