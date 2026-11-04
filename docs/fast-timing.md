# Timing the 302-mammal example

The [example data](../examples/preprint_302mammal/input) contain one reference species tree and two sets of trees for the same 2,275 genes. One set follows the reference topology (fixed); the other allows each gene its own topology (free). The timings below come from separate benchmarks with different workloads and computers.

## Perl: paired workflow

On 11 April 2026, the single-threaded Perl workflow processed both sets and calculated branch support in **9 min 52 s** on a 12-core Apple M3 Pro (macOS 14.4, Perl 5.34.1).

| Step | Elapsed time |
| --- | ---: |
| Build the free-tree matrix (2,275 genes) | 301.79 s |
| Build the fixed-tree matrix (the same 2,275 genes) | 286.92 s |
| Finalize the paired matrices and calculate branch support | 3.57 s |
| Total | 592.28 s |

These are elapsed times from one serial run of the April implementation, including intermediate file output. The trees had already been inferred.

## Fast: one tree set

Each Fast run below processed **one set** of 2,275 trees.

| Gene trees | Write the branch matrix (`--scope matrix`) | Write full output (`--scope full`) |
| --- | ---: | ---: |
| Fixed | 3.75 s | 7.84 s |
| Free | 3.76 s | 8.49 s |

Each entry is the median elapsed time of three runs on a shared Linux machine (Intel Xeon Gold 6326, one thread pinned to one logical CPU). GNU time measured the whole `saf` command, from reading the existing tree files through writing and checksumming the outputs. It excludes compiling the program and inferring the gene trees. The filesystem cache was not reset between runs.

These measurements are from the FAST-002 benchmark (CLI build `67de67fb742a`). The current CLI adds output options and reports recovery of fused branches, so its timings may differ.
