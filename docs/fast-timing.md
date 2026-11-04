# Timing the 302-mammal example

The [example data](../examples/preprint_302mammal/input) contain one reference species tree and two sets of trees for the same 2,275 genes. One set follows the reference topology (fixed); the other allows each gene its own topology (free). Each time below is for **one set**.

| Gene trees | Write the branch matrix (`--scope matrix`) | Write full output (`--scope full`) |
| --- | ---: | ---: |
| Fixed | 3.75 s | 7.84 s |
| Free | 3.76 s | 8.49 s |

Each entry is the median elapsed time of three runs on a shared Linux machine (Intel Xeon Gold 6326, one thread pinned to one logical CPU). GNU time measured the whole `saf` command, from reading the existing tree files through writing and checksumming the outputs. It excludes compiling the program and inferring the gene trees. The filesystem cache was not reset between runs.

These measurements are from the FAST-002 benchmark (CLI build `67de67fb742a`). The current CLI adds output options and reports recovery of fused branches, so its timings may differ.
