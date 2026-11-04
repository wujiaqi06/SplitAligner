# Choosing a SplitAligner implementation

| What do you want to do? | Start here |
| --- | --- |
| Follow the paper's reference workflow, including paired fixed/free finalization | [SplitAligner in Perl](../README.md#quick-start) |
| Map a batch of trees on Linux and write matrices or per-edge mapping records | [SplitAlignerFast command line](../README.md#splitalignerfast-on-linux) |
| Reuse an edge mapping across branch-length draws in C, C++ or Python | [Fast library and Python examples](../fast/docs/LIBRARY.md) |
| Align trees, work with the results and pair fixed/free inputs in R | [SplitAlignerR package and installation](https://github.com/wujiaqi06/SplitAlignerR#installation) |

The Fast CLI processes one collection of gene trees per run. It writes raw mapping states, numeric values and, in full mode, records that identify the original edges behind a match. The Perl and R paired workflows take fixed and free collections and classify the final branch coordinates. To interpret a Fast matrix column, look up its split in that run's `context.json`. [Fast output guide](../fast/docs/OUTPUT.md)

Fast builds independently from `fast/`. From a repository checkout, change into that directory to follow its build instructions. If you downloaded the dedicated Fast source archive, start in the extracted directory. [Fast build instructions](../fast/README.md)

SplitAlignerR has its own releases, installation instructions and R documentation. Fast-backed versions of the R package ship the matching core sources and compile them during package installation. [R installation](https://github.com/wujiaqi06/SplitAlignerR#installation) · [R documentation](https://wujiaqi06.github.io/SplitAlignerR/)

In the layered releases, an absent shared split is recorded once on the free group. Read `<final_label>.free.projected_groups.tsv` in Perl, `composites.jsonl` in a Fast free full/support run, or `free_fusion_groups(free)` in R. The original reference branches remain fusion members.
