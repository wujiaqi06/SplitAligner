# Representation / input-variant regression battery

- SplitAligner commit: `ca09e46`
- Date: 2026-10-08
- Data: `examples/302mammal` (toy: 302-taxon species tree, 7 free + 7 fix gene trees)
- Generator: `tests/representation_variant_regression/variant_battery.py`; inputs in `variants/`

Naming: `SxG` = species-tree variant S x gene-tree variant G. A1 = bundled baseline (species tree and all gene trees are already written unrooted, with a trifurcating root).

| Group | Variant | Result | Evidence |
|---|---|---|---|
| A | A1xA1 | PASS | baseline equals examples/302mammal/expected |
| A | A2xA1 | PASS | identical-after-mapping (B aliases remapped) |
| A | A3xA1 | PASS | identical-after-mapping (B aliases remapped) |
| A | A4xA1 | PASS | identical-after-mapping (B aliases remapped) |
| A | A5xA1 | PASS | identical-after-mapping (B aliases remapped) |
| A | A6xA1 | PASS | identical-after-mapping (B aliases remapped) |
| A | A1xA2 | PASS | byte-identical |
| A | A1xA3 | PASS | byte-identical |
| A | A1xA3b | PASS | byte-identical |
| A | A1xA4 | PASS | byte-identical |
| A | A1xA5 | PASS | byte-identical |
| A | A1xA6 | PASS | byte-identical |
| A | A5xA5 | PASS | identical-after-mapping (B aliases remapped) |
| A | A2xA4 | PASS | identical-after-mapping (B aliases remapped) |
| A | A6xA6 | PASS | identical-after-mapping (B aliases remapped) |
| A | A3xA2 | PASS | identical-after-mapping (B aliases remapped) |
| B | B1 no branch lengths | PASS | state sidecars byte-identical to baseline=True; numeric cells in free matrix=0; finalize ran, support rows=601, e.g. ['B1', 'terminal', '5', '0', '0', '0.0000000000', '0.0000000000'] |
| B | B2 no internal-node labels | PASS | byte-identical  |
| B | B3 multifurcating gene trees | PASS | 445 internal nodes collapsed across 7 free trees; every changed cell is a collapsed split -> U/NA |
| B | B4a degree-2 root, both half-lengths numeric | PASS | same inputs as A1xA3: byte-identical |
| B | B4b degree-2 root, one half-length missing | PASS | every changed cell is the root split (or a fused column synthesized from it) -> NA, in 18 gene/file pairs; states unchanged |
| B | B5 gene taxa strict subset of species taxa | PASS | species tips=302; gene tips range 255-300 over 14 trees (all strict subsets; covered by every run above) |

## Variant definitions

- A2 root on a tip: degree-2 root placed on a terminal edge; that edge's length is kept on the tip half and the other half is written `:0`.
- A3 root on an internal edge: as A2, on an internal edge chosen with a fixed seed.
- A3b: as A3 but each half carries half the length (`%.15g` of x/2), so the merged cell is re-summed in floating point.
- A4 trifurcation on another centre: the unrooted tree re-written with a different internal node as the trifurcating root; no length string changes.
- A5 ladderize flip: every child list reversed.
- A6 random child order: every child list shuffled with a fixed seed.
- B1: every `:length` removed from free and fix gene trees.
- B2: every internal-node label (support value) removed.
- B3: a seeded 25% of non-root internal nodes of each free tree collapsed into polytomies.
- B4a/B4b: degree-2 root on an internal edge with both half-lengths numeric (B4a = A3) or with one half-length omitted (B4b).
- B5: the bundled gene trees each sample a strict subset of the 302 species.

## Reading the results

- `byte-identical`: all ten matrix/state/classified outputs and the support table are identical byte for byte.
- `identical-after-mapping`: the species-tree serialization changed, so the textual `B` aliases differ; after mapping each alias to its canonical split key via `primitive_axis.tsv`, every cell is identical.
- `PASS (numeric)`: the only differences are floating-point re-summation of a degree-2 root's two half-lengths (equal to 1e-12 relative); identity, states and classification are unchanged.

## Caveats recorded by the battery

- A2/A3 place the new degree-2 root on an edge whose length string is `%.15g`-stable (568 of the 7874 length strings in the toy gene trees are not, e.g. `0.0301476780` or `0.0000023606`). Re-summing the two half-lengths of a degree-2 root goes through Perl's float stringification, so on a non-stable edge the merged cell would be numerically equal but not byte-identical (`PASS (numeric)`). This is a formatting property of the summed cell, not a change of identity, state or classification. A3b (halved lengths) happened to re-sum to byte-identical strings on this data.
- B1 (no branch lengths): `matrix` and `finalize` both complete; every state sidecar is byte-identical to the baseline and every value cell is `NA`. `final.support_b.txt` then reports `n_fix_non_na = n_free_non_na = 0` and prints `support_percent = 0.0000000000` for a 0/0 ratio rather than `NA`. Downstream readers should treat Support with a zero denominator as undefined.
- B4b: with one root half-length missing, exactly the root split's primitive cell (and any fused column synthesized from it) becomes `NA`; its state stays `D`, matching the documented rule that states never inspect lengths.
