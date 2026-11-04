"""Apply two branch-length draws to the same unrooted edge."""
from splitaligner_fast import Library

lib = Library()
parent = [-1, 0, 0]
taxon = [-1, 0, 1]
with lib.reference(2, "A,B", parent, taxon) as ref:
    with ref.coverage([0, 1]) as cov:
        with cov.mapping(parent, taxon, gene_id="gene1") as mapping:
            for lengths in ([0, 0.25, 0.50], [0, 0.50, 1.00]):
                values, valid = mapping.apply(lengths)["primitive"]
                print(values[0] if valid[0] else "NA")
