# Read the branch table for six species and three genes

All three genes lack species C. In its absence, reference branches B7 and B9
both separate A+B from D+E+F. A gene tree can supply one length for that split.

Run this block in Bash from the SplitAligner repository root. It creates a
fresh temporary directory and prints the two output tables used below.

```bash
(
  set -eu
  repo="$PWD"
  run_dir=$(mktemp -d "${TMPDIR:-/tmp}/splitaligner-six.XXXXXXXX")
  cd "$run_dir"
  cat > species.nwk <<'EOF'
((A,B),C,(D,(E,F)));
EOF
  cat > free.nwk <<'EOF'
g1((A:0.1,B:0.1):0.31,D:0.1,(E:0.1,F:0.1):0.17);
g2((A:0.1,B:0.1):0.26,E:0.1,(D:0.1,F:0.1):0.09);
g3((A:0.1,D:0.1):0.12,B:0.1,(E:0.1,F:0.1):0.23);
EOF
  cat > fix.nwk <<'EOF'
g1((A:0.1,B:0.1):0.2,(D:0.1,(E:0.1,F:0.1):0.15):0.1);
g2((A:0.1,B:0.1):0.15,(D:0.1,(E:0.1,F:0.1):0.12):0.1);
g3((A:0.1,B:0.1):0.18,(D:0.1,(E:0.1,F:0.1):0.2):0.1);
EOF
  perl "$repo/SplitAligner.pl" --mode matrix --species species.nwk --gene free.nwk --label free
  perl "$repo/SplitAligner.pl" --mode matrix --species species.nwk --gene fix.nwk --label fix
  perl "$repo/SplitAligner.pl" --mode finalize --free free.matrix_with_fuse.txt \
    --fix fix.matrix_with_fuse.txt --final_label final --species_tree species_tree.forSplit.nwk
  cat final.free.na_classified.txt final.free.projected_groups.tsv
  printf '\nAll outputs: %s\n' "$run_dir"
)
```

In `final.free.na_classified.txt`, each column is an original reference branch.
The table below shows four of its branch columns.

| Gene | B3 (C) | B7 | B8 (E+F) | B9 |
|---|---|---|---|---|
| g1 | NA_struct | NA_fuse | 0.17 | NA_fuse |
| g2 | NA_struct | NA_fuse | NA_topo | NA_fuse |
| g3 | NA_struct | NA_fuse | 0.23 | NA_fuse |

B3 is structurally absent because C is missing. The free tree for g2 places D
next to F, so the E+F split at B8 is absent. B7 and B9 keep `NA_fuse` in
every gene. The shared split has one row per gene in `final.free.projected_groups.tsv`.

| Gene | Members | Group status | Shared length |
|---|---|---|---|
| g1 | B7 and B9 | recovered | 0.31 |
| g2 | B7 and B9 | recovered | 0.26 |
| g3 | B7 and B9 | NA_fuse_topo | NA |

In g3 the shared split is absent. The group gets one `NA_fuse_topo` record,
and the two original cells still read `NA_fuse`.

Across the **original branch columns** there are 3 genes × 9 branches = 27
cells. **27 = 17 Mapped + 3 NA_struct + 6 NA_fuse + 1 NA_topo**.
The six fusion-member cells are B7 and B9 in each of the three genes.

In `final.support_b.txt`, all three fixed trees supply a length for B8, and
two free trees do. Its support is therefore 2/3 (66.7%). B3, B7 and B9 each
have a zero denominator because no gene can be evaluated individually for that branch.
Perl prints zero in the support column for those rows.
