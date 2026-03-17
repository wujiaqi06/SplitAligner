export PATH=$PATH:SplitAligner
cd examples/302mammal
SplitAligner.pl --mode matrix --species input/speciesTree302.nwk --gene input/free_tree.examples.nwk --label free
SplitAligner.pl --mode matrix --species input/speciesTree302.nwk --gene input/fix_tree.examples.nwk --label fix
SplitAligner.pl --mode finalize --free free.matrix_with_fuse.txt --fix fix.matrix_with_fuse.txt --final_label final