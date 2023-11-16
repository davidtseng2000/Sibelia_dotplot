1. This tool creates dotplot graphs based on Sibelia's output file 'blocks_coords.txt'.
2. There are two versions of the output graphs: the first uses the sequences (contigs) data from 'blocks_coords.txt', while the second orders the sequences (contigs) before drawing the dotplots.
3. In each version of the output graphs, there are two output images: the first labels axes by the seq_ids, and the second labels axes by the seq_lengths.
4. Therefore, a total of 4 different graphs are outputted.

command: python Sibelia_dotplot.py <blocks_coords.txt path> <reference contigs number> <output file> <# of contigs in reference genome>