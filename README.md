1. This is a tool of making dotplot graphs according to Sibelia's output file "blocks_coords.txt".
2. There are two versions of output graphs, the first one uses the sequences(contigs) data from "blocks_coords.txt", while the second one orders the sequences(contigs) before drawing the dotplot.
3. In each version of output graphs, there are two output images, the first one labels axes by the seq_ids, the second one labels axes by the seq_lengths.
4. So there are a total of 4 different graphs outputted.

command: python Sibelia_dotplot.py <blocks_coords.txt path> <reference contigs number> <output file> <# of contigs in reference genome>