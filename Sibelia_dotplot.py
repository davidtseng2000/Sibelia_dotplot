# command: python Sibelia_dotplot.py <blocks_coords.txt path> <reference contigs number> <output file>

import os
# os.environ['MPLCONFIGDIR'] = os.getcwd() + "/configs"
import sys
import matplotlib.pyplot as plt


# extract_mapping is to extract seqs from blocks_coords.txt
# Function to extract mapping from blocks_coords
def extract_mapping(file_path): 
    name_id_map = {}
    id_size_map = {}
    with open(file_path, 'r') as f:
        for line in f:
            if not line.startswith("Seq_id"):
                if line.startswith('-'):
                    break
                seq_id, size , gene_name = line.strip().split()
                name_id_map[gene_name] = seq_id
                id_size_map[seq_id] = int(size)
    return name_id_map, id_size_map


# extract_mapping_sort is to extract seqs from blocks_coords.txt and order the seqs by their lengths
order_ref_sort = [0]
order_target_sort = [0]
id_order_map_tar_sort = {}
id_order_map_ref_sort = {}
sorted_tar = [] 

def extract_mapping_sort(file_path):

    global sorted_tar
    global ref_num

    line_idx = 1
    id_size_map = {}

    tar_idx = 1
    ref_idx = 1

    all_tar = [(sys.maxsize, sys.maxsize)]


    with open(file_path, 'r') as f:
        for line in f:
            if line.strip() == '':
                break

            if not line.startswith("Seq_id"):
                if line.startswith('-'):
                    break
                seq_id, size , gene_name = line.strip().split()
                seq_id = int(seq_id)
                size = int(size)

                if(line_idx <= ref_num):
                    id_size_map[seq_id] = size
                    order_ref_sort.append(seq_id)
                    id_order_map_ref_sort[seq_id] = ref_idx
                    ref_idx = ref_idx + 1
                else:
                    all_tar.append((seq_id, size))
            
                line_idx = line_idx + 1

        sorted_tar = sorted(all_tar, key=lambda x: x[1], reverse=True)
        sorted_tar = sorted_tar[1:] # discard the dummy one

        for ele in sorted_tar:
            id_size_map[ele[0]] = ele[1]
            order_target_sort.append(ele[0])
            id_order_map_tar_sort[ele[0]] = tar_idx
            tar_idx = tar_idx + 1

    return  id_size_map




# Function to extract data based on filter
# ref_num = 0
def extract_data(file_path, mapping, blocks_coords, tar_or_ref):
    global ref_num
    results = []
    with open(file_path, 'r') as f:
        gene_name = ''
        for line in f:
            if line.startswith('>'):
                gene_name = line.strip()[1:]
                if(tar_or_ref == 1):
                    ref_num = ref_num + 1
            else:
                block_num, seq_num = map(int, line.strip().split())
                block_idx = 0
                while blocks_coords[block_num][block_idx][0] != mapping[gene_name]:
                    block_idx = block_idx + 1

                block_data = blocks_coords[block_num][block_idx + seq_num - 1]
                if block_data[1] == '-':
                    results.append((mapping[gene_name], block_data[3], block_data[2]))                    
                else:
                    results.append((mapping[gene_name], block_data[2], block_data[3]))
    return results



def refine_blocks_coords(file_path):
    blocks = {}
    start_reading = False
    with open(file_path, 'r') as f:
        current_block = None
        for line in f:
            if line.startswith("-"):
                start_reading = True
                continue
            if start_reading:
                if line.startswith("Block"):
                    current_block = int(line.split('#')[1].strip())
                    blocks[current_block] = []
                elif line.strip() and not line.startswith('Seq_id'):
                    if((line.strip().split())[1] == '+'):
                        blocks[current_block].append(((line.strip().split())[0], (line.strip().split())[2], (line.strip().split())[3], (line.strip().split())[1]))
                    else:
                        blocks[current_block].append(((line.strip().split())[0], (line.strip().split())[3], (line.strip().split())[2], (line.strip().split())[1]))
    return blocks


# Main function
ref_num = 0
def main():

    global ref_num

    if len(sys.argv) != 4:
        print("Usage: python Sibelia_dotplot.py <blocks_coords.txt path> <reference contigs number> <output file>")
        sys.exit()
    
    blocks_coords_path = sys.argv[1]
    ref_num = int(sys.argv[2])
    outputfile = sys.argv[3]

    
    #############################
    # Loading some necessary data
    #############################
    name_id_map, id_size_map = extract_mapping(blocks_coords_path)
    
    #####################################################
    # Draw dotplot (original data from blocks_coords.txt)
    #####################################################
    refine_blocks = refine_blocks_coords(blocks_coords_path)
    name_id_map, id_size_map = extract_mapping(blocks_coords_path)    

    seq_ids = [0]
    sizes = [0]
    for id, size in id_size_map.items():
        seq_ids.append(id)
        sizes.append(size)

    new_ids_tar = [0] 
    id_nickname_map_tar = {} #contig nickname for target contigs
    for id in seq_ids[ref_num+1:]:
        new_ids_tar.append(int(id) - ref_num)
        id_nickname_map_tar[int(id)] = int(id) - ref_num
    new_ids_tar.append(" ")

    new_ids_ref = [0] #contig nickname for ref contigs including last " "
    id_nickname_map_ref = {} #contig nickname for ref contigs
    for id in seq_ids[1:ref_num+1]:
        new_ids_ref.append(int(id))
        id_nickname_map_ref[int(id)] = int(id)
    new_ids_ref.append(" ")

    # y-axis (target contigs)
    y_positions = [0]  # starting point
    for idx in range(ref_num+1, len(sizes)):
        size = sizes[idx]
        y_positions.append(y_positions[-1] + int(size))
    
    # x-axis (ref contigs)
    x_positions = [0]  # starting point
    for idx in range(1, ref_num+1):
        size = sizes[idx]
        x_positions.append(x_positions[-1] + int(size))
        
    # (1.) draw dotplot with seq_id
    plt.figure(figsize=(9, 10))
    plt.xlabel('contigs in reference genome (name)')
    plt.ylabel('contigs in target genome (name)')   
    plt.title('dotplot')

    plt.xticks(x_positions, new_ids_ref[1:], rotation=90)
    plt.yticks(y_positions, new_ids_tar[1:], rotation=90)
    for x in x_positions:
        plt.axvline(x=x, color=(0.8, 0.8, 0.8), linestyle='--', zorder=1)
    for y in y_positions:
        plt.axhline(y=y, color=(0.8, 0.8, 0.8), linestyle='--', zorder=1)


    for i in range(1, len(refine_blocks)):    
        for j in range(len(refine_blocks[i])):
            for k in range(j+1, len(refine_blocks[i])):

                # if the marker is not from reference-target 
                if not (int(refine_blocks[i][j][0]) in range(1, ref_num+1) and not (int(refine_blocks[i][k][0]) in range(1, ref_num+1))):
                    continue

                if (refine_blocks[i][k][3] == '+' and refine_blocks[i][j][3] == '+') or (refine_blocks[i][k][3] == '-' and refine_blocks[i][j][3] == '-'):
                    y_coords = [int(refine_blocks[i][k][1])+y_positions[int(refine_blocks[i][k][0])-2], int(refine_blocks[i][k][2])+y_positions[int(refine_blocks[i][k][0])-2]]
                else:
                    y_coords = [int(refine_blocks[i][k][2])+y_positions[int(refine_blocks[i][k][0])-2], int(refine_blocks[i][k][1])+y_positions[int(refine_blocks[i][k][0])-2]]
                
                x_coords = [int(refine_blocks[i][j][1])+x_positions[int(refine_blocks[i][j][0])-1], int(refine_blocks[i][j][2])+x_positions[int(refine_blocks[i][j][0])-1]]                 
                
                if (refine_blocks[i][k][3] == '+' and refine_blocks[i][j][3] == '+') or (refine_blocks[i][k][3] == '-' and refine_blocks[i][j][3] == '-'):
                    plt.scatter(x=x_coords, y=y_coords, facecolors='none', edgecolors='red', zorder=2)
                    plt.plot(x_coords, y_coords, color='red', marker=None, zorder=2)
                else:
                    plt.scatter(x=x_coords, y=y_coords, facecolors='none', edgecolors='blue', zorder=2)
                    plt.plot(x_coords, y_coords, color='blue', marker=None, zorder=2)
    
    ax = plt.gca()
    ax.spines['bottom'].set_position(('data',0))
    ax.spines['left'].set_position(('data',0))  
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    
    plt.tight_layout()
    plt.savefig(str(outputfile) + "/seq_name_origin.svg")
    plt.savefig(str(outputfile) + "/seq_name_origin.png")
    plt.show()

    plt.clf()
    plt.cla()
    plt.close()

    # (2.) draw dotplot without seq_id (but with length)
    plt.figure(figsize=(9, 10))
    plt.xlabel('reference genome (length)')
    plt.ylabel('target genome (length)')    
    plt.title('dotplot')

    for x in x_positions:
        plt.axvline(x=x, color=(0.8, 0.8, 0.8), linestyle='--', zorder=1)
    for y in y_positions:
        plt.axhline(y=y, color=(0.8, 0.8, 0.8), linestyle='--', zorder=1)


    for i in range(1, len(refine_blocks)):    
        for j in range(len(refine_blocks[i])):
            for k in range(j+1, len(refine_blocks[i])):

                if not (int(refine_blocks[i][j][0]) in range(1, ref_num+1) and not (int(refine_blocks[i][k][0]) in range(1, ref_num+1))):
                    continue

                if (refine_blocks[i][k][3] == '+' and refine_blocks[i][j][3] == '+') or (refine_blocks[i][k][3] == '-' and refine_blocks[i][j][3] == '-'):
                    y_coords = [int(refine_blocks[i][k][1])+y_positions[int(refine_blocks[i][k][0])-2], int(refine_blocks[i][k][2])+y_positions[int(refine_blocks[i][k][0])-2]]
                else:
                    y_coords = [int(refine_blocks[i][k][2])+y_positions[int(refine_blocks[i][k][0])-2], int(refine_blocks[i][k][1])+y_positions[int(refine_blocks[i][k][0])-2]]
                
                x_coords = [int(refine_blocks[i][j][1])+x_positions[int(refine_blocks[i][j][0])-1], int(refine_blocks[i][j][2])+x_positions[int(refine_blocks[i][j][0])-1]]
                            
                if (refine_blocks[i][k][3] == '+' and refine_blocks[i][j][3] == '+') or (refine_blocks[i][k][3] == '-' and refine_blocks[i][j][3] == '-'):
                    plt.scatter(x=x_coords, y=y_coords, facecolors='none', edgecolors='red', zorder=2)
                    plt.plot(x_coords, y_coords, color='red', marker=None, zorder=2)
                else:
                    plt.scatter(x=x_coords, y=y_coords, facecolors='none', edgecolors='blue', zorder=2)
                    plt.plot(x_coords, y_coords, color='blue', marker=None, zorder=2)

    ax = plt.gca()
    ax.spines['bottom'].set_position(('data',0))
    ax.spines['left'].set_position(('data',0))  
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)

    plt.tight_layout()
    plt.savefig(str(outputfile) + "/seq_length_origin.svg")
    plt.savefig(str(outputfile) + "/seq_length_origin.png")
    plt.show()

    plt.clf()
    plt.cla()
    plt.close()

    #########################################################
    # draw dotplot where contigs are ordered by their lengths
    #########################################################

    id_size_map = extract_mapping_sort(blocks_coords_path)

    refine_blocks = refine_blocks_coords(blocks_coords_path)

    seq_ids = [0]
    sizes = [0]
    for id, size in id_size_map.items():
        seq_ids.append(id)
        sizes.append(size)


    new_ids_tar = [0] #contig nickname for target contigs including last " "
    for id in order_target_sort[1:]:
        new_ids_tar.append(id_nickname_map_tar[int(id)])
    new_ids_tar.append(" ")

    new_ids_ref = [0] #contig nickname for ref contigs including last " "
    for id in order_ref_sort[1:]:
        new_ids_ref.append(id_nickname_map_ref[int(id)])
    new_ids_ref.append(" ")


    y_positions = [0]   # starting point
    for idx in range(ref_num+1, len(sizes)):
        size = sizes[idx]
        y_positions.append(y_positions[-1] + int(size))
    
    x_positions = [0]  # starting point
    for idx in range(1, ref_num+1):
        size = sizes[idx]
        x_positions.append(x_positions[-1] + int(size))
        
    # (1.) draw dotplot with seq_id
    plt.figure(figsize=(9, 10))
    plt.xlabel('contigs in reference genome (name)')
    plt.ylabel('contigs in target genome (name)')    
    plt.title('dotplot')

    plt.xticks(x_positions, new_ids_ref[1:], rotation=90)
    plt.yticks(y_positions, new_ids_tar[1:], rotation=90)
    for x in x_positions:
        plt.axvline(x=x, color=(0.8, 0.8, 0.8), linestyle='--', zorder=1)
    for y in y_positions:
        plt.axhline(y=y, color=(0.8, 0.8, 0.8), linestyle='--', zorder=1)

    for i in range(1, len(refine_blocks)):  
        for j in range(len(refine_blocks[i])):
            for k in range(j+1, len(refine_blocks[i])):

                if not (int(refine_blocks[i][j][0]) in range(1, ref_num+1) and not (int(refine_blocks[i][k][0]) in range(1, ref_num+1))):
                    continue
                
                if (refine_blocks[i][k][3] == '+' and refine_blocks[i][j][3] == '+') or (refine_blocks[i][k][3] == '-' and refine_blocks[i][j][3] == '-'):
                    y_coords = [int(refine_blocks[i][k][1])+y_positions[id_order_map_tar_sort[int(refine_blocks[i][k][0])]-1], int(refine_blocks[i][k][2])+y_positions[id_order_map_tar_sort[int(refine_blocks[i][k][0])]-1]]
                else:
                    y_coords = [int(refine_blocks[i][k][2])+y_positions[id_order_map_tar_sort[int(refine_blocks[i][k][0])]-1], int(refine_blocks[i][k][1])+y_positions[id_order_map_tar_sort[int(refine_blocks[i][k][0])]-1]]
                
                x_coords = [int(refine_blocks[i][j][1])+x_positions[id_order_map_ref_sort[int(refine_blocks[i][j][0])]-1], int(refine_blocks[i][j][2])+x_positions[id_order_map_ref_sort[int(refine_blocks[i][j][0])]-1]]
                                
                if (refine_blocks[i][k][3] == '+' and refine_blocks[i][j][3] == '+') or (refine_blocks[i][k][3] == '-' and refine_blocks[i][j][3] == '-'):
                    plt.scatter(x=x_coords, y=y_coords, facecolors='none', edgecolors='red', zorder=2)
                    plt.plot(x_coords, y_coords, color='red', marker=None, zorder=2)
                else:
                    plt.scatter(x=x_coords, y=y_coords, facecolors='none', edgecolors='blue', zorder=2)
                    plt.plot(x_coords, y_coords, color='blue', marker=None, zorder=2)

    ax = plt.gca()
    ax.spines['bottom'].set_position(('data',0))
    ax.spines['left'].set_position(('data',0))
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    
    plt.tight_layout()
    plt.savefig(str(outputfile) + "/seq_name_ordered.svg")
    plt.savefig(str(outputfile) + "/seq_name_ordered.png")
    plt.show()

    plt.clf()
    plt.cla()
    plt.close()

    # (2.) draw dotplot without seq_id (but with length)
    plt.figure(figsize=(9, 10))
    plt.xlabel('reference genome (length)')
    plt.ylabel('target genome (length)')   
    plt.title('dotplot')

    for x in x_positions: 
        plt.axvline(x=x, color=(0.8, 0.8, 0.8), linestyle='--', zorder=1)
    for y in y_positions:
        plt.axhline(y=y, color=(0.8, 0.8, 0.8), linestyle='--', zorder=1)

    for i in range(1, len(refine_blocks)):  
        for j in range(len(refine_blocks[i])):
            for k in range(j+1, len(refine_blocks[i])):

                if not (int(refine_blocks[i][j][0]) in range(1, ref_num+1) and not (int(refine_blocks[i][k][0]) in range(1, ref_num+1))):
                    continue

                if (refine_blocks[i][k][3] == '+' and refine_blocks[i][j][3] == '+') or (refine_blocks[i][k][3] == '-' and refine_blocks[i][j][3] == '-'):
                    y_coords = [int(refine_blocks[i][k][1])+y_positions[id_order_map_tar_sort[int(refine_blocks[i][k][0])]-1], int(refine_blocks[i][k][2])+y_positions[id_order_map_tar_sort[int(refine_blocks[i][k][0])]-1]]
                else:
                    y_coords = [int(refine_blocks[i][k][2])+y_positions[id_order_map_tar_sort[int(refine_blocks[i][k][0])]-1], int(refine_blocks[i][k][1])+y_positions[id_order_map_tar_sort[int(refine_blocks[i][k][0])]-1]]
                
                x_coords = [int(refine_blocks[i][j][1])+x_positions[id_order_map_ref_sort[int(refine_blocks[i][j][0])]-1], int(refine_blocks[i][j][2])+x_positions[id_order_map_ref_sort[int(refine_blocks[i][j][0])]-1]]
                                
                if (refine_blocks[i][k][3] == '+' and refine_blocks[i][j][3] == '+') or (refine_blocks[i][k][3] == '-' and refine_blocks[i][j][3] == '-'):
                    plt.scatter(x=x_coords, y=y_coords, facecolors='none', edgecolors='red', zorder=2)
                    plt.plot(x_coords, y_coords, color='red', marker=None, zorder=2)
                else:
                    plt.scatter(x=x_coords, y=y_coords, facecolors='none', edgecolors='blue', zorder=2)
                    plt.plot(x_coords, y_coords, color='blue', marker=None, zorder=2)

    ax = plt.gca()
    ax.spines['bottom'].set_position(('data',0)) 
    ax.spines['left'].set_position(('data',0))
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    
    plt.tight_layout()
    plt.savefig(str(outputfile) + "/seq_length_ordered.svg")
    plt.savefig(str(outputfile) + "/seq_length_ordered.png")
    plt.show()

    plt.clf()
    plt.cla()
    plt.close()

    # Ending
    print("Sibelia_dotplot.py done!")
    

if __name__ == '__main__':
    main()