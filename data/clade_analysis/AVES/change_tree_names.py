import argparse
import re
import os
from Bio import SeqIO

def change_tree_names(fasta_name, tree_name):

    with open(tree_name) as tree_file:
        tree_line = tree_file.readline()
    tree_file.close()


    seqs = SeqIO.to_dict(SeqIO.parse(fasta_name, 'fasta'))
    seq_names = list(seqs.keys())

    #k=0
    for seq in seq_names:
        #k += 1
        seq_ac = seq.split('/')[0]
        print(seq_ac)
        #if seq_ac == 'NC':
        #    seq_ac = '_'.join(seq.split('_')[:2])
        search_seq = re.search(seq_ac+r'[\w\n_\.\?\/-]+', tree_line)
        if search_seq.group() == seq:
            continue
        else:
            #print(search_seq.group(), seq)
            tree_line = re.sub(seq_ac+r'[\w\n_\.\?\/-]+', seq, tree_line)
        #if k ==1591:
        #    print(seq_ac)
        #    break
    tree_file_name_new = os.path.splitext(tree_name)[0] + '_new.nwk'
    print(tree_file_name_new)

    with open(tree_file_name_new, 'w') as tree_file:
        tree_file.write(tree_line)
    tree_file.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-in_fasta", "--input_fasta", type=str,
                        help="Input fasta file with updated names", required=True)
    parser.add_argument("-in_tree", "--input_file_tree", type=str,
                        help="Input tree file in nwk format", required=True)
    args = parser.parse_args()
    change_tree_names(args.input_fasta, args.input_file_tree)