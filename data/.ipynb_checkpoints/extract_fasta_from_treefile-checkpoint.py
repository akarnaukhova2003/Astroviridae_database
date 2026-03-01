#!/usr/bin/env python3

import argparse
import re
import pandas as pd
from Bio import SeqIO


def extract_taxa_id(file_tree, file_clusters, target_color, out_file):

    taxa = []

    with open(file_tree) as f:
        text = f.read()

    block = re.search(r"taxlabels(.*?);", text, re.S | re.I).group(1)

    for line in block.splitlines():
        line = line.strip()
        if not line:
            continue

        name = re.split(r"/", line)[0].strip()
        m = re.search(r"color\s*=\s*(#[0-9A-Fa-f]{6})", line)
        tax_color = m.group(1) if m else None

        taxa.append((name[1:], tax_color))

    classes = [name for name, color in taxa if color == target_color]

    df = pd.read_csv(file_clusters, sep="\t")

    res = []

    for i in range(len(df)):
        name_cluster = df.iloc[i, 0].split('_')[0]
        if name_cluster in classes:
            same_id = df.iloc[i, 1].split('_')[0]
            res.append(same_id)


    #частный случай для подмены лидера кластера 
    if 'MF033385' in classes: 
        print('yeees')
        extra_ids = [
            'MF033386',
            'KT946732',
            'KT946733',
            'KT946734',
            'MT549858',
            'MT549859',
            'PQ678009',
            'MF033385'
        ]
        res.extend(extra_ids)

    res = list(dict.fromkeys(res))

    with open(out_file, "w") as f:
        for el in res:
            f.write(el + "\n")


def extract_to_fasta_from_fasta(fasta_file, id_list_file, output_fasta):

    with open(id_list_file) as f:
        id_list = set(line.strip() for line in f if line.strip())

    with open(output_fasta, "w") as out:
        for record in SeqIO.parse(fasta_file, "fasta"):
            acc = record.id.split('_')[0]
            if acc in id_list:
                out.write(f">{record.id}\n")
                seq = str(record.seq)
                for i in range(0, len(seq), 60):
                    out.write(seq[i:i+60] + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Extract sequences by cluster color from tree and FASTA"
    )

    parser.add_argument("--tree", required=True, help="Tree file (with taxlabels)")
    parser.add_argument("--clusters", required=True, help="Clusters TSV file")
    parser.add_argument("--color", required=True, help="Target color, e.g. #ff0000")
    parser.add_argument("--fasta", required=True, help="Input FASTA")
    parser.add_argument("--out_ids", required=True, help="Output file with IDs")
    parser.add_argument("--out_fasta", required=True, help="Output FASTA")

    args = parser.parse_args()

    extract_taxa_id(
        args.tree,
        args.clusters,
        args.color,
        args.out_ids
    )

    extract_to_fasta_from_fasta(
        args.fasta,
        args.out_ids,
        args.out_fasta
    )


if __name__ == "__main__":
    main()