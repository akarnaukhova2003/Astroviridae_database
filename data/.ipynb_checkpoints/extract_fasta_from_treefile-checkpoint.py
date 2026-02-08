#!/usr/bin/env python3

import argparse
from Bio import SeqIO


def extract_taxa_id(file_tree, file_clusters, target_color, out_file):
    import re
    import pandas as pd

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

    classes = []
    for el in taxa:
        if el[1] == target_color:
            classes.append(el[0])

    import pandas as pd
    df = pd.read_csv(file_clusters, sep="\t")
    res = []
    for i in range(len(df)):
        if df.iloc[i, 0][0:3] != 'NC_':
            name_cluster = df.iloc[i, 0].split('_')[0]
        else:
            name_cluster = df.iloc[i, 0].split('_')[0] + '_' + df.iloc[i, 0].split('_')[1]
        if name_cluster in classes:
            if df.iloc[i, 1][0:3] != 'NC_':
                same_id = df.iloc[i, 1].split('_')[0]
            else:
                same_id = df.iloc[i, 1].split('_')[0] + '_' + df.iloc[i, 1].split('_')[1]
            res.append(same_id)

    with open(out_file, "w") as f:
        for el in res:
            f.write(el + "\n")


def extract_to_fasta_from_fasta(fasta_file, id_list_file, output_fasta):

    with open(id_list_file) as f:
        id_list = set(line.strip() for line in f if line.strip())
    with open(output_fasta, "w") as out:
        for record in SeqIO.parse(fasta_file, "fasta"):
            name = record.id.split()[0]
            if name[0:3] != 'NC_':
                acc = name.split('_')[0]
            else:
                acc = name.split('_')[0] + '_' + name.split('_')[1]
            if acc in id_list:
                header = f">{record.id}"
                out.write(header + "\n")
                seq = str(record.seq)
                for i in range(0, len(seq), 60):
                    out.write(seq[i:i+60] + "\n")


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument("--tree", required=True)
    parser.add_argument("--clusters", required=True)
    parser.add_argument("--color", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--out_ids", required=True)
    parser.add_argument("--out_fasta", required=True)

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
