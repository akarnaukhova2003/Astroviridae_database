#!/usr/bin/env python3

import argparse
import re
import pandas as pd
from Bio import SeqIO


def extract_taxa_id(file_tree, file_clusters, target_colors, out_file):

    taxa = []

    with open(file_tree) as f:
        text = f.read()

    match = re.search(
        r"taxlabels(.*?);",
        text,
        re.S | re.I
    )

    if not match:
        raise ValueError("Cannot find taxlabels block in tree file")

    block = match.group(1)

    for line in block.splitlines():
        line = line.strip()
        if not line:
            continue

        # берём accession из NEXUS до "/"
        name = re.split(r"/", line)[0]
        name = name.replace("'", "")

        # ищем цвет FigTree
        m = re.search(
            r"color\s*=\s*(#[0-9A-Fa-f]{6})",
            line
        )
        tax_color = m.group(1) if m else None

        taxa.append((name, tax_color))

    # ID всех выбранных цветов
    classes = [
        name
        for name, color in taxa
        if color in target_colors
    ]

    res = []

    # расширяем через кластеры
    if file_clusters:
        df = pd.read_csv(
            file_clusters,
            sep="\t",
            header=None
        )

        for _, row in df.iterrows():
            # Извлекаем Accession из первого столбца до "/"
            cluster_str = str(row.iloc[0])
            cluster_id = cluster_str.split("/")[0] if "/" in cluster_str else cluster_str

            if cluster_id in classes:
                seq_str = str(row.iloc[1])
                seq_id = seq_str.split("/")[0] if "/" in seq_str else seq_str
                res.append(seq_id)

    # если кластеров нет
    else:
        res = classes.copy()

    # специальный случай
    if "MF033385" in classes:
        print("Adding MF033385 extra IDs")
        extra_ids = [
            "MF033386",
            "KT946732",
            "KT946733",
            "KT946734",
            "MT549858",
            "MT549859",
            "PQ678009",
            "MF033385"
        ]
        res.extend(extra_ids)

    # удаляем дубли
    res = list(dict.fromkeys(res))

    with open(out_file, "w") as f:
        for el in res:
            f.write(el + "\n")

    print(f"IDs extracted from tree: {len(res)}")


def extract_to_fasta_from_fasta(
        fasta_file,
        id_list_file,
        output_fasta,
        separator):

    with open(id_list_file) as f:
        id_list = {
            line.strip()
            for line in f
            if line.strip()
        }

    count = 0

    with open(output_fasta, "w") as out:
        for record in SeqIO.parse(fasta_file, "fasta"):
            # извлекаем accession
            if separator:
                acc = record.id.split(separator)[0]
            else:
                acc = record.id

            if acc in id_list:
                count += 1
                out.write(f">{record.id}\n")
                seq = str(record.seq)
                for i in range(0, len(seq), 60):
                    out.write(seq[i:i+60] + "\n")

    print(f"Sequences extracted from FASTA: {count}")


def main():
    parser = argparse.ArgumentParser(
        description="Extract sequences by colors from NEXUS tree and FASTA"
    )

    parser.add_argument(
        "--tree",
        required=True,
        help="NEXUS tree file"
    )

    parser.add_argument(
        "--clusters",
        required=False,
        default=None,
        help="Clusters TSV file (optional)"
    )

    parser.add_argument(
        "--colors",
        nargs="+",
        required=True,
        help="Colors to extract, e.g. #3399ff #9933ff"
    )

    parser.add_argument(
        "--fasta",
        required=True,
        help="Input FASTA file"
    )

    parser.add_argument(
        "--separator",
        required=False,
        default="/",
        help="Separator in FASTA header for accession extraction, e.g. / _ ."
    )

    parser.add_argument(
        "--out_ids",
        required=True,
        help="Output IDs file"
    )

    parser.add_argument(
        "--out_fasta",
        required=True,
        help="Output FASTA file"
    )

    args = parser.parse_args()

    extract_taxa_id(
        file_tree=args.tree,
        file_clusters=args.clusters,
        target_colors=args.colors,
        out_file=args.out_ids
    )

    extract_to_fasta_from_fasta(
        fasta_file=args.fasta,
        id_list_file=args.out_ids,
        output_fasta=args.out_fasta,
        separator=args.separator
    )


if __name__ == "__main__":
    main()