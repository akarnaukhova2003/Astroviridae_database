#!/usr/bin/env python3

import argparse
import pandas as pd
from Bio import SeqIO
import numpy as np


DOMTBL_COLUMNS = [
    'gene_id',
    'pfam',
    'length',
    'sequence_id',
    'strand',
    'alignment_length',
    'evalue',
    'bit_score',
    'identity',
    'query_start',
    'query_end',
    'subject_start',
    'subject_end',
    'coverage',
    'some_ratio',
    'gaps',
    'mismatches',
    'q_len',
    's_len',
    'frame',
    'alignment_identity',
    'score',
    'category',
    'description_1',
    'description_2',
    'description_3'
]


def parse_domtblout(domtblout_file: str):
    df = pd.read_csv(
        domtblout_file,
        comment="#",
        sep=r"\s+",
        names=DOMTBL_COLUMNS,
        engine="python",
        on_bad_lines="skip"
    )
    return df


def get_ids_with_rdrp(df: pd.DataFrame):
    filtered_df = df.loc[df['gene_id'].str.contains("RdRP", na=False)]
    all_ids = df['sequence_id'].dropna().unique()
    bad_ids = list(set(df['sequence_id']) - set(filtered_df['sequence_id']))
    prefix = pd.Series(
    np.where(
        df["sequence_id"].str.startswith("NC_"),
        df["sequence_id"].str.split("_").str[:2].str.join("_"),
        df["sequence_id"].str.split("_").str[0]
    ),
    index=df.index)
    print(bad_ids)
    filtered_df_1a = df[df['sequence_id'].isin(bad_ids)]
    return set(filtered_df['sequence_id'].dropna().unique()), bad_ids, filtered_df, filtered_df_1a


def filter_fasta(good_ids: set, input_fasta: str, output_fasta: str):
    with open(output_fasta, "w") as out_f:
        for record in SeqIO.parse(input_fasta, "fasta"):
            if record.id in good_ids:
                SeqIO.write(record, out_f, "fasta")


def main():
    parser = argparse.ArgumentParser(
        description="Remove sequences from FASTA based on Pfam domtblout (non-RdRP hits)"
    )

    parser.add_argument(
        "-d", "--domtblout",
        required=True,
        help="Pfam domtblout file"
    )
    parser.add_argument(
        "-i", "--input-fasta",
        required=True,
        help="Input FASTA file"
    )
    parser.add_argument(
        "-o", "--output-fasta",
        required=True,
        help="Output FASTA file (RdRP-only sequences kept)"
    )
    parser.add_argument(
        "--csv",
        help="Optional: save filtered domtblout to CSV"
    )
    parser.add_argument(
        "--ids",
        help="Optional: save list of removed sequence IDs"
    )

    args = parser.parse_args()

    df = parse_domtblout(args.domtblout)
    good_ids, bad_ids, filtered_df, filtered_df_1a = get_ids_with_rdrp(df)

    if args.csv:
        filtered_df_1a.to_csv(args.csv, index=False)

    if args.ids:
        with open(args.ids, "w") as f:
            for seq_id in sorted(bad_ids):
                f.write(seq_id + "\n")

    filter_fasta(
        good_ids=good_ids,
        input_fasta=args.input_fasta,
        output_fasta=args.output_fasta
    )


if __name__ == "__main__":
    main()
