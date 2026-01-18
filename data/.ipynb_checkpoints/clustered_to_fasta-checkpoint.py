#!/usr/bin/env python3

import argparse
from Bio import SeqIO


def extract_sequences_by_id(
    file_with_names: str,
    file_with_sequences: str,
    output_file: str
):
    names_to_extract = set()

    # читаем ID из первого FASTA
    with open(file_with_names) as f:
        for record in SeqIO.parse(f, "fasta"):
            names_to_extract.add(record.id)

    # фильтруем последовательности
    with open(output_file, "w") as out_f:
        for record in SeqIO.parse(file_with_sequences, "fasta"):
            if record.id in names_to_extract:
                SeqIO.write(record, out_f, "fasta")


def main():
    parser = argparse.ArgumentParser(
        description="Extract sequences from FASTA by IDs listed in another FASTA"
    )
    parser.add_argument(
        "-n", "--names",
        required=True,
        help="FASTA file with sequence IDs to extract"
    )
    parser.add_argument(
        "-s", "--sequences",
        required=True,
        help="FASTA file with sequences"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output FASTA file"
    )

    args = parser.parse_args()

    extract_sequences_by_id(
        file_with_names=args.names,
        file_with_sequences=args.sequences,
        output_file=args.output
    )


if __name__ == "__main__":
    main()

