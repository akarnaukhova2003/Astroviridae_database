import argparse

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment


def split_alignment(input_file, output_file, start, end):
    """
    Extract a region from a multiple sequence alignment.

    Coordinates are 1-based and inclusive:
        start=1, end=3060
    corresponds to Python slicing:
        [0:3060]
    """

    alignment = AlignIO.read(input_file, "fasta")

    part = MultipleSeqAlignment(
        record[start - 1:end]
        for record in alignment
    )

    AlignIO.write(
        part,
        output_file,
        "fasta"
    )

    print(
        f"Saved: {output_file} "
        f"({start}-{end}, length={end - start + 1})"
    )


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="Extract a region from a nucleotide alignment."
    )

    parser.add_argument(
        "--input",
        required=True,
        help="Input FASTA alignment"
    )

    parser.add_argument(
        "--output",
        required=True,
        help="Output FASTA alignment"
    )

    parser.add_argument(
        "--start",
        type=int,
        required=True,
        help="Start coordinate, 1-based"
    )

    parser.add_argument(
        "--end",
        type=int,
        required=True,
        help="End coordinate, 1-based and inclusive"
    )

    args = parser.parse_args()

    split_alignment(
        input_file=args.input,
        output_file=args.output,
        start=args.start,
        end=args.end
    )
