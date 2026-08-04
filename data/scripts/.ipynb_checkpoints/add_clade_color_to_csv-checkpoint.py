#!/usr/bin/env python3

import argparse
import pandas as pd


def read_ids(file):
    """
    Читает файл с ID (один ID на строку)
    """

    with open(file) as f:
        return {
            line.strip()
            for line in f
            if line.strip()
        }



def assign_colors(id_files, colors):

    id_to_color = {}

    for file, color in zip(id_files, colors):

        ids = read_ids(file)

        print(f"{file}: {len(ids)} IDs")


        for id_ in ids:

            # если ID встретился только в одном кладе
            if id_ not in id_to_color:

                id_to_color[id_] = color

            else:
                print(
                    f"Warning: {id_} appears in multiple files"
                )


    return id_to_color



def add_color_column(
        input_csv,
        id_column,
        id_to_color,
        output_csv):


    df = pd.read_csv(input_csv)


    df["clade_color"] = (
        df[id_column]
        .map(id_to_color)
        .fillna("NA")
    )


    df.to_csv(
        output_csv,
        index=False
    )


    print(
        f"Saved: {output_csv}"
    )



def main():

    parser = argparse.ArgumentParser(
        description=
        "Add clade_color based on ID lists"
    )


    parser.add_argument(
        "--input",
        required=True,
        help="Input CSV file"
    )


    parser.add_argument(
        "--id_column",
        required=True,
        help="Column with accession IDs"
    )


    parser.add_argument(
        "--id_files",
        nargs=5,
        required=True,
        help="Five files with IDs"
    )


    parser.add_argument(
        "--colors",
        nargs=5,
        required=True,
        help=
        "Five clade labels, e.g. B G R PY O"
    )


    parser.add_argument(
        "--output",
        required=True,
        help="Output CSV"
    )


    args = parser.parse_args()


    id_to_color = assign_colors(
        args.id_files,
        args.colors
    )


    add_color_column(
        args.input,
        args.id_column,
        id_to_color,
        args.output
    )



if __name__ == "__main__":
    main()