#!/usr/bin/env python3

import argparse
import pandas as pd
import re


def normalize_id(id_):
    """
    Нормализация accession ID:
    NC-123 == NC_123 == NC123
    """

    if pd.isna(id_):
        return None

    id_ = str(id_).strip()

    # убираем разделители
    id_ = re.sub(r"[-_]", "", id_)

    return id_



def read_ids(file):
    """
    Читает файл с ID (один ID на строку)
    """

    ids = set()

    with open(file) as f:
        for line in f:
            line = line.strip()

            if line:
                ids.add(
                    normalize_id(line)
                )

    return ids



def assign_colors(id_files, colors):

    id_to_color = {}

    for file, color in zip(id_files, colors):

        ids = read_ids(file)

        print(f"{file}: {len(ids)} IDs")


        for id_ in ids:

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


    # если CSV разделён ;
    df = pd.read_csv(
        input_csv,
        sep=";"
    )


    # создаём нормализованный ID
    df["_normalized_id"] = (
        df[id_column]
        .apply(normalize_id)
    )


    df["clade_color"] = (
        df["_normalized_id"]
        .map(id_to_color)
        .fillna("NA")
    )


    # удаляем временный столбец
    df.drop(
        columns=["_normalized_id"],
        inplace=True
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
        nargs="+",
        required=True,
        help="Files with IDs"
    )


    parser.add_argument(
        "--colors",
        nargs="+",
        required=True,
        help="Color labels corresponding to ID files"
    )


    parser.add_argument(
        "--output",
        required=True,
        help="Output CSV"
    )


    args = parser.parse_args()


    if len(args.id_files) != len(args.colors):
        raise ValueError(
            "Number of id_files must equal number of colors"
        )


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