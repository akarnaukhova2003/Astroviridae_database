#!/usr/bin/env python3

import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Merge Astroviridae tables and add host class information"
    )

    parser.add_argument(
        "--lineage",
        required=True,
        help="CSV file with lineage rank information (e.g. Astroviridae_lineage_rank.csv)"
    )
    parser.add_argument(
        "--results",
        required=True,
        help="CSV file with main results table (e.g. Astroviridae_15102025.csv)"
    )
    parser.add_argument(
        "--coords",
        required=True,
        help="CSV file with ORF coordinates (e.g. Astroviridae_26112025_orf-coords.csv)"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output CSV file name"
    )

    args = parser.parse_args()

    df = pd.read_csv(args.lineage)
    df_res = pd.read_csv(args.results)
    df_coords = pd.read_csv(args.coords)

    result = df[["host", "class"]]
    df_res["class"] = None

    for i in range(len(df_res)):
        host1 = df_res.loc[i, "Host"]
        for j in range(len(result)):
            host2 = result.loc[j, "host"]
            if host1 == host2:
                df_res.loc[i, "class"] = result.loc[j, "class"]
                break

    df_res.rename(columns={df_res.columns[0]: "ID"}, inplace=True)

    df_res = df_res.merge(
        df_coords,
        how="left",
        left_on="ID",
        right_on="GBAC"
    )

    df_res.drop(columns=["GBAC"], inplace=True)

    df_res.to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
