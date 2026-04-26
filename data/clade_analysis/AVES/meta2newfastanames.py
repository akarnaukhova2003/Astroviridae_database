import argparse
import os
import re
from Bio import SeqIO
import pandas as pd


def create_mapping_dict(mapping_file_path):
    mapping_df = pd.read_csv(mapping_file_path, sep='\t', header=None,
                             names=['regex_pattern', 'short_code'])

    compiled_patterns = []
    for _, row in mapping_df.iterrows():
        try:
            pattern = re.compile(str(row['regex_pattern']), re.IGNORECASE)
            compiled_patterns.append((pattern, row['short_code']))
        except re.error as e:
            print(f"Warning: Invalid regex pattern '{row['regex_pattern']}': {e}")

    return compiled_patterns


def map_value_to_short_code(value, compiled_patterns, default="NA"):
    if pd.isna(value):
        return default

    for pattern, short_code in compiled_patterns:
        if pattern.search(str(value)):
            return short_code

    print(f"{value} was not found in mapping file")
    return str(value).replace(" ", "-")


def process_column(meta_df, mapping_file_path, column):
    compiled_patterns = create_mapping_dict(mapping_file_path)

    meta_df[column + "_original"] = meta_df[column]
    meta_df[column] = meta_df[column].apply(
        lambda x: map_value_to_short_code(x, compiled_patterns)
    )

    return meta_df


def extract_year_from_date(date_string):
    if pd.isna(date_string):
        return "NA"

    matches = re.findall(r'\b\d{4}\b', str(date_string))
    for match in matches:
        year = int(match)
        if 1900 <= year <= 2099:
            return match

    print("Year not found:", date_string)
    return "NA"


def clean_value(value):
    if pd.isna(value):
        return "NA"

    return str(value)\
        .replace('/', '-')\
        .replace(' ', '-')\
        .replace(',', '-')\
        .replace(')', '-')\
        .replace('(', '-')\
        .strip()


def create_string(row, columns):
    return "/".join(str(row[col]) for col in columns)


def normalize_id(x):
    return re.split(r"[/_]", str(x))[0]


def meta2fastaname(fasta_path, meta_path, hostmap_path, countrymap_path, meta_columns):
    meta = pd.read_csv(meta_path)

    meta["ID"] = meta["ID"].apply(normalize_id)

    meta = process_column(meta, hostmap_path, column="Host")
    meta = process_column(meta, countrymap_path, column="Geo location")

    meta['Year'] = meta['Collection date'].apply(extract_year_from_date)

    columns_to_include = meta_columns.split(",")

    for col in columns_to_include:
        meta[col] = meta[col].apply(clean_value)

    meta['info_string'] = meta.apply(
        lambda row: create_string(row, columns_to_include),
        axis=1
    )

    new_seq_names = dict(zip(meta["ID"], meta["info_string"]))

    seqs = list(SeqIO.parse(fasta_path, 'fasta'))
    print(f'Total number of records in FASTA: {len(seqs)}')

    new_seqs = []
    missing = 0

    for seq in seqs:
        ac = normalize_id(seq.name)

        if ac not in new_seq_names:
            print(f"WARNING: {ac} not found in metadata")
            missing += 1
            continue

        new_name = new_seq_names[ac]

        seq.id = new_name
        seq.name = new_name
        seq.description = ''

        new_seqs.append(seq)

    print(f'Total processed: {len(new_seqs)}')
    print(f'Missing in metadata: {missing}')

    out_path = os.path.splitext(fasta_path)[0] + '_ren.fasta'
    with open(out_path, 'w', encoding='utf-8') as out:
        SeqIO.write(new_seqs, out, 'fasta')

    print(f"Saved to: {out_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-in_fasta", required=True)
    parser.add_argument("-in_meta", required=True)
    parser.add_argument("-country_map", required=True)
    parser.add_argument("-host_map", required=True)
    parser.add_argument("-columns", required=True)

    args = parser.parse_args()

    meta2fastaname(
        args.in_fasta,
        args.in_meta,
        args.host_map,
        args.country_map,
        args.columns
    )