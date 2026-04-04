import argparse
import os
import re
from Bio import SeqIO
import pandas as pd


def create_mapping_dict(mapping_file_path):
    """
    Create a dictionary mapping regex patterns to short codes.
    """
    mapping_df = pd.read_csv(mapping_file_path, sep='\t', header=None,
                             names=['regex_pattern', 'short_code'])

    pattern_to_code = dict(zip(mapping_df['regex_pattern'], mapping_df['short_code']))

    compiled_patterns = []
    for _, row in mapping_df.iterrows():
        try:
            pattern = re.compile(str(row['regex_pattern']), re.IGNORECASE)
            compiled_patterns.append((pattern, row['short_code']))
        except re.error as e:
            print(f"Warning: Invalid regex pattern '{row['regex_pattern']}': {e}")

    return pattern_to_code, compiled_patterns


def map_value_to_short_code(value, compiled_patterns, default="NA"):
    if pd.isna(value):
        return default

    for pattern, short_code in compiled_patterns:
        if pattern.search(value):
            return short_code

    print("{} was not found in mapping file".format(value))
    return value.replace(" ", "-")


def process_column(meta_df, mapping_file_path, column="Host"):
    pattern_to_code, compiled_patterns = create_mapping_dict(mapping_file_path)

    meta_processed = meta_df.copy()

    meta_processed[column + "_original"] = meta_processed[column]
    meta_processed[column] = meta_processed[column].apply(
        lambda x: map_value_to_short_code(x, compiled_patterns)
    )

    return meta_processed


# EXTRACTING YEAR FROM COLLECTION DATE
def extract_year_from_date(date_string):
    if pd.isna(date_string):
        return "NA"

    flexible_pattern = r'\b\d{4}\b'
    matches = re.findall(flexible_pattern, date_string)

    if matches:
        for match in matches:
            year_int = int(match)
            if 1900 <= year_int <= 2099:
                return match
        return matches[0]

    print("Year was not found")
    print(date_string)
    return date_string


def clean_value(value):
    """Clean a single value: handle NaN and replace unwanted characters"""
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
    parts = [str(row[col]) for col in columns]
    return "/".join(parts)


def meta2fastaname(fasta_path, meta_path, hostmap_path, countrymap_path, meta_columns):
    meta = pd.read_csv(meta_path)

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
    print('Total number of records: {}'.format(len(seqs)))

    new_seqs = []
    k = 0

    for seq in seqs:
        ac = seq.name.split('_')[0]
        seq.name = new_seq_names[ac]
        seq.id = new_seq_names[ac]
        seq.description = ''
        new_seqs.append(seq)
        k += 1

    print('Total number of processed records: {}'.format(k))

    with open(os.path.splitext(fasta_path)[0] + '_ren.fasta', 'w', encoding='utf-8') as output_handle:
        SeqIO.write(new_seqs, output_handle, 'fasta')


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-in_fasta", type=str,
                        help="Input alignment with sequences in fasta format.", required=True)
    parser.add_argument("-in_meta", type=str,
                        help="Input metafile", required=True)
    parser.add_argument("-country_map", type=str,
                        help="TSV file with regular expressions for countries and their ISO codes", required=True)
    parser.add_argument("-host_map", type=str,
                        help="TSV file with regular expressions for host names and their short names", required=True)
    parser.add_argument("-columns", type=str,
                        help="Columns to be included in new sequence names, separated by comma",
                        required=True)

    args = parser.parse_args()
    meta2fastaname(args.in_fasta, args.in_meta,
                   args.host_map, args.country_map, args.columns)