
from pathlib import Path
import pandas as pd

base_dir = Path(
    "/Users/abagavetdinova/Desktop/lab/Astroviridae_database/"
    "data/snakemake/ORF_shift_analysis/prodigal"
)

gff_file = base_dir / "genes.gff"
csv_file = (
    base_dir.parent.parent
    / "constants"
    / "Astroviridae_Aves_Amphibia_14072026.csv"
)
output_file = base_dir / "first_long_orfs.csv"

def normalize_accession(value):
    accession = str(value).strip().split("/", 1)[0]
    if accession.startswith("NC_"):
        accession = accession.replace("NC_", "NC-", 1)
    return accession

metadata = pd.read_csv(csv_file)
metadata["ID"] = metadata["ID"].apply(normalize_accession)
metadata["1A-strand"] = pd.to_numeric(
    metadata["1A-strand"], errors="coerce"
)

strand_by_accession = (
    metadata.drop_duplicates("ID")
    .set_index("ID")["1A-strand"]
    .to_dict()
)

columns = [
    "accession", "start", "end", "strand",
    "1A-strand", "length_nt", "stop_start",
    "stop_end", "attributes"
]

rows = []

with open(gff_file) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue

        fields = line.rstrip().split("\t")
        if len(fields) != 9 or fields[2] != "CDS":
            continue

        accession = normalize_accession(fields[0])
        start = int(fields[3])
        end = int(fields[4])
        strand = fields[6]
        attributes = fields[8]
        length = end - start + 1

        if length <= 2000:
            continue

        strand_1a = strand_by_accession.get(accession)

        if strand == "+":
            stop_start = end - 2
            stop_end = end
        else:
            stop_start = start
            stop_end = start + 2

        rows.append({
            "accession": accession,
            "start": start,
            "end": end,
            "strand": strand,
            "1A-strand": strand_1a,
            "length_nt": length,
            "stop_start": stop_start,
            "stop_end": stop_end,
            "attributes": attributes
        })

df = pd.DataFrame(rows, columns=columns)

if not df.empty:
    # При 1A-strand = -1.0 выбираем ORF с максимальной координатой.
    # В остальных случаях — ORF с минимальной координатой.
    df["_sort_coord"] = df["start"]
    reverse_mask = df["1A-strand"].eq(-1.0)
    df.loc[reverse_mask, "_sort_coord"] = (
        -df.loc[reverse_mask, "end"]
    )

    df = df.sort_values(
        ["accession", "_sort_coord"],
        kind="stable"
    )

    df = df.drop_duplicates(
        subset="accession",
        keep="first"
    )

    df = df.drop(columns="_sort_coord")

df.to_csv(output_file, index=False)

print(f"Всего геномов с ORF > 2000 нт: {len(df)}")
print(
    "Геномов с 1A-strand = -1.0: "
    f"{(df['1A-strand'] == -1.0).sum()}"
)
print(
    "Геномов без определённого направления: "
    f"{df['1A-strand'].isna().sum()}"
)
print(f"Результат сохранён: {output_file}")
