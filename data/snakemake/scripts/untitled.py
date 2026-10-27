from pathlib import Path
import pandas as pd

base_dir = Path(
    "/Users/abagavetdinova/Desktop/lab/Astroviridae_database/"
    "data/snakemake/ORF_shift_analysis/prodigal"
)

gff_file = base_dir / "genes.gff"
output_file = base_dir / "first_long_orfs.csv"

columns = [
    "accession", "start", "end", "strand",
    "length_nt", "stop_start", "stop_end", "attributes"
]

rows = []

with open(gff_file) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue

        fields = line.rstrip().split("\t")
        if len(fields) != 9 or fields[2] != "CDS":
            continue

        accession = fields[0]
        start = int(fields[3])
        end = int(fields[4])
        strand = fields[6]
        attributes = fields[8]
        length = end - start + 1

        if length <= 2000:
            continue

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
            "length_nt": length,
            "stop_start": stop_start,
            "stop_end": stop_end,
            "attributes": attributes
        })

df = pd.DataFrame(rows, columns=columns)

if not df.empty:
    df = df.sort_values(
        ["accession", "start"],
        kind="stable"
    )
    df = df.drop_duplicates(
        subset="accession",
        keep="first"
    )

df.to_csv(output_file, index=False)

print(f"Всего геномов с ORF > 2000 нт: {len(df)}")
print(f"Результат сохранён: {output_file}")