import argparse
from pathlib import Path

import pandas as pd
from Bio import SeqIO

parser = argparse.ArgumentParser(
    description="Выбрать по одной последовательности CYTB на вид и сохранить FASTA."
)
parser.add_argument(
    "results_dir",
    nargs="?",
    type=Path,
    default=Path(__file__).resolve().parent / "CYTB_results",
    help="Папка с исходными CSV и FASTA. По умолчанию CYTB_results рядом со скриптом.",
)
results_dir = parser.parse_args().results_dir
candidates_path = results_dir / "CYTB_all_candidates.csv"
source_fasta = results_dir / "CYTB_all_sequences.fasta"
output_fasta = results_dir / "CYTB_selected_sequences.fasta"

for path in (candidates_path, source_fasta):
    if not path.is_file():
        parser.exit(1, f"Файл не найден: {path}\n")

candidates = pd.read_csv(candidates_path, sep=None, engine="python", encoding="utf-8-sig")
species_column = "Species" if "Species" in candidates.columns else "NCBI_organism"
required = {species_column, "Accession", "Length_aa", "Partial"}
missing_columns = required - set(candidates.columns)
if missing_columns:
    parser.exit(1, "В CSV отсутствуют столбцы: " + ", ".join(sorted(missing_columns)) + "\n")

candidates[species_column] = candidates[species_column].astype("string").str.strip()
candidates["Accession"] = candidates["Accession"].astype("string").str.strip()
if candidates[[species_column, "Accession"]].fillna("").eq("").any().any():
    parser.exit(1, "В CSV есть пустые названия видов или идентификаторы записей.\n")

candidates["Partial"] = (
    candidates["Partial"].astype(str).str.strip().str.lower().isin(["true", "1", "yes"])
)
candidates["Length_difference"] = (pd.to_numeric(candidates["Length_aa"]) - 380).abs()

selected = candidates.sort_values(
    by=["Partial", "Length_difference"], ascending=[True, True]
).drop_duplicates(subset=species_column, keep="first")

wanted = set(selected["Accession"])
with source_fasta.open(encoding="utf-8") as handle:
    records = {
        record.id: record
        for record in SeqIO.parse(handle, "fasta")
        if record.id in wanted
    }

missing = wanted - records.keys()
if missing:
    parser.exit(1, "В исходном FASTA отсутствуют записи: " + ", ".join(sorted(missing)) + "\n")

output_records = []
for accession, species in zip(selected["Accession"], selected[species_column]):
    record = records[accession][:]
    record.id = "_".join(species.split())
    record.name = record.id
    record.description = ""
    output_records.append(record)

count = SeqIO.write(output_records, output_fasta, "fasta")
print(f"Всего кандидатов: {len(candidates)}")
print(f"Записано последовательностей: {count}")
print(f"FASTA сохранён: {output_fasta.resolve()}")
