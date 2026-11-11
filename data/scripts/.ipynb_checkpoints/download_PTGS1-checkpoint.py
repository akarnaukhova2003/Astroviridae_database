
import csv
import time
from pathlib import Path

from Bio import Entrez, SeqIO


# =========================
# SETTINGS
# =========================

Entrez.email = "adelina.kzn@gmail.com"

OUTPUT_DIR = Path("PTGS1_results")
OUTPUT_DIR.mkdir(exist_ok=True)

FASTA_FILE = OUTPUT_DIR / "PTGS1_all_sequences.fasta"
CSV_FILE = OUTPUT_DIR / "PTGS1_all_candidates.csv"
MISSING_FILE = OUTPUT_DIR / "PTGS1_missing_species.txt"

BATCH_SIZE = 200
SEARCH_PAGE_SIZE = 5000
REQUEST_DELAY = 0.4


# =========================
# SPECIES
# =========================

species_list = [
    "Gallus gallus",
    "Anser anser",
    "Meleagris gallopavo",
    "Anas platyrhynchos",
    "Columba livia",
    "Neogale vison",
    "Cygnus olor",
    "Gerygone igata",
    "Spilopelia chinensis",
    "Arenaria interpres",
    "Leucogeranus leucogeranus",
    "Canis lupus",
    "Rhipidura fuliginosa",
    "Corvus splendens",
    "Recurvirostra novaehollandiae",
    "Sus scrofa",
    "Cairina moschata",
    "Eolophus roseicapilla",
    "Alisterus scapularis",
    "Mareca falcata",
    "Heteronotia planiceps",
    "Pelophylax nigromaculatus",
    "Montifringilla taczanowskii",
    "Phrynocephalus erythrurus",
    "Prunella montanella",
    "Petroica macrocephala",
    "Macaca mulatta",
    "Hemidactylus bowringii",
    "Gekko similignum",
    "Duttaphrynus melanostictus"
]


# =========================
# SEARCH TERMS
# =========================

SEARCH_TERMS = (
    '("prostaglandin-endoperoxide synthase 1"[Title] '
    'OR "prostaglandin G/H synthase 1"[Title] '
    'OR "cyclooxygenase-1"[Title] '
    'OR "cyclooxygenase 1"[Title] '
    'OR "COX-1"[Title] '
    'OR PTGS1[Gene Name])'
)


PTGS1_TERMS = [
    "prostaglandin-endoperoxide synthase 1",
    "prostaglandin g/h synthase 1",
    "cyclooxygenase-1",
    "cyclooxygenase 1",
    "cox-1"
]


# =========================
# HELPERS
# =========================

def get_record_text(record):
    """Collect text from protein record annotations and features."""
    parts = [record.description]

    for feature in record.features:
        for key, values in feature.qualifiers.items():
            parts.extend(values)

    return " ".join(str(part) for part in parts).lower()


def get_gene_names(record):
    """Extract gene qualifiers."""
    genes = []

    for feature in record.features:
        genes.extend(feature.qualifiers.get("gene", []))

    return [gene.lower().strip() for gene in genes]


def is_ptgs2(record, text, genes):
    """Exclude PTGS2 / COX-2 records."""
    if any(gene == "ptgs2" for gene in genes):
        return True

    if "ptgs2" in text:
        return True

    if "cyclooxygenase-2" in text or "cyclooxygenase 2" in text:
        return True

    if "cox-2" in text:
        return True

    return False


def is_ptgs1(record):
    """Check whether the record corresponds to PTGS1."""
    text = get_record_text(record)
    genes = get_gene_names(record)

    if is_ptgs2(record, text, genes):
        return False, False, False

    gene_match = any(gene == "ptgs1" for gene in genes)

    product_match = any(term in text for term in PTGS1_TERMS)

    return gene_match or product_match, gene_match, product_match


def get_organism(record):
    return record.annotations.get("organism", "")


def is_partial(record):
    """Mark records that appear to be partial."""
    text = get_record_text(record)

    partial_terms = [
        "partial",
        "fragment",
        "incomplete"
    ]

    return any(term in text for term in partial_terms)


def search_protein_ids(species):
    """Retrieve all Protein IDs returned by NCBI for a species."""
    query = f'({SEARCH_TERMS}) AND "{species}"[Organism]'

    with Entrez.esearch(
        db="protein",
        term=query,
        retmax=0
    ) as handle:
        search_result = Entrez.read(handle)

    total = int(search_result["Count"])
    ids = []

    for start in range(0, total, SEARCH_PAGE_SIZE):
        with Entrez.esearch(
            db="protein",
            term=query,
            retstart=start,
            retmax=SEARCH_PAGE_SIZE
        ) as handle:
            result = Entrez.read(handle)

        ids.extend(result["IdList"])
        time.sleep(REQUEST_DELAY)

    return ids


def fetch_records(ids):
    """Download GenBank records in batches."""
    for start in range(0, len(ids), BATCH_SIZE):
        batch = ids[start:start + BATCH_SIZE]

        with Entrez.efetch(
            db="protein",
            id=",".join(batch),
            rettype="gb",
            retmode="text"
        ) as handle:
            records = list(SeqIO.parse(handle, "genbank"))

        yield from records
        time.sleep(REQUEST_DELAY)


# =========================
# MAIN
# =========================

all_records = []
metadata = []
missing_species = []

seen_accessions = set()

for index, species in enumerate(species_list, start=1):
    print(f"\n[{index}/{len(species_list)}] Searching: {species}")

    try:
        protein_ids = search_protein_ids(species)

        print(f"NCBI records found: {len(protein_ids)}")

        species_records = []

        for record in fetch_records(protein_ids):
            accession = record.id

            if accession in seen_accessions:
                continue

            organism = get_organism(record)

            if organism.lower() != species.lower():
                continue

            matched, gene_match, product_match = is_ptgs1(record)

            if not matched:
                continue

            seen_accessions.add(accession)
            species_records.append(record)

            description = record.description
            length = len(record.seq)

            metadata.append({
                "Host_species": species,
                "Accession": accession,
                "Description": description,
                "NCBI_organism": organism,
                "Length_aa": length,
                "Gene_match": gene_match,
                "Product_match": product_match,
                "RefSeq": accession.startswith(("NP_", "XP_")),
                "Predicted": accession.startswith("XP_"),
                "Partial": is_partial(record)
            })

        all_records.extend(species_records)

        print(f"PTGS1 sequences found: {len(species_records)}")

        if not species_records:
            missing_species.append(species)

    except Exception as error:
        print(f"ERROR for {species}: {error}")
        missing_species.append(species)

    time.sleep(REQUEST_DELAY)


# =========================
# SAVE FASTA
# =========================

if all_records:
    SeqIO.write(all_records, FASTA_FILE, "fasta")

print(f"\nFASTA saved: {FASTA_FILE}")


# =========================
# SAVE METADATA CSV
# =========================

fieldnames = [
    "Host_species",
    "Accession",
    "Description",
    "NCBI_organism",
    "Length_aa",
    "Gene_match",
    "Product_match",
    "RefSeq",
    "Predicted",
    "Partial"
]

with open(CSV_FILE, "w", newline="", encoding="utf-8-sig") as file:
    writer = csv.DictWriter(file, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(metadata)

print(f"Metadata saved: {CSV_FILE}")


# =========================
# SAVE MISSING SPECIES TXT
# =========================

with open(MISSING_FILE, "w", encoding="utf-8") as file:
    file.write("\n".join(missing_species))

print(f"Missing species saved: {MISSING_FILE}")


# =========================
# SUMMARY
# =========================

print("\n========== SUMMARY ==========")

print(f"Species searched: {len(species_list)}")
print(f"Species with PTGS1: {len(species_list) - len(missing_species)}")
print(f"Species without PTGS1: {len(missing_species)}")
print(f"Total PTGS1 sequences: {len(all_records)}")

print("\nMissing species:")
for species in missing_species:
    print(f"- {species}")

print("\nOutput directory:", OUTPUT_DIR.resolve())
