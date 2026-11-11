
import re
import time
from pathlib import Path

import pandas as pd
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

# =========================
# SETTINGS
# =========================

Entrez.email = "adelina.kzn15@example.com"
# Если есть NCBI API key, укажи его:
# Entrez.api_key = "YOUR_NCBI_API_KEY"

OUTPUT_DIR = Path("PTGS1_results")
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

MAX_RESULTS = 100
REQUEST_DELAY = 0.4
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



def make_query(species):
    protein_query = (
        '("prostaglandin-endoperoxide synthase 1"[Title] '
        'OR "prostaglandin G/H synthase 1"[Title] '
        'OR "cyclooxygenase-1"[Title] '
        'OR "cyclooxygenase 1"[Title] '
        'OR PTGS1[Gene Name])'
    )

    return f'({protein_query}) AND "{species}"[Organism]'

def get_annotations(record):
    genes = []
    products = []
    notes = []

    for feature in record.features:
        genes.extend(feature.qualifiers.get("gene", []))
        products.extend(feature.qualifiers.get("product", []))
        notes.extend(feature.qualifiers.get("note", []))

    return genes, products, notes


def evaluate_record(record, requested_species):
    description = record.description.lower()

    genes, products, notes = get_annotations(record)

    genes = [x.lower() for x in genes]
    products = [x.lower() for x in products]
    notes = [x.lower() for x in notes]

    annotation = " ".join(
        [description] + genes + products + notes
    )

    # Exclude PTGS2
    if re.search(r"\bptgs2\b", annotation):
        return None

    # Exact PTGS1 gene annotation
    gene_match = "ptgs1" in genes

    # Protein name annotation
    product_match = any(
        term in annotation
        for term in [
            "prostaglandin-endoperoxide synthase 1",
            "prostaglandin g/h synthase 1",
            "cyclooxygenase-1",
            "cyclooxygenase 1",
            "cox-1"
        ]
    )

    if not gene_match and not product_match:
        return None

    organism = record.annotations.get("organism", "")

    # Check organism
    if organism.lower() != requested_species.lower():
        return None

    length = len(record.seq)
    accession = record.id

    partial = bool(
        re.search(r"\bpartial\b", annotation)
        or re.search(r"\bfragment\b", annotation)
    )

    refseq = accession.startswith("NP_")
    predicted = accession.startswith("XP_")

    # Ranking
    score = 0

    if gene_match:
        score += 100

    if product_match:
        score += 50

    if refseq:
        score += 30

    if not predicted:
        score += 10

    if not partial:
        score += 20

    # Prefer full-length proteins around 600 aa
    if 550 <= length <= 750:
        score += 30
    elif length >= 500:
        score += 10

    return {
        "species": requested_species,
        "accession": accession,
        "description": record.description,
        "organism": organism,
        "length": length,
        "gene_match": gene_match,
        "product_match": product_match,
        "refseq": refseq,
        "predicted": predicted,
        "partial": partial,
        "score": score,
        "record": record
    }


# =========================
# NCBI SEARCH
# =========================

def search_protein(species):
    query = make_query(species)

    try:
        with Entrez.esearch(
            db="protein",
            term=query,
            retmax=MAX_RESULTS,
            sort="relevance"
        ) as handle:
            search_result = Entrez.read(handle)

        ids = search_result["IdList"]

        if not ids:
            return []

        time.sleep(REQUEST_DELAY)

        with Entrez.efetch(
            db="protein",
            id=",".join(ids),
            rettype="gb",
            retmode="text"
        ) as handle:
            records = list(SeqIO.parse(handle, "genbank"))

        candidates = []

        for record in records:
            result = evaluate_record(record, species)

            if result is not None:
                candidates.append(result)

        candidates.sort(
            key=lambda x: (
                x["score"],
                x["length"]
            ),
            reverse=True
        )

        return candidates

    except Exception as error:
        print(f"ERROR: {species}: {error}")
        return []

def main():
    fasta_records = []
    metadata = []
    all_candidates = []
    missing = []

    print("=" * 60)
    print("PTGS1 SEARCH IN NCBI PROTEIN")
    print("=" * 60)

    for i, species in enumerate(species_list, 1):
        print(f"\n[{i}/{len(species_list)}] {species}")

        candidates = search_protein(species)

        if not candidates:
            print("NOT FOUND")
            missing.append(species)
            time.sleep(REQUEST_DELAY)
            continue

        # Save all candidates for review
        for candidate in candidates:
            all_candidates.append({
                "Host_species": species,
                "Accession": candidate["accession"],
                "Description": candidate["description"],
                "NCBI_organism": candidate["organism"],
                "Length_aa": candidate["length"],
                "Gene_match": candidate["gene_match"],
                "Product_match": candidate["product_match"],
                "RefSeq": candidate["refseq"],
                "Predicted": candidate["predicted"],
                "Partial": candidate["partial"],
                "Score": candidate["score"]
            })

        # Select best candidate
        best = candidates[0]
        record = best["record"]

        species_id = species.replace(" ", "_")
        fasta_id = f"{species_id}|{best['accession']}"

        fasta_record = SeqRecord(
            record.seq,
            id=fasta_id,
            description="PTGS1"
        )

        fasta_records.append(fasta_record)

        metadata.append({
            "Host_species": species,
            "Accession": best["accession"],
            "Description": best["description"],
            "NCBI_organism": best["organism"],
            "Length_aa": best["length"],
            "Gene_match": best["gene_match"],
            "Product_match": best["product_match"],
            "RefSeq": best["refseq"],
            "Predicted": best["predicted"],
            "Partial": best["partial"],
            "Score": best["score"],
            "Number_of_candidates": len(candidates)
        })

        print(f"FOUND: {best['accession']}")
        print(f"Length: {best['length']} aa")
        print(f"Candidates: {len(candidates)}")
        print(f"RefSeq: {best['refseq']}")
        print(f"Partial: {best['partial']}")

        time.sleep(REQUEST_DELAY)


    fasta_file = OUTPUT_DIR / "PTGS1_all_hosts.fasta"

    SeqIO.write(
        fasta_records,
        fasta_file,
        "fasta"
    )



    metadata_file = OUTPUT_DIR / "PTGS1_metadata.csv"

    pd.DataFrame(metadata).to_csv(
        metadata_file,
        index=False
    )


    candidates_file = OUTPUT_DIR / "PTGS1_all_candidates.csv"

    pd.DataFrame(all_candidates).to_csv(
        candidates_file,
        index=False
    )


    missing_file = OUTPUT_DIR / "PTGS1_missing_species.txt"

    with open(missing_file, "w") as file:
        file.write("\n".join(missing))


    print("\n" + "=" * 60)
    print("SEARCH COMPLETED")
    print("=" * 60)

    print(f"Total species: {len(species_list)}")
    print(f"Found: {len(fasta_records)}")
    print(f"Not found: {len(missing)}")

    print("\nOUTPUT FILES:")

    print(f"FASTA: {fasta_file}")
    print(f"Metadata: {metadata_file}")
    print(f"All candidates: {candidates_file}")
    print(f"Missing: {missing_file}")

    if missing:
        print("\nSPECIES WITHOUT PTGS1:")

        for species in missing:
            print(f"- {species}")


if __name__ == "__main__":
    main()
