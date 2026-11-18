#!/usr/bin/env python3
"""Download the Open Tree of Life induced subtree for the 32 species below."""

import csv
import json
import re
import sys
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

BASE_URL = "https://api.opentreeoflife.org/v3"
OUTPUT_DIR = Path("OpenTree_results")
TREE_FILE = OUTPUT_DIR / "OpenTree_subtree.nwk"
MAPPING_FILE = OUTPUT_DIR / "OpenTree_taxa_mapping.csv"
NOT_MATCHED_FILE = OUTPUT_DIR / "OpenTree_not_matched.txt"
NOT_IN_TREE_FILE = OUTPUT_DIR / "OpenTree_not_in_tree.txt"

SPECIES = [
    "Alisterus scapularis",
    "Anas platyrhynchos",
    "Anser anser",
    "Anser cygnoides",
    "Arenaria interpres",
    "Cairina moschata",
    "Canis lupus",
    "Columba livia",
    "Corvus splendens",
    "Cygnus olor",
    "Duttaphrynus melanostictus",
    "Eolophus roseicapilla",
    "Gallus gallus",
    "Gekko similignum",
    "Gerygone igata",
    "Hemidactylus bowringii",
    "Heteronotia planiceps",
    "Leucogeranus leucogeranus",
    "Mareca falcata",
    "Meleagris gallopavo",
    "Montifringilla taczanowskii",
    "Neogale vison",
    "Pelophylax nigromaculatus",
    "Periparus ater",
    "Petroica macrocephala",
    "Phrynocephalus erythrurus",
    "Prunella montanella",
    "Recurvirostra novaehollandiae",
    "Macaca mulatta",
    "Rhipidura fuliginosa",
    "Spilopelia chinensis",
    "Sus scrofa"
]

# Match quoted labels and comments before delimiters, so commas inside them
# are never counted as separators.
NEWICK_TOKEN = re.compile(r"'(?:[^']|'')*'|\[[^\]]*\]|[(),:;]|[^\s'\[\](),:;]+")


def post(endpoint, payload, timeout):
    request = Request(
        f"{BASE_URL}/{endpoint}",
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"},
    )
    try:
        with urlopen(request, timeout=timeout) as response:
            return response.status, json.load(response)
    except HTTPError as error:
        with error:
            body = error.read().decode("utf-8", errors="replace")
        try:
            return error.code, json.loads(body)
        except ValueError:
            raise RuntimeError(f"{endpoint}: HTTP {error.code}") from error


def match_names():
    status, data = post(
        "tnrs/match_names",
        {"names": SPECIES, "context_name": "Animals", "do_approximate_matching": False},
        timeout=120,
    )
    if status != 200:
        raise RuntimeError(f"TNRS HTTP {status}: {data.get('message', data)}")

    results = {result["name"].casefold(): result["matches"] for result in data["results"]}
    mapping = []
    for name in SPECIES:
        matches = [
            match for match in results.get(name.casefold(), [])
            if not match.get("is_approximate_match", False)
        ]
        best = {}
        if matches:
            best_score = max(match["score"] for match in matches)
            best = {
                match["taxon"]["ott_id"]: match["taxon"]
                for match in matches if match["score"] == best_score
            }

        row = {"input_name": name, "matched_name": "", "ott_id": "", "status": "not matched"}
        if len(best) == 1:
            ott_id, taxon = next(iter(best.items()))
            row.update(matched_name=taxon["name"], ott_id=ott_id, status="matched")
            print(f"{name} -> {taxon['name']} (OTT {ott_id})")
        elif best:
            row["status"] = "ambiguous"
        mapping.append(row)
    return mapping


def get_subtree(ott_ids, excluded):
    remaining = set(ott_ids)
    while len(remaining) >= 2:
        status, data = post(
            "tree_of_life/induced_subtree",
            {"ott_ids": sorted(remaining), "label_format": "name_and_id"},
            timeout=180,
        )
        missing = set()
        # The v3 API uses OTT keys such as "ott123" in unknown and broken.
        for field in ("ott_ids_not_in_tree", "ott_ids_not_in_graph", "unknown", "broken"):
            for value in data.get(field) or ():
                match = re.fullmatch(r"(?:ott)?(\d+)", str(value))
                if match:
                    missing.add(int(match.group(1)))

        missing &= remaining
        if missing and status in (200, 400, 404):
            excluded.update(missing)
            remaining.difference_update(missing)
            continue
        if status != 200:
            raise RuntimeError(f"Induced subtree HTTP {status}: {data.get('message', data)}")
        newick = data.get("newick", "").strip()
        if not newick.endswith(";"):
            raise RuntimeError("OpenTree did not return a complete Newick tree.")
        return newick, remaining
    raise RuntimeError("Fewer than two distinct OTT IDs available for the subtree.")


def leaf_labels(newick):
    expecting_leaf = True
    labels = []
    for token in NEWICK_TOKEN.findall(newick):
        if token.startswith("["):
            continue
        if token in ("(", ","):
            expecting_leaf = True
            continue
        if expecting_leaf:
            if token in (")", ":", ";"):
                raise RuntimeError("Newick contains an unlabeled leaf.")
            labels.append(token.strip("'"))
        expecting_leaf = False
    return labels


def write_names(path, names):
    path.write_text("".join(f"{name}\n" for name in names), encoding="utf-8")


def main():
    OUTPUT_DIR.mkdir(exist_ok=True)
    print("Matching taxonomic names...")
    mapping = match_names()
    with MAPPING_FILE.open("w", encoding="utf-8", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=mapping[0], delimiter="\t")
        writer.writeheader()
        writer.writerows(mapping)

    unmatched = [row["input_name"] for row in mapping if row["ott_id"] == ""]
    write_names(NOT_MATCHED_FILE, unmatched)
    matched = {row["ott_id"] for row in mapping if row["ott_id"] != ""}
    excluded = set()
    print("Building induced subtree...")
    try:
        newick, expected_ids = get_subtree(matched, excluded)
        labels = leaf_labels(newick)
        leaf_ids = set()
        for label in labels:
            match = re.search(r"_ott(\d+)$", label)
            if match:
                leaf_ids.add(int(match.group(1)))

        excluded.update(expected_ids - leaf_ids)
        if len(labels) != len(expected_ids) or leaf_ids != expected_ids:
            raise RuntimeError(
                f"Newick validation failed: {len(labels)} leaves, "
                f"{len(expected_ids)} expected. "
                f"Missing OTT IDs: {sorted(expected_ids - leaf_ids)}. "
                f"Unexpected OTT IDs: {sorted(leaf_ids - expected_ids)}."
            )

        TREE_FILE.write_text(newick + "\n", encoding="utf-8")
        print(f"Saved {TREE_FILE} with {len(labels)}/{len(SPECIES)} requested leaves.")
        if len(labels) != len(SPECIES):
            print("Warning: the saved subtree is incomplete. See the reports.", file=sys.stderr)
    finally:
        write_names(
            NOT_IN_TREE_FILE,
            [row["input_name"] for row in mapping if row["ott_id"] in excluded],
        )


if __name__ == "__main__":
    try:
        main()
    except (RuntimeError, URLError, OSError, ValueError) as error:
        sys.exit(f"Error: {error}")
