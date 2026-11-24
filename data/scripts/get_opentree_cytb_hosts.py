#!/usr/bin/env python3
"""Download an Open Tree of Life induced subtree for species listed in a text file."""

import argparse
import csv
import json
import re
import sys
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen

BASE_URL = "https://api.opentreeoflife.org/v3"
# Query names used by Open Tree Taxonomy. Input names are kept in the report.
QUERY_NAMES = {
    "leucogeranus leucogeranus": "Grus leucogeranus",
    "neogale vison": "Neovison vison",
    "rhesus monkey": "Macaca mulatta",
}


def read_species(path):
    species = {}
    for line in path.read_text(encoding="utf-8-sig").splitlines():
        name = line.strip()
        if name and not name.startswith("#"):
            species.setdefault(name.casefold(), name)
    if len(species) < 2:
        raise ValueError("The input file must contain at least two distinct species.")
    return list(species.values())


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


def match_names(species):
    query_names = {name: QUERY_NAMES.get(name.casefold(), name) for name in species}
    status, data = post(
        "tnrs/match_names",
        {"names": list(dict.fromkeys(query_names.values())), "context_name": "Animals", "do_approximate_matching": False},
        timeout=120,
    )
    if status != 200:
        raise RuntimeError(f"TNRS HTTP {status}: {data.get('message', data)}")

    results = {result["name"].casefold(): result["matches"] for result in data["results"]}
    mapping = []
    for name in species:
        matches = [
            match for match in results.get(query_names[name].casefold(), [])
            if not match.get("is_approximate_match", False)
        ]
        row = {
            "input_name": name, "matched_name": "", "ott_id": "",
            "status": "not matched", "score": "", "equal_score_candidates": "",
        }
        if matches:
            best_score = max(match["score"] for match in matches)
            candidates = {}
            for match in matches:
                if match["score"] == best_score:
                    ott_id = match["taxon"]["ott_id"]
                    if ott_id not in candidates or not match["is_synonym"]:
                        candidates[ott_id] = match

            # A smaller OTT ID provides a reproducible choice after a tie.
            ordered = sorted(candidates, key=lambda ott_id: (candidates[ott_id]["is_synonym"], ott_id))
            ott_id = ordered[0]
            selected = candidates[ott_id]
            taxon = selected["taxon"]
            status = "matched"
            if len(ordered) > 1:
                status = "matched: direct preferred"
                if selected["is_synonym"] == candidates[ordered[1]]["is_synonym"]:
                    status = "ambiguous: provisional"
                    print(f"WARNING: {name}: provisional OTT {ott_id}; review the mapping table.",
                          file=sys.stderr)
                row["equal_score_candidates"] = "; ".join(
                    f"{candidates[candidate]['taxon']['name']} (OTT {candidate}, "
                    f"{'synonym' if candidates[candidate]['is_synonym'] else 'direct'})"
                    for candidate in ordered
                )
            row.update(matched_name=taxon["name"], ott_id=ott_id, status=status, score=best_score)
            print(f"{name} -> {taxon['name']} (OTT {ott_id})")
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


def figtree_newick(newick, mapping):
    """Render the OpenTree topology with input species names and display lengths."""
    from io import StringIO

    from Bio import Phylo

    input_names = {}
    for row in mapping:
        if row["ott_id"] != "":
            input_names.setdefault(row["ott_id"], row["input_name"])

    tree = Phylo.read(StringIO(newick), "newick")

    def render(clade):
        if clade.is_terminal():
            ott_id = int(re.search(r"_ott(\d+)$", clade.name)[1])
            name = input_names[ott_id].replace("'", "''")
            return f"'{name}'"

        # Unary nodes do not define a new grouping of leaves.
        if len(clade.clades) == 1:
            return render(clade.clades[0])

        children = ",".join(f"{render(child)}:1" for child in clade.clades)
        return f"({children})"

    # OpenTree supplies topology only. The unit lengths are for display.
    return render(tree.root) + ";\n"


def write_names(path, names):
    path.write_text("".join(f"{name}\n" for name in names), encoding="utf-8")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("species_file", type=Path, help="UTF-8 text file, one species per line")
    parser.add_argument("--output-dir", type=Path, default=Path("OpenTree_CYTB_hosts"),
                        help="Directory for the tree and reports (default: OpenTree_CYTB_hosts)")
    args = parser.parse_args()

    species = read_species(args.species_file)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    tree_file = args.output_dir / "OpenTree_subtree.nwk"
    mapping_file = args.output_dir / "OpenTree_taxa_mapping.csv"
    not_matched_file = args.output_dir / "OpenTree_not_matched.txt"
    not_in_tree_file = args.output_dir / "OpenTree_not_in_tree.txt"
    print("Matching taxonomic names...")
    mapping = match_names(species)
    with mapping_file.open("w", encoding="utf-8", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=mapping[0], delimiter="\t")
        writer.writeheader()
        writer.writerows(mapping)

    unmatched = [row["input_name"] for row in mapping if row["ott_id"] == ""]
    write_names(not_matched_file, unmatched)
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

        tree_file.write_text(figtree_newick(newick, mapping), encoding="utf-8")
        print(f"Saved {tree_file} with {len(labels)}/{len(species)} requested leaves.")
        if len(labels) != len(species):
            print("Warning: the saved subtree is incomplete. See the reports.", file=sys.stderr)
    finally:
        write_names(
            not_in_tree_file,
            [row["input_name"] for row in mapping if row["ott_id"] in excluded],
        )


if __name__ == "__main__":
    try:
        main()
    except (RuntimeError, URLError, OSError, ValueError) as error:
        sys.exit(f"Error: {error}")
