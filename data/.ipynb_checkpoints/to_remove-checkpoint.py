#!/usr/bin/env python3

def load_ids(ids_file):
    ids = set()
    with open(ids_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            clean_id = line.split("/")[0]
            ids.add(clean_id)
    return ids


def filter_fasta(input_fasta, ids_to_remove, output_fasta):
    with open(input_fasta) as infile, open(output_fasta, "w") as outfile:
        write_block = True
        header = None
        sequence_lines = []

        for line in infile:
            line = line.rstrip()

            if line.startswith(">"):
                # если есть предыдущая запись — записываем её
                if header is not None and write_block:
                    outfile.write(header + "\n")
                    outfile.write("\n".join(sequence_lines) + "\n")

                header = line
                fasta_id = line[1:].split("_")[0]  # до первого "_"
                write_block = fasta_id not in ids_to_remove
                sequence_lines = []
            else:
                sequence_lines.append(line)

        if header is not None and write_block:
            outfile.write(header + "\n")
            outfile.write("\n".join(sequence_lines) + "\n")


if __name__ == "__main__":
    ids_to_remove = load_ids("ids.txt")
    filter_fasta("input.fasta", ids_to_remove, "output.fasta")