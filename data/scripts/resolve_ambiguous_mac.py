import argparse
import copy
import os
import pandas as pd
import re
import subprocess
import sys
from Bio import SeqIO

ambig_nt = [
    'n','r','y','k','m','s','w','b','d','h','v','total'
]

def resolve_ambiguos(input_file, output_dir, window, path_to_blast,
                     evalue, word_size):

    if sys.platform.startswith("win"):
        output_dir = output_dir.rstrip("\\") + "\\"
    else:
        output_dir = output_dir.rstrip("/") + "/"

    fasta_al = list(SeqIO.parse(input_file, "fasta"))
    window = 100

    print('------Finding sequences with ambiguous characters------')

    fasta_al_less_amb = []
    list_slices = []

    for rec in fasta_al.copy():
        amb_total = len(re.findall(r"[nrykmswbdhv]", str(rec.seq)))

        if amb_total == 0:
            fasta_al_less_amb.append(rec)
            continue

        rec_seq_len = len(re.sub("-", "", str(rec.seq)))
        print(f'{rec.id}: {amb_total} ({round(amb_total/rec_seq_len,2)})')

        if amb_total / rec_seq_len > 0.01:
            print(rec.id, 'exceeded threshold')
            continue

        if re.search(r"nnnnn", str(rec.seq)):
            print(rec.id, 'has >=5 n in a row')
            continue

        fasta_al_less_amb.append(rec)

        starts = [m.start() for m in re.finditer(r"[nrykmswbdhv]", str(rec.seq))]
        i = 0

        while i < len(starts):
            if starts[i] < window // 2:
                st, e = 0, window
            elif len(rec.seq) - starts[i] < window // 2:
                st = len(rec.seq) - window
                e = len(rec.seq)
            else:
                st = starts[i] - window // 2
                e = starts[i] + window // 2

            cur_slice_rec = copy.deepcopy(rec[st:e])
            cur_slice_rec.description = ''

            cur_slice_rec.id = f"{rec.id}_{st+1}:{starts[i]+1}:{e}"
            list_slices.append(cur_slice_rec)

            i += 1

    file_name_slices = os.path.splitext(input_file)[0] + "_slices.fasta"
    SeqIO.write(list_slices, file_name_slices, "fasta")

    file_name_less_amb = os.path.splitext(input_file)[0] + "_less_amb.fasta"
    SeqIO.write(fasta_al_less_amb, file_name_less_amb, "fasta")

    print('------Creating BLAST database------')

    if sys.platform.startswith("win"):
        makeblast_command = (
            f"{path_to_blast}makeblastdb.exe "
            f"-in {file_name_less_amb} "
            f"-dbtype nucl "
            f"-out {output_dir}local_db"
        )

        blastn_command = (
            f"{path_to_blast}blastn.exe "
            f"-db {output_dir}local_db "
            f"-query {file_name_slices} "
            f"-outfmt 6 "
            f"-out {output_dir}blast.out "
            f"-strand plus "
            f"-evalue {evalue} "
            f"-word_size {word_size} "
            f"-max_target_seqs 30"
        )
    else:
        makeblast_command = (
            f"{path_to_blast}makeblastdb "
            f"-in {file_name_less_amb} "
            f"-dbtype nucl "
            f"-out {output_dir}local_db"
        )

        blastn_command = (
            f"{path_to_blast}blastn "
            f"-db {output_dir}local_db "
            f"-query {file_name_slices} "
            f"-outfmt 6 "
            f"-out {output_dir}blast.out "
            f"-strand plus "
            f"-evalue {evalue} "
            f"-word_size {word_size} "
            f"-max_target_seqs 30"
        )

    subprocess.call(makeblast_command, shell=True)

    with open(output_dir + "blast.err", "w") as err:
        subprocess.call(blastn_command, shell=True, stderr=err)

    print("------Done------")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("-input", "--input_file", required=True)
    parser.add_argument("-pout", "--path_out")
    parser.add_argument("-w", "--window", default=100)
    parser.add_argument("-evalue", "--evalue", type=float, default=1e-20)
    parser.add_argument("-word_size", "--word_size", type=int, default=7)
    parser.add_argument("-pb", "--path_blast", required=True)

    args = parser.parse_args()

    args.input_file = os.path.realpath(args.input_file)

    if not args.path_out:
        args.path_out = os.path.dirname(args.input_file)

    resolve_ambiguos(
        args.input_file,
        args.path_out,
        args.window,
        args.path_blast,
        args.evalue,
        args.word_size
    )
