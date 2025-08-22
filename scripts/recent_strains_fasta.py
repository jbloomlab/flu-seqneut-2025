import sys

import Bio.SeqIO

import pandas as pd


sys.stderr = sys.stdout = open(snakemake.log[0], "w")

df = pd.read_csv(snakemake.input.csv).query("strain_type == 'circulating_2025'")

for f, seq_type in [
    (snakemake.output.nt_fasta, "nt_sequence_HA_ectodomain"),
    (snakemake.output.prot_fasta, "protein_sequence_HA_ectodomain"),
]:
    assert df[seq_type].notnull().all()
    with open(f, "w") as f_obj:
        for head, seq in df[["strain", seq_type]].drop_duplicates().itertuples(index=False):
            f_obj.write(f">{head}\n{seq}\n")
