import sys

import Bio.SeqIO

import pandas as pd


sys.stderr = sys.stdout = open(snakemake.log[0], "w")

recent_vaccine_strains = snakemake.params.recent_vaccine_strains

df = pd.read_csv(snakemake.input.csv)

assert set(recent_vaccine_strains).issubset(df["strain"])

df = df[
    (df["strain_type"] == "circulating_2025")
    | df["strain"].isin(recent_vaccine_strains)
]

for f, seq_type in [
    (snakemake.output.nt_fasta, "nt_sequence_HA_ectodomain"),
    (snakemake.output.prot_fasta, "protein_sequence_HA_ectodomain"),
]:
    assert df[seq_type].notnull().all()
    with open(f, "w") as f_obj:
        for head, seq in df[["strain", seq_type]].drop_duplicates().itertuples(index=False):
            if head in recent_vaccine_strains:
                head = head + " " + recent_vaccine_strains[head].replace(" ", "_")
            f_obj.write(f">{head}\n{seq}\n")
