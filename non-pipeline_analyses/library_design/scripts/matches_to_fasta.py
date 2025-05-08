"""Convert CSV of matches to FASTA, making sure all have matches."""


import pandas as pd

input_csv = "results/strains_for_library/match_prot_to_genbank_nt.csv"
print(f"Reading matches from {input_csv=}")
df = pd.read_csv(input_csv)

if (
    any(df["nt_sequence"].isnull())
    or any(df["prot_sequence"].map(len) * 3 != df["nt_sequence"].map(len))
):
    raise ValueError(f"Not a match for all sequences")

seqs = {"H3": [], "H1": []}
for tup in df.itertuples():
    s = f">{tup.strain} protein identical to {tup.accession_w_aa_muts_added}\n{tup.nt_sequence}\n"
    if tup.strain.split("_")[-1] == "H3N2":
        seqs["H3"].append(s)
    elif tup.strain.split("_")[-1] == "H1N1":
        seqs["H1"].append(s)
    else:
        raise ValueError(f"Cannot determine subtype from {tup.strain=}")

for desc, seqlist in [
    ("h3", seqs["H3"]), ("h1", seqs["H1"]), ("h3_and_h1", seqs["H3"] + seqs["H1"])
]:
    fastafile = f"results/strains_for_library/{desc}_nt_seqs_for_library.fasta"
    print(f"Writing in FASTA format to {fastafile=}")
    with open(fastafile, "w") as f:
        f.write("\n".join(seqlist))
