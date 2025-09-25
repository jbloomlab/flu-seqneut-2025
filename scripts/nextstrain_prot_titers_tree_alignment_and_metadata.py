import datetime
import sys

import numpy

import pandas as pd


sys.stdout = sys.stderr = open(snakemake.log[0], "w")

subtypes = snakemake.params.subtypes
circulating_strain_type = snakemake.params.circulating_strain_type
recent_vaccine_strains = snakemake.params.recent_vaccine_strains
frac_below_cols = snakemake.params.frac_below_cols
prefix_alignment = snakemake.params.prefix_alignment

viruses = (
    pd.read_csv(snakemake.input.viral_libraries_csv)
    [
        [
            "strain",
            "subtype",
            "strain_type",
            "protein_sequence_HA_ectodomain",
            "subclade",
            "collection_date",
        ]
    ]
    .drop_duplicates()
)

assert len(viruses) == viruses["strain"].nunique()

assert set(recent_vaccine_strains).issubset(viruses["strain"])

df = (
    viruses[
        (viruses["strain_type"] == circulating_strain_type) 
        | viruses["strain"].isin(recent_vaccine_strains)
    ]
    .assign(
        strain_type=lambda x: x["strain"].map(recent_vaccine_strains).where(
            x["strain"].isin(recent_vaccine_strains), x["strain_type"]
        )
    )
)

print(f"{len(df)=} of {len(viruses)} are {circulating_strain_type=} or {recent_vaccine_strains=}")
    
# ensure collection_date is in valid format
year = datetime.datetime.now().year
if all((df["collection_date"] > year - 100) & (df["collection_date"] < year + 1)):
    df = df.rename(columns={"collection_date": "date"})
else:
    raise ValueError(f"Not valid numerical dates in {df['collection_date']=}")

titers = (
    pd.read_csv(snakemake.input.summarized_titers_csv)
    .rename(columns={"virus": "strain"})
)
assert set(frac_below_cols).issubset(titers.columns), f"{frac_below_cols=}, {titers.columns=}"
assert set(titers["strain"]).issubset(df["strain"])
for col in ["median_titer"] + frac_below_cols:
    assert col not in df.columns, col
    df = df.merge(
        (
            titers
            .assign(serum_group=lambda x: f"{col}_" + x["serum_group"] + "_sera")
            .pivot_table(index="strain", values=col, columns="serum_group")
            .reset_index()
        ),
        on="strain",
        how="left",
        validate="one_to_one",
    )

titers = (
    pd.read_csv(snakemake.input.titers_csv)
    [["serum", "virus", "titer"]]
    .rename(columns={"virus": "strain"})
    .assign(log2_titer=lambda x: numpy.log(x["titer"]) / numpy.log(2))
)
sera_metadata = pd.read_csv(snakemake.input.sera_metadata_csv)[
    ["serum", "group", "collection_date_numerical", "age_years", "sex"]
]
assert len(sera_metadata) == sera_metadata["serum"].nunique()
assert set(titers["serum"]).issubset(sera_metadata["serum"])
assert len(titers) == len(titers[["serum", "strain"]].drop_duplicates())
titers = titers.merge(sera_metadata, on="serum", validate="many_to_one")
print(f"Read {len(titers)=} for {titers['serum'].nunique()=}")

for subtype in subtypes:
    print(f"\nProcessing {subtype=}")
    subtype_df = df[df["subtype"] == subtype].drop(columns="subtype")
    print(f"{len(subtype_df)=} of {len(df)=} are {subtype=}")

    subtype_titers = titers[titers["strain"].isin(subtype_df["strain"])]
    print(f"{len(subtype_titers)=} of {len(titers)=} are {subtype=} retained strains")

    # remove subtype suffix if present 
    strain_rename = {
        s: (s[: -len(subtype) - 1] if s.endswith(f"_{subtype}") else s)
        for s in subtype_df["strain"]
    }
    subtype_df["strain"] = subtype_df["strain"].map(strain_rename)
    assert len(subtype_df) == subtype_df["strain"].nunique()
    assert len(subtype_df) > 0, f"{len(subtype_df)=} for {subtype=}"

    alignment_file = snakemake.output[f"alignment_{subtype}"]
    metadata_file = snakemake.output[f"metadata_{subtype}"]
    titers_file = snakemake.output[f"titers_{subtype}"]

    print(f"Writing alignment to {alignment_file=}")
    with open(alignment_file, "w") as f:
        for tup in subtype_df.itertuples():
            seq = prefix_alignment[subtype] + tup.protein_sequence_HA_ectodomain
            f.write(f">{tup.strain}\n{seq}\n")

    print(f"Writing metadata to {metadata_file=}")
    (
        subtype_df
        .drop(columns=["protein_sequence_HA_ectodomain"])
        .to_csv(metadata_file, index=False, sep="\t", float_format="%.6g")
    )

    print(f"Writing titers to {titers_file=}")
    print(strain_rename)
    (
        subtype_titers
        .assign(strain=lambda x: x["strain"].map(strain_rename))
        .to_csv(titers_file, sep="\t", index=False, float_format="%.6g")
    )
