"""Aggregate human sera metadata and titers."""

import datetime
import re
import sys

import pandas as pd


sys.stdout = sys.stderr = open(snakemake.log[0], "w")

# get variables from `snakemake
groups_metadata = snakemake.input.groups_metadata
group_titers = snakemake.input.group_titers
groups = snakemake.params.groups
human_sera_to_exclude = snakemake.params.human_sera_to_exclude
metadata_csv = snakemake.output.metadata_csv
titers_csv = snakemake.output.titers_csv


def age_years(age):
    if isinstance(age, (int, float)):
        return age
    elif re.fullmatch(r"\d+(?:\.\d*)", age):
        return float(age)
    elif re.fullmatch(r"\d+to\d+", age):
        return sum(map(float, age.split("to"))) / 2
    elif re.fullmatch(r"\d+\-\d+ y", age):
        return sum(map(float, age.split()[0].split("-"))) / 2
    elif age == "<6 m":
        return 0.25
    elif age == "6 m-2 y":
        return 1.25
    else:
        raise ValueError(f"Cannot parse {age=}")


def parse_month_year(s):
    if " to " in s:
        s = s.split(" to ")[-1]
    for fmt in ("%b-%Y", "%B-%Y"):
        try:
            return datetime.datetime.strptime(s, fmt)
        except ValueError:
            pass
    raise ValueError(f"Not a recognized month-year: {s}")


dfs = []
for group, group_csv in zip(groups, groups_metadata, strict=True):
    print(f"\nReading metadata: {group=}, {group_csv=}")

    df = pd.read_csv(group_csv).assign(
        age_years=lambda x: x["age"].map(age_years),
        collection_date_numerical=lambda x: x["collection_date"].map(parse_month_year),
    )

    print(f"Read metadata for {len(df)=} sera")
    assert (df["species"] == "human").all(), f"{group=}, {df['species'].unique()=}"
    assert df["sex"].isin({"M", "F"}).all(), f"{group=}, {df['sex'].unique()=}"
    assert (df["cohort"] == group).all(), f"{group=}, {df['cohort'].unique()=}"
    assert len(df) == df["bloom_lab_id"].nunique(), f"non-unique serum IDs for {group=}"
    assert "exposure_info" in df.columns

    dfs.append(
        df.rename(
            columns={
                "cohort": "group",
                "bloom_lab_id": "serum",
                "age": "age_string",
                "collection_date": "collection_date_string",
            }
        )[
            [
                "group",
                "serum",
                "collection_date_numerical",
                "age_years",
                "sex",
                "collection_date_string",
                "age_string",
                "exposure_info",
            ]
        ]
    )

metadata_df = pd.concat(dfs, ignore_index=True)
print(f"\nOverall read metadata for {len(metadata_df)=} sera")
assert len(metadata_df) == metadata_df["serum"].nunique()

metadata_df


titers_df = []
for group, group_csv in zip(groups, group_titers, strict=True):
    print(f"\nReading titers: {group=}, {group_csv=}")
    df = pd.read_csv(group_csv)
    print(f"Read {len(df)=} titers")
    titers_df.append(df)

titers_df = pd.concat(titers_df, ignore_index=True)

print(f"\nOverall read {len(titers_df)=} titers")

titers_df


# Get non-excluded sera with both metadata and titers, and write to CSVs

# first some sanity checks to make sure the serum and groups are matched
# as expected between metadata and titers
extra_sera_w_titers = set(titers_df["serum"]) - set(metadata_df["serum"])
if extra_sera_w_titers:
    raise ValueError(f"Some sera w titers lack metadata:\n{extra_sera_w_titers=}")
assert set(titers_df.columns).intersection(metadata_df.columns) == {"group", "serum"}
assert (
    len(titers_df[["serum", "group"]].drop_duplicates()) == titers_df["serum"].nunique()
)
sera_groups = metadata_df[["serum", "group"]].merge(
    titers_df[["serum", "group"]].drop_duplicates(),
    on="serum",
    how="outer",
    suffixes=["_metadata", "_titers"],
)
if len(sera_groups) != metadata_df["serum"].nunique():
    raise ValueError(
        "inconsistent serum-group assignments:\n"
        + str(
            sera_groups.query("group_metadata != group_titers")
            .query("group_metadata.notnull()")
            .query("group_titers.notnull()")
        )
    )


print(f"Excluding these sera:\n{human_sera_to_exclude=}")
if not set(human_sera_to_exclude).issubset(metadata_df["serum"]):
    raise ValueError(
        "Excluded sera lack metadata, maybe you specified names wrong?\n"
        f"{set(human_sera_to_exclude) - set(metadata_df['serum'])=}"
    )
else:
    filtered_titers_df = titers_df.query(
        "serum not in @human_sera_to_exclude"
    ).sort_values(["group", "serum", "virus"])
    assert len(filtered_titers_df) == len(
        filtered_titers_df.groupby(["serum", "virus"])
    )

filtered_metadata_df = metadata_df[
    metadata_df["serum"].isin(set(filtered_titers_df["serum"]))
].sort_values(["group", "serum"])

assert (
    filtered_metadata_df["serum"].nunique()
    == len(filtered_metadata_df)
    == filtered_titers_df["serum"].nunique()
)

print(f"\n{len(filtered_titers_df)=} titers for {len(filtered_metadata_df)=} sera")

print(f"Writing titers to {titers_csv=}")
filtered_titers_df.to_csv(titers_csv, index=False)

print(f"Writing sera metadata to {metadata_csv=}")
filtered_metadata_df.to_csv(metadata_csv, index=False)
