"""Top-level ``snakemake`` file that runs analysis."""

import pandas as pd
from os.path import join

configfile: "config.yml"

include: "seqneut-pipeline/seqneut-pipeline.smk"


rule all:
    input:
        # output from seqneut-pipeline
        seqneut_pipeline_outputs,
        # outputs from custom analyses
        "results/viral_strain_seqs/circulating_2025_HA_ectodomain_nts.fa",
        "results/viral_strain_seqs/circulating_2025_HA_ectodomain_prots.fa",
        "results/aggregated_analyses/human_sera_metadata.csv",    

# --- Everything below here is a custom analysis not part of seqneut-pipeline -----------


rule recent_strains_fasta:
    """Generate FASTA of recent strains in library."""
    input:
        csv=config["viral_libraries"]["flu-seqneut-2025_library_actual"],
    output:
        nt_fasta="results/viral_strain_seqs/circulating_2025_HA_ectodomain_nts.fa",
        prot_fasta="results/viral_strain_seqs/circulating_2025_HA_ectodomain_prots.fa",
    params:
        recent_vaccine_strains=config["recent_vaccine_strains"],  # also include these
    conda:
        "seqneut-pipeline/environment.yml"
    log:
        "results/logs/recent_strains_fasta.txt",
    script:
        "scripts/recent_strains_fasta.py"


# human sera groups to use in aggregated analyses of human sera titers
human_sera_groups_to_analyze = [
    g for g in groups if g not in config["human_sera_groups_to_exclude"]
]

rule aggregate_human_sera_metadata_and_titers:
    """Aggregate metadata and titers on human sera."""
    input:
        groups_metadata=expand(
            "data/sera_metadata/{group}_metadata.csv", group=human_sera_groups_to_analyze
        ),
        group_titers=expand(
            "results/aggregated_titers/titers_{group}.csv",
            group=human_sera_groups_to_analyze,
        ),
    params:
        groups=human_sera_groups_to_analyze,
        human_sera_to_exclude=config["human_sera_to_exclude"],
    output:
        metadata_csv="results/aggregated_analyses/human_sera_metadata.csv",
        titers_csv="results/aggregated_analyses/human_sera_titers.csv",
    conda:
        "seqneut-pipeline/environment.yml"
    log:
        notebook="results/aggregated_analyses/aggregate_human_sera_metadata_and_titers.ipynb",
    notebook:
        "notebooks/aggregate_human_sera_metadata_and_titers.py.ipynb"
