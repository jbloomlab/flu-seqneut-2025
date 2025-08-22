"""Top-level ``snakemake`` file that runs analysis."""

import pandas as pd
from os.path import join

configfile: "config.yml"

include: "seqneut-pipeline/seqneut-pipeline.smk"


rule all:
    input:
        seqneut_pipeline_outputs,  # outputs from pipeline
        "results/viral_strain_seqs/circulating_2025_HA_ectodomain_nts.fa",
        "results/viral_strain_seqs/circulating_2025_HA_ectodomain_prots.fa",


rule recent_strains_fasta:
    """Generate FASTA of recent strains in library."""
    input:
        csv=config["viral_libraries"]["flu-seqneut-2025_library_actual"],
    output:
        nt_fasta="results/viral_strain_seqs/circulating_2025_HA_ectodomain_nts.fa",
        prot_fasta="results/viral_strain_seqs/circulating_2025_HA_ectodomain_prots.fa",
    conda:
        "seqneut-pipeline/environment.yml"
    log:
        "results/logs/recent_strains_fasta.txt",
    script:
        "scripts/recent_strains_fasta.py"
