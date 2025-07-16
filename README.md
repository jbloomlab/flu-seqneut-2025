# Sequencing based neutralization assays using a library of pdmH1N1 and H3N2 human influenza strains
Experiments and analysis performed by Caroline Kikawa, using method and analysis developed by the [Bloom lab](https://jbloomlab.github.io/) and described in [Loes et al (2024)](https://journals.asm.org/doi/10.1128/jvi.00689-24) and [Kikawa et al (2025)](https://www.biorxiv.org/content/10.1101/2025.03.04.641544v1).

## Quick summary
* The viruses included in the library are outlined here:
    * Nucleotide and protein **HA ectodomain only** sequences in CSV format in [data/ha_sequences/flu-seqneut-2025-library.csv](data/ha_sequences/flu-seqneut-2025-library.csv)
    * Protein **HA ctodomain only** sequences FASTA format are placed in [data/ha_sequences/library_2025_HA_ectodomain_protein_sequences.fasta](data/ha_sequences/library_2025_HA_ectodomain_protein_sequences.fasta)
* The library design is outlined in detail in [non-pipeline_analyses/library_design/](non-pipeline_analyses/library_design/)
* The sera we will assay and their associated metadata are placed in [data/sera_metadata/](data/sera_metadata/)
    * The Seattle Children's Hospital (`SCH`) cohort in Seattle, Washington, United States of America [data/sera_metadata/SCH_metadata.csv](data/sera_metadata/SCH_metadata.csv)
    * The University of Washington Medical Center (`UWMC`) cohort in Seattle, Washington, United States of America [data/sera_metadata/UWMC_metadata.csv](data/sera_metadata/UWMC_metadata.csv)
    * The National Institutes of Infectious Disease (`NIID`) cohort in Tokyo, Japan [data/sera_metadata/NIID_metadata.csv](data/sera_metadata/NIID_metadata.csv)
    * The EPI-HK cohort (`EPIHK`) at Hong Kong University in Hong Kong [data/sera_metadata/EPIHK_metadata.csv](data/sera_metadata/EPIHK_metadata.csv)
    * Innoculated ferrets from studies at the Francis Crick Institute (`FCI`) in London, United Kingdom [data/sera_metadata/FCI_metadata.csv](data/sera_metadata/FCI_metadata.csv)

## Running the pipeline
This repository contains an analysis of the data using the Bloom lab software [`seqneut-pipeline`](https://github.com/jbloomlab/seqneut-pipeline) as a submodule. See that repository for intstructions on how to use Github submodules, including `seqneut-pipeline`. 

The configuration for the analysis is in [config.yml](config.yml) and the analysis itself is run by `snakemake` using [Snakefile](Snakefile).
Again, see [`seqneut-pipeline`](https://github.com/jbloomlab/seqneut-pipeline) for more description of how the pipeline works.

To run the pipeline, build the `seqneut-pipeline` conda environment from the [environment.yml](https://github.com/jbloomlab/seqneut-pipeline/blob/main/environment.yml) in `seqneut-pipeline`.
Then run the pipeline using:

    snakemake -j <n_jobs> --software-deployment-method conda

To run on the Hutch cluster, you can use the Bash script [run_Hutch_cluster.bash](run_Hutch_cluster.bash)
