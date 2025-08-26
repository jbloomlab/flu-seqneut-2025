# Sequencing based neutralization assays using a library of pdmH1N1 and H3N2 human influenza strains
Experiments and analysis performed by Caroline Kikawa, using method and analysis developed by the [Bloom lab](https://jbloomlab.github.io/) and described in [Loes et al (2024)](https://journals.asm.org/doi/10.1128/jvi.00689-24) and [Kikawa et al (2025)](https://www.biorxiv.org/content/10.1101/2025.03.04.641544v1).

## Summary of key data

### Viral strains

Details about the viruses included in the library are in [data/viral_libraries/flu-seqneut-2025-barcode-to-strain_actual.csv](data/viral_libraries/flu-seqneut-2025-barcode-to-strain_actual.csv).
That file contains an entry for each barcoded variant (note most strains have multiple barcodes) giving:
  - strain name
  - subtype
  - strain_type: either a strain representative of this circulating in 2025 designed for this study, or a past vaccine strain
  - barcode
  - nucleotide sequence for the HA **ectodomain only**
  - protein sequence for the HA **ectodomain only**
  - Genbank accession

Because the above file lists each barcode, strains are included more than once. If you want just FASTA files of the unique HA ectodomain sequences for the circulating 2025 strain and recent vaccine strains, see:
 - [results/viral_strain_seqs/circulating_2025_HA_ectodomain_prots.fa](results/viral_strain_seqs/circulating_2025_HA_ectodomain_prots.fa)
 - [results/viral_strain_seqs/circulating_2025_HA_ectodomain_nts.fa](results/viral_strain_seqs/circulating_2025_HA_ectodomain_nts.fa)

The above files contain the barcoded variants that actually ended up in the library at adequate representation to get high-quality data.
A small fraction of the strains that we attempted to design into the library did not actually end up at sufficient representation are so are not in the actual library; the full set of designed strains are in [data/viral_libraries/flu-seqneut-2025-barcode-to-strain_designed.csv](data/viral_libraries/flu-seqneut-2025-barcode-to-strain_designed.csv)

### Sera
Aggregated metadata about all relevant **human** sera for which titers were measured is in [results/aggregated_analyses/human_sera_metadata.csv](results/aggregated_analyses/human_sera_metadata.csv).
Note that these metadata are only for human sera not specified for exclusion (see relevant portion of [config.yml](config.yml) for details); in general the sera included here are the ones relevant to the goals of the study.
Note that the ages and collection dates are often ranges, so the numerical values in that CSV are often rounded to 1-5 year timespans for the ages and to the nearest month for the collection date.

The source data about the sera at the per-cohort (group) level and associated metadata are in [data/sera_metadata/](data/sera_metadata/); note these source data contain all sera that were planned for use in the experiments potentially including some for which no titers were measured, as well as for some non-human sera.
Specifically, the sera are from the following cohorts:
 - The Seattle Children's Hospital (`SCH`) cohort in Seattle, Washington, United States of America [data/sera_metadata/SCH_metadata.csv](data/sera_metadata/SCH_metadata.csv)
 - The University of Washington Medical Center (`UWMC`) cohort in Seattle, Washington, United States of America [data/sera_metadata/UWMC_metadata.csv](data/sera_metadata/UWMC_metadata.csv)
 - The National Institutes of Infectious Disease (`NIID`) cohort in Tokyo, Japan [data/sera_metadata/NIID_metadata.csv](data/sera_metadata/NIID_metadata.csv)
 - The EPI-HK cohort (`EPIHK`) at Hong Kong University in Hong Kong [data/sera_metadata/EPIHK_metadata.csv](data/sera_metadata/EPIHK_metadata.csv)
 - Innoculated ferrets from studies at the Francis Crick Institute (`FCI`) in London, United Kingdom [data/sera_metadata/FCI_metadata.csv](data/sera_metadata/FCI_metadata.csv)

### Titers
The aggregated titers for all relevant **human sera** are in [results/aggregated_analyses/human_sera_titers.csv](results/aggregated_analyses/human_sera_titers.csv).
Note that these titers are only for human sera not specified for exclusion (see relevant portion of [config.yml](config.yml) for details); in general the sera included here are the ones relevant to the goals of the study.

## Running the pipeline
This repository contains an analysis of the data using the Bloom lab software [seqneut-pipeline](https://github.com/jbloomlab/seqneut-pipeline) as a submodule. See that repository for instructions on how to use Github submodules, including `seqneut-pipeline`. 

The configuration for the analysis is in [config.yml](config.yml) and the analysis itself is run by `snakemake` using [Snakefile](Snakefile).
Input data are in [./data/](data), and all results created by the pipeline are placed in [./results/](results).

[Snakefile](Snakefile) both runs [seqneut-pipeline](https://github.com/jbloomlab/seqneut-pipeline) and runs custom rules within [Snakefile](Snakefile) that process the input data in [./data/](data) and the results in [./results/](results) to make some summary files and plots.

To run the pipeline, build the `seqneut-pipeline` conda environment from the [environment.yml](https://github.com/jbloomlab/seqneut-pipeline/blob/main/environment.yml) in `seqneut-pipeline`.
Then run the pipeline using:

    snakemake -j <n_jobs> --software-deployment-method conda

To run on the Hutch cluster, you can use the Bash script [run_Hutch_cluster.bash](run_Hutch_cluster.bash)
