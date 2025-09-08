# Near real-time on the human neutralizing antibody landscape to influenza virus to inform vaccine-strain selection in September 2025
This repository has the data and computer code for the study described in [Kikawa et al (2025)](https://doi.org/10.1101/2025.09.06.674661).

Briefly, this study used sequencing-based neutralization assays to measure titers to influenza viruses with HAs from seasonal H3N2 and H1N1 viruses representative of those circulating in the summer of 2025 against a set of human sera collected in late 2024 to spring of 2025.

## Summary of key data

### Viral strains

Details about the viruses included in the library are in [data/viral_libraries/flu-seqneut-2025-barcode-to-strain_actual.csv](data/viral_libraries/flu-seqneut-2025-barcode-to-strain_actual.csv).
That file contains an entry for each barcoded variant (note most strains have multiple barcodes) giving:
  - *strain*: strain name
  - *subtype*: H1N1 or H3N2
  - *strain_type*: either a strain representative of this circulating in 2025 designed for this study, or a past vaccine strain
  - *barcode*: nucleotide barcode
  - *nt_sequence_HA_ectodomain*: nucleotide sequence for the HA **ectodomain only**
  - *protein_sequence_HA_ectodomain*: protein sequence for the HA **ectodomain only**
  - *accession*: Genbank accession
  - *subclade*: subclade of strain
  - *collection_date*: collection date of strain
  - *vaccine_type*: if a vaccine strain, is it the cell- or egg-based strain

Because the above file lists each barcode, strains are included more than once. If you want just FASTA files of the unique HA ectodomain sequences for the circulating 2025 strain and recent vaccine strains, see:
 - [results/viral_strain_seqs/recent_HA_ectodomain_prots.fa](results/viral_strain_seqs/recent_HA_ectodomain_prots.fa)
 - [results/viral_strain_seqs/recent_HA_ectodomain_nts.fa](results/viral_strain_seqs/recent_HA_ectodomain_nts.fa)

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
 - Innoculated ferrets from studies at the Francis Crick Institute (`FCI`) in London, United Kingdom [data/sera_metadata/FCI_metadata.csv](data/sera_metadata/FCI_metadata.csv). Note that there were some issues with coagulation of these ferret sera that may affect those titers.

### Titers
The aggregated titers for all relevant **human sera** are in [results/aggregated_analyses/human_sera_titers.csv](results/aggregated_analyses/human_sera_titers.csv).
Note that these titers are only for human sera not specified for exclusion (see relevant portion of [config.yml](config.yml) for details); in general the sera included here are the ones relevant to the goals of the study.

### Interactive displays of results
See the HTML documentation rendered at [https://jbloomlab.github.io/flu-seqneut-2025/](https://jbloomlab.github.io/flu-seqneut-2025/) for interactive plots summarizing the results (at the bottom of the page), as well as notebooks showing all neutralization curves and details on per-plate and per-serum quality control.

## Running the pipeline
This repository contains an analysis of the data using the Bloom lab software [seqneut-pipeline](https://github.com/jbloomlab/seqneut-pipeline) as a submodule. See that repository for instructions on how to use Github submodules, including `seqneut-pipeline`. 

The configuration for the analysis is in [config.yml](config.yml) and the analysis itself is run by `snakemake` using [Snakefile](Snakefile).
Input data are in [./data/](data), and all results created by the pipeline are placed in [./results/](results).

[Snakefile](Snakefile) both runs [seqneut-pipeline](https://github.com/jbloomlab/seqneut-pipeline) and runs custom rules within [Snakefile](Snakefile) that process the input data in [./data/](data) and the results in [./results/](results) to make some summary files and plots.

To run the pipeline, build the `seqneut-pipeline` conda environment from the [environment.yml](https://github.com/jbloomlab/seqneut-pipeline/blob/main/environment.yml) in `seqneut-pipeline`.
Then run the pipeline using:

    snakemake -j <n_jobs> --software-deployment-method conda

To run on the Hutch cluster, you can use the Bash script [run_Hutch_cluster.bash](run_Hutch_cluster.bash)
