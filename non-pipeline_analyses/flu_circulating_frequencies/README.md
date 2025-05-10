# Analyze circulating H3 HA protein sequences in 2022-2024
Analysis by Caroline Kikawa and Jesse Bloom. 

## Overview 
We want to be able to compare our selected library sequences to those circulating in nature.
To do this, we download all recent available sequneces from the database GISAID.
Then we compare protein sequences and identify exact and close matches. 

## Input data 
The configuration for the analysis is in [config.yaml](config.yaml) and the input data are in [./data/](data).

* [GISAID](https://gisaid.org/) sequences are in [data/](data/)
    * Downloaded from EpiFlu in April 2025 selecting the following options:
       - Influenza A
          - H3N2 or H1N1pdm09
       - Host of human
       - *Collection Date* from *April 1, 2024* to *April 15, 2025*
       - Select *Required Segments* HA and select *only complete* 
       - Only keeping *original* sequences (excluding lab passaged)
       - Only keeping complete sequences
       - Downloading just HA protein sequences
    * **Due to GISAID data sharing rules, this file is not tracked in the GitHub repository**

## Workflow
First, build and activate the conda environment with:

    conda env create -f environment.yml
    conda activate flu_circulating_frequencies

Configure the analysis in [config.yml](config.yml).
Then run the pipeline with:
        
    snakemake -j 16 --software-deployment-method conda

### Output
All the results are placed in [results](results) and are organized by the `group` flag designated in the analysis configuration. 
Analysis notebooks are rendered in [results/notebooks](results/notebooks).
We generally default to looking at the HA1 sequences that are within 1 amino acid mutation of currently circulating sequences.
The corresponding analysis notebooks for H1 and H3 are here:
* [results/notebooks/strain_counts_h1_h1-gisaid-ha1-within1.html](results/notebooks/strain_counts_h1_h1-gisaid-ha1-within1.html)
* [results/notebooks/strain_counts_h3_h3-gisaid-ha1-within1.html](results/notebooks/strain_counts_h3_h3-gisaid-ha1-within1.html)
    