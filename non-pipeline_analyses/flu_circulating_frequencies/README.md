# Analyze circulating H3 HA protein sequences in 2022-2024
Analysis by Caroline Kikawa and Jesse Bloom. 

## Overview 

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

Run with:
`snakemake -j 16 --software-deployment-method conda`

### Get frequencies of strains matching (or nearly matching) library strains

###
    