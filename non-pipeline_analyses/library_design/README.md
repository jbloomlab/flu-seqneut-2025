# Choosing HA sequences and designing barcoded constructs for the `seqneut` library
This subdirectory describes the creation of the 2025 pdmH1N1 and H3N2 combined subtype library, including details on the barcoded HA expression plasmid and the rationale for the strains chosen for the library. 
Author: Caroline Kikawa

## Quick summary
This directory contains the analysis performed to select strains and the design of the barcoded HA constructs for the 2025 pdmH1N1 and H3N2 library. 

The rest of this README contains the **overview (*section 1*)** and the details of **strain choice (*section 2*)**, **barcoded construct design and ordering (*section 3*)**, and **Nextstrain tree generation (*section 4*)** involved in designing a library for `seqneut` assays.

## 1. Library design overview
There were several approaches to subselect pdmH1N1 and H3N2 strains for the HA variant libraries, keeping in mind an overall library size is limited to ~100 HA variants. 
* Select high "branchiness" strains listed in [Nextstrain](https://nextstrain.org/) 6-month builds.
* Antigenically drifted strains (c/o John Huddleston from the Bedford lab and Andrea Loes from the Bloom lab)
* For H3's, strains with mutations at higher observations than expected based on site-specific mutation rates, and that are at antigenic sites (c/o Sam Turner from the Smith lab)
* Also strains used for ferret infections.

## 2. Choosing strains -- the detailed approach

### 2.1. Download and process the Nexstrain 6-month pdmH1N1 and H3N2 HA trees to get accessions of interest
Build and activate the conda environment to perform these analysis placed here [environments/nextstrain.yaml](environments/nextstrain.yaml) with the following commands:

        conda env create -f environments/nextstrain.yaml
        conda activate nextstrain

Download the trees and extract relevant attributes with the bash script (note that for reproducibility [scripts/download_trees.sh](scripts/download_trees.sh) specifies the date of the JSON using the `@2025-03-28` type syntax at the end of the URLs, update this if you want newer trees):

        bash scripts/download_trees.sh

These commands create the files in [data/nextstrain_trees](data/nextstrain_trees) and [data/auspice_tables](data/auspice_tables).

Next, select nodes with high "branchiness" defined by the Local Branching Index (LBI), and compare accessions to those suggested by collaborators, and save a non-redundant list of all strain accessions to [data/accessions_to_download/](data/accessions_to_download/). The strains suggested by collaborators are in [data/strain_lists_from_collaborators](data/strain_lists_from_collaborators), and were selected by various means. The script run by the command below combines those strains with the ones with high branchiness to to create the list of accessions to download. Note that the script has various file paths and parameters (in particular, the cutoff dates for the haplotypes and the cutoff for the number of strains to choose) hardcoded in the script, so those parameters need to be adjusted if you want to do things differently. Use this command:

        python scripts/select_compare_high_lbi_sequences.py

In addition, the following files were manually added to [data/accessions_to_download/](data/accessions_to_download/):

 - [data/accessions_to_download/ST_strains.csv](data/accessions_to_download/ST_strains.csv): accessions suggested by Sam Turner
 - [data/accessions_to_download/H1_ferret_strains.csv](data/accessions_to_download/H1_ferret_strains.csv): H1 strains used for ferret infections by Nicola Lewis.
 - [data/accessions_to_download/H3_ferret_strains.csv](data/accessions_to_download/H3_ferret_strains.csv): H3 strains used for ferret infections by Nicola Lewis.
 - [data/accessions_to_download/2025-04-25_additional_H1s.csv](data/accessions_to_download/2025-04-25_additional_H1s.csv): Additional H1 strains representing haplotypes present in April 2025 that were not noted in the original analyses done by Caroline Kikawa, Andrea Loes and John Huddleston in February-March 2025.

### 2.2. Download and process sequences from GISAID
Next, download from GISAID the sequences of interest using their accessions as compiled by the previous step in [data/accessions_to_download/](data/accessions_to_download/).
To do this:
  * Navigate to the EpiFlu tab at [GISAID](https://gisaid.org/)
  * Paste tab-delimited accessions and search to get the sequences to download.
  * The downloaded sequences are in [./data/sequences](./data/sequences)
  * **Due to GISAID data sharing rules, these files are not tracked in the GitHub repository**

In addition, the following files were manually added to [data/sequences/](data/sequences/):
- [data/sequences/h1_id_map.csv](data/sequences/h1_id_map.csv): the matched isolate GISAID accession and HA gene GISAID accession for each H1 strain, necessary for mapping strain IDs provided by collaborators
- [data/sequences/h3_id_map.csv](data/sequences/h3_id_map.csv): the matched isolate GISAID accession and HA gene GISAID accession for each H3 strain, necessary for mapping strain IDs provided by collaborators

Align the sequences with the bash script:

        bash scripts/align_sequences.sh

The alignments go into `results/gisaid_circulating_alignments/`, which is not tracked in this repo due to GISAID data sharing rules.

### 2.3. Run interactive Jupyter Notebook to choose final protein sequences
First, build and activate the conda environment placed in [environments/library_design.yaml](environments/library_design.yaml) with:

        conda env create -f environments/library_design.yaml
        conda activate library_design
    
Then, run the Jupyter Notebook [notebooks/design_base_library.ipynb](notebooks/design_base_library.ipynb) interactively to analyze all Nextstrain and Nextclade haplotypes together. 
That notebook is annotated with details on how we finalized the set of sequences included in the library. 

The output of the notebook is placed in `results/selected_library_strains/`, and includes FASTA files of the selected H1 and H3 proteins and CSVs summarizing those proteins.
The key output file is `results/selected_library_strains/h3_and_h1_prots.fasta`, which contains the H1 and H3 proteins.
Due to GISAID data sharing rules, these files are not tracked in the repo.

### 2.4. Get Genbank accessions / sequences corresponding to each sequence.
The next step is to get Genbank accessions and nucleotide sequences corresponding to the protein for each strain listed in `results/selected_library_strains/h3_and_h1_prots.fasta`.
We do this as Genbank sequences (unlike GISAID sequences) can be publicly shared.
To do this, we use the script at [https://github.com/jbloom/match_prot_to_genbank_nt](https://github.com/jbloom/match_prot_to_genbank_nt).
Specifically:

First add [https://github.com/jbloom/match_prot_to_genbank_nt](https://github.com/jbloom/match_prot_to_genbank_nt) as a git submodule with `git submodule add https://github.com/jbloom/match_prot_to_genbank_nt`.

Then build the conda environment in [match_prot_to_genbank_nt/environment.yaml](match_prot_to_genbank_nt/environment.yaml) with:

    conda env create -f match_prot_to_genbank_nt/environment.yaml

Then activate that conda environment with:

    conda activate match_prot_to_genbank_nt

Now use the script to find the matches with:

    python match_prot_to_genbank_nt/match_prot_to_genbank_nt.py \
        --query-prots results/selected_library_strains/h3_and_h1_prots.fasta \
        --taxon "Influenza A virus" \
        --outdir results/strains_for_library

This creates the file [results/strains_for_library/match_prot_to_genbank_nt.csv](results/strains_for_library/match_prot_to_genbank_nt.csv) which contains the matches for all of the proteins.

To create a FASTA of the nucleotide sequences, then run:

    python scripts/matches_to_fasta.py

This creates from [results/strains_for_library/match_prot_to_genbank_nt.csv](results/strains_for_library/match_prot_to_genbank_nt.csv) the following output files:

 - [results/strains_for_library/h3_nt_seqs_for_library.fasta](results/strains_for_library/h3_nt_seqs_for_library.fasta): the H3 sequences
 - [results/strains_for_library/h1_nt_seqs_for_library.fasta](results/strains_for_library/h1_nt_seqs_for_library.fasta): the H1 sequences
 - [results/strains_for_library/h3_and_h1_nt_seqs_for_library.fasta](results/strains_for_library/h3_and_h1_nt_seqs_for_library.fasta): the H3 and H1 sequences

## 3. Barcoded construct design

The barcoded construct used in these experiments is similar to those described in [Loes et al (2024)](https://journals.asm.org/doi/10.1128/jvi.00689-24) and [Kikawa et al. 2025](https://www.biorxiv.org/content/10.1101/2025.03.04.641544v1).

Constructs are designed in an interactive notebook [notebooks/design_HA_inserts.ipynb](notebooks/design_HA_inserts.ipynb). 
First, reactivate the `library_design` conda environment with:

    conda activate library_design

Then run [notebooks/design_HA_inserts.ipynb](notebooks/design_HA_inserts.ipynb).

The final output of the above notebook summarized here:
* A spreadsheet [results/ordersheets/ordersheet.csv](results/ordersheets/ordersheet.csv) containing tables including columns `name` and `sequence` required for ordering constructs from Twist
* Spreadsheets of [results/ordersheets/h1_inserts.csv](results/ordersheets/h1_inserts.csv) and [results/ordersheets/h3_inserts.csv](results/ordersheets/h3_inserts.csv) which include the abbreviated `name` for Twist ordering, the `strain_name` and the `sequence` for each construct
