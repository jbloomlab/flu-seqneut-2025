# calculate_shannon_entropy.smk

# Description: Calculates Shannon entropy within date range
# Author: Caroline Kikawa

rule calculate_shannon_entropy:
    """
    Get site-wise Shannon entropy from alignment.
    Subselect sequences based on collection date.
    Specify dates in params.
    """
    input:
#        aligned_ncbi_fasta = config["aligned_ncbi_prots"],
        aligned_gisaid_fasta = config["aligned_gisaid_prots"],
    output:
#        ncbi_shannon_entropy_csv = config["shannon_entropy_ncbi_prots"],
        gisaid_shannon_entropy_csv = config["shannon_entropy_gisaid_prots"],
    params:
        date_start = "2023-05-01",
        date_end = "2023-11-01",
    script:
        "../scripts/calculate_sitewise_entropy.py"
