# get_metadata.smk

# Description: Get sequence metadata from FASTA headers
# Author: Caroline Kikawa

rule merge_library_strains:
    """Add library strain ectodomain sequences to the reduced GISAID sequences."""
    input:
        gisaid_fasta = config["gisaid_prots"],
    output:
        metadata = config["gisaid_metadata"],
    script:
        "../scripts/get_metadata_from_headers.py"

