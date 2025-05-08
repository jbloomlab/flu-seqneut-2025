# align_fastas.smk

# Description: Aligns FASTA format protein sequences
# Author: Caroline Kikawa





rule reduce_gisaid_protein_sequences:
    """Reduce number of similar sequences in 2-year GISAID download."""
    input:
        gisaid_fasta = config["gisaid_prots"],
    output:
        reduced_gisaid_fasta = config["reduced_gisaid_prots"],
    shell:
        # Set sequence identity threshold flag -c to 99%
        "cd-hit -i {input.gisaid_fasta} -o {output.reduced_gisaid_fasta} -c 0.99"


rule merge_library_strains:
    """Add library strain ectodomain sequences to the reduced GISAID sequences."""
    input:
        reduced_gisaid_fasta = config["reduced_gisaid_prots"],
        library_chimeras_fasta = config["strain_prots"],
    output:
        reduced_gisaid_plus_library_fasta = config["reduced_gisaid_library_prots"],
    script:
        "../scripts/merge_fasta_files.py"


rule align_reduced_gisaid_protein_sequences:
    """Align reduced HA sequences from GISAID with MAFFT"""
    input:
        reduced_gisaid_plus_library_fasta = config["reduced_gisaid_library_prots"],
    output:
        aligned_reduced_gisaid_plus_library_fasta = config["aligned_reduced_gisaid_library_prots"],
    shell:
        "mafft --auto {input.reduced_gisaid_plus_library_fasta} > {output.aligned_reduced_gisaid_plus_library_fasta}"

