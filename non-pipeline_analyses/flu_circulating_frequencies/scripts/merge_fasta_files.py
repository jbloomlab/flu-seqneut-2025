"""
merge_fasta_files.py

Description: Takes library and GISAID FASTA files, compares them, and makes output FASTA of
sequences from both files, keeping the library sequence over the GISAID sequence if they are 
an exact match. 
Author: Caroline Kikawa
"""

# Imports
from Bio import SeqIO

# # Functions
# def get_sequences_from_fasta(fasta_file):
#     sequences = {}
#     for record in SeqIO.parse(fasta_file, "fasta"):
#         sequences[str(record.seq)] = record
#     return sequences

# def merge_fasta_files(fasta_a, fasta_b, output_file):
#     """Merges sequences from two FASTA files, prioritizing sequences from library_fasta if they are in both files."""
#     # Load sequences from both FASTA files
#     library_seq = get_sequences_from_fasta(library_fasta)
#     fasta_seq = get_sequences_from_fasta(gisaid_fasta)

#     # Merge sequences from both library and GISAID, prioritizing library
#     combined_seq = list(library_seq.values())  # Start with sequences from library
#     for seq, record in gisaid_fasta.items():
#         if seq not in library_seq:  # Add GISIAD sequence only if sequence is not already in library
#             combined_seq.append(record)

#     # Write the combined sequences to output.fasta
#     with open(output_fasta, "w") as output_handle:
#         SeqIO.write(combined_seq, output_handle, "fasta")

# # Main method
# if __name__ == "__main__":
#     library_fasta = snakemake.input.library_chimeras_fasta  # Input library FASTA
#     gisaid_fasta = snakemake.input.reduced_gisaid_fasta  # Input GISAID FASTA
#     output_fasta = snakemake.output.reduced_gisaid_plus_library_fasta  # Output FASTA

#     merge_fasta_files(library_fasta, gisaid_fasta, output_fasta)


# Functions
def get_sequences_from_fasta(fasta_file):
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequences[str(record.seq)] = record
    return sequences

def merge_fasta_files(library_fasta, gisaid_fasta, output_file):
    """Merges sequences from two FASTA files, prioritizing sequences from library_fasta if they are in both files."""
    # Load sequences from both FASTA files
    library_seq = get_sequences_from_fasta(library_fasta)
    fasta_seq = get_sequences_from_fasta(gisaid_fasta)  # Correctly load GISAID sequences into fasta_seq

    # Merge sequences from both library and GISAID, prioritizing library
    combined_seq = list(library_seq.values())  # Start with sequences from library
    for seq, record in fasta_seq.items():  # Iterate over fasta_seq, not gisaid_fasta
        if seq not in library_seq:  # Add GISAID sequence only if sequence is not already in library
            combined_seq.append(record)

    # Write the combined sequences to output.fasta
    with open(output_file, "w") as output_handle:
        SeqIO.write(combined_seq, output_handle, "fasta")

# Main method
if __name__ == "__main__":
    library_fasta = snakemake.input.library_chimeras_fasta  # Input library FASTA
    gisaid_fasta = snakemake.input.reduced_gisaid_fasta  # Input GISAID FASTA
    output_fasta = snakemake.output.reduced_gisaid_plus_library_fasta  # Output FASTA

    merge_fasta_files(library_fasta, gisaid_fasta, output_fasta)