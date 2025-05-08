"""
get_metadata_from_headers.py

Description: Strip '|' from FASTA format headers and make a metadata file.
Author: Caroline Kikawa
"""

# Imports
from Bio import SeqIO
import pandas as pd

# Function
def fasta_headers_to_csv(fasta_file, output_csv, columns):
    # List to store the split headers
    headers = []

    # Parse the FASTA file
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Split the header on '|'
        splits = record.description.split('|')
        header_list.append(splits)

    # Create a DataFrame from the list of headers
    df = pd.DataFrame(header_list, columns = columns)

    # Save the DataFrame to a CSV file
    df.to_csv(output_csv, index=False, header=False)

# Main method
if __name__ == "__main__":
    fasta_file = snakemake.input.gisaid_fasta # Input FASTA file path
    columns = ['accession', 'strain', 'submitting_institute', 'date'] # Depends on header format
    output_csv = snakemake.output.metadata  # Output CSV file path
    
    fasta_headers_to_csv(fasta_file, output_csv)

