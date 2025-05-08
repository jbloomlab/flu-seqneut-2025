"""
Get site-wise Shannon entropy from a protein alignment.

"""

# Imports
from Bio import SeqIO
import numpy as np
import pandas as pd
import collections


# Parameters
fasta_file = snakemake.input.aligned_gisaid_fasta
output_csv = snakemake.output.gisaid_shannon_entropy_csv
date_start = snakemake.params.date_start
date_end = snakemake.params.date_end

# Functions
def extract_date(header):
    """Extract the date from the FASTA header."""
    date = header.split('|')[3]  
    return date.strip()

def filter_sequences_by_date(fasta_file, date_start, date_end):
    """Filter sequences based on date range."""
    filtered_sequences = []
    
    # Ensure the FASTA file is read correctly
    with open(fasta_file, 'r') as handle:
        for record in SeqIO.parse(handle, "fasta"):
            date = extract_date(record.description)
            if date_start <= date <= date_end:
                filtered_sequences.append(str(record.seq))
    
    return filtered_sequences

def get_string(val):
    str_val = str(val)
    return str_val

def get_int(val):
    int_val = int(val)
    return int_val

def calculate_shannon_entropy(sequences): 
    """ 
    Calculate Shannon entropy for each site in alignment. 
    Here, calculating Shannon entropy as: 
    H = -sum(pk * log(pk))
    """
    # Number of sites based on the length of the first sequence
    num_sites = len(sequences[0]) if sequences else 0
    entropy_list = []
    
    for i in range(num_sites):
        # Get the residues at this site across all sequences
        residues = [seq[i] for seq in sequences if i < len(seq)]
        
        # Calculate frequency of each residue
        if residues:
            unique, counts = np.unique(residues, return_counts=True)
            
            # Calculate Shannon entropy
            probabilities = counts / len(residues)
            entropy = -np.sum(probabilities * np.log2(probabilities))

            # Get dictionary of residues
            temp_dict = dict(zip(unique, counts))
            # Convert to regular Python types for both keys and values
            temp_dict = {str(key) if isinstance(key, (np.str_, np.integer, np.floating)) else key: 
                              (str(value) if isinstance(value, np.str_) else
                               int(value) if isinstance(value, np.integer) else
                               float(value) if isinstance(value, np.floating) else value)
                              for key, value in temp_dict.items()}

            # Get consensus residue
            consensus_res = max(temp_dict, key=temp_dict.get)

            # Save entropy, consensus residue and residue/counts dictionary
            entropy_list.append([entropy, consensus_res, temp_dict])
            
        else:
            # Save entropy, consensus residue and residue/counts dictionary
            entropy_list.append([0, np.Nan, np.Nan])  # If no residues, entropy is 0
    
    return entropy_list


# Process
filtered_sequences = filter_sequences_by_date(fasta_file, date_start, date_end)
entropy_values = calculate_shannon_entropy(filtered_sequences)

# Create a DataFrame and save results to CSV
entropy_df = (pd.DataFrame(entropy_values, columns = ['entropy', 'consensus_res', 'unique_res_counts'])
              # Add site column to DataFrame
              # .insert(0, 'site', list(range(1, len(entropy_values) + 1)))  # Sites are 1-indexed for easier interpretation
             )

# Add site column to DataFrame
entropy_df['site'] = list(range(1, len(entropy_values) + 1))  # Sites are 1-indexed for easier interpretation

# Reorder columns
entropy_df = entropy_df[['site', 'entropy', 'consensus_res', 'unique_res_counts']]

# Save DataFrame to CSV
entropy_df.to_csv(output_csv, index=False)