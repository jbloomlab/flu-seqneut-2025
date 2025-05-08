#!/bin/bash
set -euo pipefail

# Align downloaded sequences

mkdir -p ./results/alignments
mkdir -p ./results/gisaid_circulating_alignments

# Function to run MAFFT only if output file doesn't exist
run_mafft_if_needed() {
    local input_file=$1
    local output_file=$2

    if [[ ! -f "$output_file" ]]; then
        echo "Aligning $input_file -> $output_file"
        mafft --auto "$input_file" > "$output_file"
    else
        echo "Alignment already exists: $output_file"
    fi
}

run_mafft_if_needed ./data/sequences/2025-04-07_h1n1pdm_ha.fasta ./results/alignments/2025-04-07_h1n1pdm_ha_aligned.fasta
run_mafft_if_needed ./data/sequences/2025-04-07_h3n2_ha.fasta ./results/alignments/2025-04-07_h3n2_ha_aligned.fasta
run_mafft_if_needed ./data/sequences/2025-04-08_ferret_strains_h3.fasta ./results/alignments/2025-04-08_ferret_strains_h3.fasta
run_mafft_if_needed ./data/sequences/2025-04-08_ferret_strains_h1.fasta ./results/alignments/2025-04-08_ferret_strains_h1.fasta
run_mafft_if_needed ./data/sequences/ST_h3n2_metadata.fasta ./results/alignments/ST_h3n2_metadata.fasta
run_mafft_if_needed ./data/sequences/2025-04-25_additional_H1s.fasta ./results/alignments/2025-04-25_additional_H1s.fasta
