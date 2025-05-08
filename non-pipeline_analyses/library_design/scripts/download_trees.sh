#!/bin/bash
set -euo pipefail

# Download the trees

mkdir -p ./data/nextstrain_trees

curl https://nextstrain.org/seasonal-flu/h3n2/ha/6m@2025-03-28 \
    --header 'Accept: application/vnd.nextstrain.dataset.main+json' \
    --compressed > ./data/nextstrain_trees/h3n2_ha_6m.json

curl https://nextstrain.org/seasonal-flu/h1n1pdm/ha/6m@2025-03-28 \
    --header 'Accept: application/vnd.nextstrain.dataset.main+json' \
    --compressed > ./data/nextstrain_trees/h1n1pdm_ha_6m.json


# Extract relevant attributes from the tree JSON in a data frame, see https://gist.github.com/huddlej/5d7bd023d3807c698bd18c706974f2db for more details.

mkdir -p ./data/auspice_tables

python scripts/auspice_tree_to_table.py \
    --tree ./data/nextstrain_trees/h3n2_ha_6m.json \
    --output-metadata ./data/auspice_tables/h3n2_ha_6m.tsv \
    --attributes div num_date clade_membership subclade accession_ha country haplotype mutations lbi

python scripts/auspice_tree_to_table.py \
    --tree ./data/nextstrain_trees/h1n1pdm_ha_6m.json \
    --output-metadata ./data/auspice_tables/h1n1pdm_ha_6m.tsv \
    --attributes div num_date clade_membership subclade accession_ha country haplotype mutations lbi
