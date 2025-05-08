''' 
select_compare_high_lbi_sequences.py
Author: Caroline Kikawa
'''

##################################################################################################################################
# Import packages
import pandas as pd
import os

accessionsdir = './data/accessions_to_download/'
os.makedirs(accessionsdir, exist_ok=True)

##################################################################################################################################
# Get high frequency haplotypes
h3n2_haploytypes = pd.read_csv('./data/auspice_tables/h3n2_ha_6m.tsv', sep='\t')
h1n1pdm_haplotypes = pd.read_csv('./data/auspice_tables/h1n1pdm_ha_6m.tsv', sep='\t')

# Only choose haplotypes sampled within last 12 months
# Today is April 2, 2025, or 2025.257, so 12 months ago is 2024.257
h3n2_representative_strains = (h3n2_haploytypes.query('num_date > 2024.257'))

# The 'N145S' specified in the J.2 strains is unecessary and weird
h3n2_representative_strains['haplotype'] = h3n2_representative_strains['haplotype'].apply(lambda x: x.replace('-N145S', '') if 'J.2:' in x else x)

# Select strains with lowest branch divergence for each haplotypes
h3n2_representative_strains = h3n2_representative_strains.loc[h3n2_representative_strains.groupby("haplotype")["div"].idxmin()]

# Choose top N strains with highest LBI
N_strains_cutoff = 20
h3n2_representative_strains = (h3n2_representative_strains
                               .sort_values(by = 'lbi', ascending=False)
                               .head(N_strains_cutoff)
                               .reset_index(drop=True)
                              )

# Only choose haplotypes sampled within last 6 months
# Today is April 2, 2025, or 2025.257, so 6 months ago is 2024.757
h1n1pdm_representative_strains = (h1n1pdm_haplotypes.query('num_date > 2024.757'))

# Select strains with lowest branch divergence for each haplotypes
h1n1pdm_representative_strains = h1n1pdm_representative_strains.loc[h1n1pdm_representative_strains.groupby("haplotype")["div"].idxmin()]

# Choose top N strains with highest LBI
N_strains_cutoff = 10
h1n1pdm_representative_strains = (h1n1pdm_representative_strains
                               .sort_values(by = 'lbi', ascending=False)
                               .head(N_strains_cutoff)
                               .reset_index(drop=True)
                              )


##################################################################################################################################
## Strains from collaborators
AL_H1_strains = pd.read_csv('./data/strain_lists_from_collaborators/250328_H1_selected_sequences_AL.csv')
AL_H3_strains = pd.read_csv('./data/strain_lists_from_collaborators/250328_H3_selected_sequences_AL.csv')

JH_H1_strains = pd.read_csv('./data/strain_lists_from_collaborators/Huddleston-strains-for-neutralization_H1.csv')
JH_H3_strains = pd.read_csv('./data/strain_lists_from_collaborators/Huddleston-strains-for-neutralization_H3.csv')

##################################################################################################################################
# Compile all H3 strains
h3n2_representative_strains['method'] = 'CK'
AL_H3_strains['method'] = 'AL'
JH_H3_strains['method'] = 'JH'

all_h3 = pd.concat([
    h3n2_representative_strains,
    AL_H3_strains,
    JH_H3_strains
])


outfile = os.path.join(accessionsdir, 'H3_accessions.tsv')
print(f'Saving accessions to {outfile}...')
all_h3[['accession_ha']].drop_duplicates().to_csv(outfile, index=False, sep='\t')
print('Done.')

outfile = os.path.join(accessionsdir, 'H3_candidate_summary.csv')
print(f'Saving strain info dataframe to {outfile}...')
all_h3.to_csv(outfile, index=False)
print('Done.')

##################################################################################################################################
# Compile all H1 strains
h1n1pdm_representative_strains['method'] = 'CK'
AL_H1_strains['method'] = 'AL'
JH_H1_strains['method'] = 'JH'

all_h1 = pd.concat([
    h1n1pdm_representative_strains,
    AL_H1_strains,
    JH_H1_strains
])

outfile = os.path.join(accessionsdir, 'H1_accessions.tsv')
print(f'Saving accessions to {outfile}...')
all_h1[['accession_ha']].drop_duplicates().to_csv(outfile, index=False, sep='\t')
print('Done.')

outfile = os.path.join(accessionsdir, 'H1_candidate_summary.csv')
print(f'Saving strain info dataframe to {outfile}...')
all_h1.to_csv(outfile, index=False)
print('Done.')
