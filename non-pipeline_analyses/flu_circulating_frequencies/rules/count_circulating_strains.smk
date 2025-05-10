# count_circulating_strains.smk

# Counts circulating strains that match library strains
# Author: Jesse Bloom, edited by Caroline Kikawa

rule count_circulating_strains:
    """Get counts of each strain."""
    input:
        strain_prots=lambda wc: config["strain_prots"][wc.group],
        protset=lambda wc: config["protsets"][wc.protset]["protset"],
    output:
        counts_overall="results/strain_counts_{group}/{protset}_counts_overall.csv",
        counts_by_date="results/strain_counts_{group}/{protset}_counts_by_date.csv",
        strain_matches="results/strain_counts_{group}/{protset}_strain_matches.csv",
        counts_unmatched_sequences="results/strain_counts_{group}/{protset}_unmatched_sequence_counts.csv",
        counts_recent_unmatched_sequences="results/strain_counts_{group}/{protset}_counts_recent_unmatched_sequences.csv",
    params:
        trim=lambda wc: config["protsets"][wc.protset]["trim"],
        maxdiff=lambda wc: config["protsets"][wc.protset]["maxdiff"],
    log:
        "results/strain_counts_{group}/logs/strain_counts_{protset}.txt",
    script:
        "../scripts/strain_counts.py"