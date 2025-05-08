# infer_phylogenetic_tree.smk

# Description: Aligns protein sequences and then builds phylogenetic tree
# Author: Caroline Kikawa





# Make alignments of HA1 and ectodomain
# Include library and all 2-year strains
# Include vaccine strains


rule infer_reduced_gisaid_protein_tree:
    """Use IQtree to infer GISAID HA tree"""
    input:
        aligned_reduced_gisaid_fasta = config["aligned_reduced_gisaid_prots"],
    output:
        config["reduced_gisaid_iqtree"],
    params:
        prefix = config["reduced_gisaid_iqtree_prefix"],
    shell:
        # 'FLU' is an amino acid substitution model for influenza proteins (https://pubmed.ncbi.nlm.nih.gov/20384985/)
        """
        iqtree \
            -s {input.aligned_reduced_gisaid_fasta} \
            -nt AUTO \
            -m FLU \
            -pre {params.prefix}
        """

rule reduced_timetree:
    """Run ``treetime`` to root and make time-scaled tree."""
    input:
        aligned_reduced_gisaid_fasta = config["aligned_reduced_gisaid_prots"],
        metadata = config['spikes_metadata'],
        reduced_tree = config["reduced_gisaid_iqtree"],
    output:
        os.path.splitext(config['divergencetree'])[0] + '.nexus',
        os.path.splitext(config['timetree'])[0] + '.nexus',
        config['root_to_tip'],
    params:
        outdir=os.path.dirname(config['timetree'])
    shell:
        """
        treetime \
            --aln {input.aligned_reduced_gisaid_fasta} \
            --dates {input.metadata} \
            --tree {input.reduced_tree} \
            --outdir {params.outdir} \
            --verbose 1
        """
