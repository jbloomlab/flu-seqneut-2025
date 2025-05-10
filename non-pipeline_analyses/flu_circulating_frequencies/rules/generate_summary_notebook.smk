# generate_summary_notebook.smk

# Runs analysis notebook for protsets and saves all notebooks to results subdirectory
# Author: Caroline Kikawa


rule generate_summary_notebook:
    input:
        counts_by_date="results/strain_counts_{group}/{protset}_counts_by_date.csv"
    output:
        nb="results/notebooks/strain_counts_{group}_{protset}.ipynb"
    params:
        protset="{protset}",
        group="{group}"
    log:
        "results/notebooks/logs/strain_counts_{group}_{protset}.log"
    shell:
        """
        papermill notebooks/plot_strain_frequencies.ipynb {output.nb} \
            -p protset {params.protset} \
            -p group {params.group} \
            -p counts_by_date {input.counts_by_date} \
            > {log} 2>&1
        """



rule convert_notebook_to_html:
    input:
        notebook="results/notebooks/{notebook}.ipynb"
    output:
        html="/fh/fast/bloom_j/computational_notebooks/ckikawa/2025/flu-seqneut-2025/non-pipeline_analyses/flu_circulating_frequencies/results/notebooks/{notebook}.html"
    shell:
        "jupyter nbconvert --to html {input.notebook} --output {output.html}"
