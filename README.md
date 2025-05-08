# Sequencing based neutralization assays using a library of pdmH1N1 and H3N2 human influenza strains
Experiments and analysis performed by Caroline Kikawa, using method and analysis developed by the [Bloom lab](https://jbloomlab.github.io/) and described in [Loes et al (2024)](https://journals.asm.org/doi/10.1128/jvi.00689-24) and [Kikawa et al (2025)](https://www.biorxiv.org/content/10.1101/2025.03.04.641544v1).

## Running the pipeline
This repository contains an analysis of the data using the Bloom lab software [`seqneut-pipeline`](https://github.com/jbloomlab/seqneut-pipeline) as a submodule. See that repository for intstructions on how to use Github submodules, including `seqneut-pipeline`. 

The configuration for the analysis is in [config.yml](config.yml) and the analysis itself is run by `snakemake` using [Snakefile](Snakefile).
Again, see [`seqneut-pipeline`](https://github.com/jbloomlab/seqneut-pipeline) for more description of how the pipeline works.

To run the pipeline, build the `seqneut-pipeline` conda environment from the [environment.yml](https://github.com/jbloomlab/seqneut-pipeline/blob/main/environment.yml) in `seqneut-pipeline`.
Then run the pipeline using:

    snakemake -j <n_jobs> --software-deployment-method conda

To run on the Hutch cluster, you can use the Bash script [run_Hutch_cluster.bash](run_Hutch_cluster.bash)
