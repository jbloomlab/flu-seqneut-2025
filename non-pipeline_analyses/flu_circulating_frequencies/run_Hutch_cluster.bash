#!/bin/bash

echo "Running snakemake..."

snakemake \
    -j 16 \ # number of jobs
    --use-conda \
    --keep-going False \
    --rerun-incomplete \ # rerun if results is considered incomplete
    --rerun-triggers mtime \ # rerun if file modification date changes
    --retries 1 \ # number of attempts to retry if job fails
    --latency-wait 20 # seconds to wait before checking for missing files
    

echo "Done."