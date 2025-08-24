#!/bin/bash

echo "Running snakemake..."

snakemake -j 16 --software-deployment-method conda

echo "Done."
