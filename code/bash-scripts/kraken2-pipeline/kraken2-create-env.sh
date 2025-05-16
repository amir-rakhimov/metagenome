#!/usr/bin/env bash
source ~/miniconda3/etc/profile.d/conda.sh
# Install all tools with conda
## Setup channels for conda installation
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# Delete the previous env
conda remove --name kraken2-tools-2.1.3 --all

## Create a conda environment
conda create -yqn kraken2-tools-2.1.3 -c conda-forge -c bioconda kraken2 krakentools bracken krona bowtie2

conda activate kraken2-tools-2.1.3