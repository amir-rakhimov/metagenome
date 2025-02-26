#!/usr/bin/env bash
source ~/miniconda3/etc/profile.d/conda.sh
# Install all tools with conda
## Setup channels for conda installation
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

## Create a conda environment
conda create -yqn mag_assembly-tools -c conda-forge -c bioconda megahit \
    prodigal metabat2 checkm2 gtdbtk busco
conda activate mag_assembly-tools
# Download CheckM2 database
checkm2 database --download --path ~/projects/metagenome/data/checkm2_db/
# Dowload GTDB-TK reference database
export GTDBTK_DATA_PATH=~/miniconda3/envs/mag_assembly-tools/share/gtdbtk-1.7.0/db/
download-db.sh

# Install QUAST
wget https://github.com/ablab/quast/releases/download/quast_5.2.0/quast-5.2.0.tar.gz -O ~/quast-5.2.0.tar.gz
tar -xzf ~/quast-5.2.0.tar.gz -C ~/