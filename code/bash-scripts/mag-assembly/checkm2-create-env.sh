#!/usr/bin/env bash
#SBATCH -t 10:00:00
#SBATCH -J 20250430-15-27-checkm2-create-env
#SBATCH --output jobreports/20250430-15-27-checkm2-create-env-out.txt
#SBATCH --error jobreports/20250430-15-27-checkm2-create-env-out.txt
source ~/miniconda3/etc/profile.d/conda.sh
# Install all tools with conda
## Setup channels for conda installation
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

## Create a conda environment
mamba env create -yqn checkm2-env -f code/bash-scripts/checkm2.yml
wget https://github.com/chklovski/CheckM2/archive/refs/heads/main.zip \
 -O ~/checkm2-main.zip

unzip ~/checkm2-main.zip -d ~/
conda activate checkm2-env
cd ~/CheckM2-main/
python setup.py install
# Download CheckM2 database
#checkm2 database --download --path ~/projects/metagenome/data/checkm2_db/

wget https://zenodo.org/records/14897628/files/checkm2_database.tar.gz \
 -O ~/projects/metagenome/data/checkm2_db/CheckM2_database/checkm2_database.tar.gz
tar -xvf ~/projects/metagenome/data/checkm2_db/CheckM2_database/checkm2_database.tar.gz \
 -C data/checkm2_db/CheckM2_database/
checkm2 testrun