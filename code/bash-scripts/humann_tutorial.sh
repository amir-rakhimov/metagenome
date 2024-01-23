# Set conda channel priority
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels biobakery

# Create a conda environment
conda create --name humann-tools-3.7 -c biobakery python=3.7 humann metaphlan

conda activate humann-tools-3.7 
# Testing installation: run HUMAnN unit tests
humann_test
# Download demo databases
humann_databases --download chocophlan DEMO humann_dbs
humann_databases --download uniref DEMO_diamond humann_dbs
# HUMAnN uses MetaPhlAn to detect organisms in the community. It should be installed
metaphlan --version

# 1. Metagenome functional profiling
mkdir data 
mkdir output
cd data
wget https://github.com/biobakery/humann/raw/master/examples/demo.fastq.gz
wget https://github.com/biobakery/humann/raw/master/examples/demo.m8
wget https://github.com/biobakery/humann/raw/master/examples/demo.sam
cd ..
# Running HUMAnN
humann --input data/demo.fastq.gz --output output/demo_fastq --threads 8