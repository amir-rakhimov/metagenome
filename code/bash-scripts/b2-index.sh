#!/usr/bin/env bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate qc-tools
bowtie2-build --large-index --threads 12 common_data/reference_genomes/all_hosts_reference.fasta \
common_data/reference_genomes/bowtie2_indices/all_hosts_reference