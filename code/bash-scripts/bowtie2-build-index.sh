#!/usr/bin/env bash
nthreads=12
source ~/miniconda3/etc/profile.d/conda.sh
conda activate qc-tools
bowtie2-build --large-index --threads ${nthreads} ~/common_data/reference_genomes/all_hosts_reference.fasta \
  ~/common_data/bowtie2_indices/all_hosts_reference

