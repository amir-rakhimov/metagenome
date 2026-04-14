#!/usr/bin/env bash
RUN_ID=$(date +%Y%m%d_%H%M%S)
# Kraken2 parameters:
######
# kmer length for kraken2-build and braken-build
kmer_len=35 
minimizer_len=31
minimizer_spaces=7
# read length for kraken2-build and bracken
read_len=150
nthreads=38
# kraken2_db_dir=data/kraken2_db/k2_large_20240307
# kraken2_db_dir=data/kraken2_db/k2_large_${date_var} 
kraken2_db_dir=data/kraken2_db/k2_large_20250919

