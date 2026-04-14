#!/usr/bin/env bash
#include <omp.h>
#include <stdio.h>
export OMP_NUM_THREADS=49
source ~/miniconda3/etc/profile.d/conda.sh

#I am requesting 4 nodes containing 49 CPUs, with 8 GB memory per CPU. Total: 392 GB
######
# This script performs the following:
# 1. Downloads NCBI taxonomy for the reference database. 

# 2. Downloads databases for taxonomic analysis. The databases are `bacteria`,
# `archaea`, `viral`, `human`, `UniVec_Core`, `fungi`, `plant`, `protozoa`,
# `plasmid`. 

# The output database is located in the `data/kraken2_db/k2_large_2025` directory

# 3. Adds the naked mole-rat reference genome to the database. The reference 
# genome path is `~/common_data/reference_genomes/Heter_glaber.v1.7_hic_pac_genomic_kraken2.fna`.

# 4. Builds the Kraken2 reference database in the same directory as the previous output
conda activate kraken2-tools-2.1.6

intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
## Download taxonomy
# echo "Downloading taxonomy"
# kraken2-build --download-taxonomy --db "${kraken2_db_dir}" --use-ftp
## Download databases
# intermediate_date_time=$(date +"%F %H:%M:%S")
# echo "${intermediate_date_time}"
# echo "Downloading bacteria"
# kraken2-build --download-library bacteria --db "${kraken2_db_dir}" --use-ftp

# intermediate_date_time=$(date +"%F %H:%M:%S")
# echo "${intermediate_date_time}"
# echo "Downloading archaea"
# kraken2-build --download-library archaea --db "${kraken2_db_dir}" --use-ftp

# intermediate_date_time=$(date +"%F %H:%M:%S")
# echo "${intermediate_date_time}"
# echo "Downloading viral"
# kraken2-build --download-library viral --db "${kraken2_db_dir}" --use-ftp

# intermediate_date_time=$(date +"%F %H:%M:%S")
# echo "${intermediate_date_time}"
# echo "Downloading human"
# kraken2-build --download-library human --db "${kraken2_db_dir}" --use-ftp

# intermediate_date_time=$(date +"%F %H:%M:%S")
# echo "${intermediate_date_time}"
# echo "Downloading UniVec_Core"
# kraken2-build --download-library UniVec_Core --db "${kraken2_db_dir}" --use-ftp

# intermediate_date_time=$(date +"%F %H:%M:%S")
# echo "${intermediate_date_time}"
# echo "Downloading fungi"
# kraken2-build --download-library fungi --db "${kraken2_db_dir}" --use-ftp

# intermediate_date_time=$(date +"%F %H:%M:%S")
# echo "${intermediate_date_time}"
# echo "Downloading plant"
# kraken2-build --download-library plant --db "${kraken2_db_dir}" --use-ftp

# intermediate_date_time=$(date +"%F %H:%M:%S")
# echo "${intermediate_date_time}"
# echo "Downloading protozoa"
# kraken2-build --download-library protozoa --db "${kraken2_db_dir}" --use-ftp

# intermediate_date_time=$(date +"%F %H:%M:%S")
# echo "${intermediate_date_time}"
# echo "Downloading plasmid"
# kraken2-build --download-library plasmid --db "${kraken2_db_dir}"

# intermediate_date_time=$(date +"%F %H:%M:%S")
# echo "${intermediate_date_time}"
# echo "Adding NMR"
# kraken2-build --add-to-library ~/common_data/reference_genomes/Heter_glaber.v1.7_hic_pac_genomic_kraken2.fna  \
#   --db "${kraken2_db_dir}"
# ### nt is for fragments of novel organisms (e.g. 16S rRNA gene) that don't have full genomes
# #kraken2-build --download-library nt --db "${kraken2_db_dir}"
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"

## Build the database: uses taxonomy and library
echo "Building the database"
kraken2-build --build --threads "${nthreads}" \
 --db "${kraken2_db_dir}" \
 --kmer-len "${kmer_len}"  \
 --minimizer-len "${minimizer_len}" \
 --minimizer-spaces "${minimizer_spaces}"
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
### 1. minimizer length must be no more than 31 for nucleotide databases, and 15 for protein databases
### 2. minimizer length must be no more than the k-mer length
### 3. if you change minimizer length, be sure to adjust minimizer spaces:
### A number s < l/4 can be chosen, and s positions in the minimizer will be
#  masked out during all comparisons.
echo "generating kmer distribution"
python ~/miniconda3/envs/kraken2-tools-2.1.6/bin/generate_kmer_distribution.py \
  -i "${kraken2_db_dir}"/database"${read_len}"mers.kraken \
  -o "${kraken2_db_dir}"/database"${read_len}"mers.kmer_distrib

intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
echo "Building the bracken index"
bracken-build -d "${kraken2_db_dir}" \
 -t "${nthreads}" \
 -k "${kmer_len}" \
 -l "${read_len}"
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
