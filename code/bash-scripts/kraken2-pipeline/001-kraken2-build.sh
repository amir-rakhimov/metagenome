#!/usr/bin/env bash
export OMP_NUM_THREADS=49
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate kraken2-tools-2.1.6
set -euo pipefail
shopt -s nullglob
source config/bash/config.sh


# This script performs the following:
# 1. Downloads NCBI taxonomy for the reference database. 

# 2. Downloads databases for taxonomic analysis. The databases are `bacteria`,
# `archaea`, `viral`, `human`, `UniVec_Core`, `fungi`, `plant`, `protozoa`,
# `plasmid`. 

# The output database is located in the `data/kraken2_db/k2_large_2025` directory

# 3. Adds the naked mole-rat reference genome to the database. The reference 
# genome path is `~/common_data/reference_genomes/Heter_glaber.v1.7_hic_pac_genomic_kraken2.fna`.

# 4. Builds the Kraken2 reference database in the same directory as the previous output

# Function that marks a sample as done (in a hidden directory)
is_done()   { [[ -e "${mark_dir}/${2}.${1}.done" ]]; }   # is_done STEP SAMPLE
mark_done() { touch "${mark_dir}/${2}.${1}.done"; }      # mark_done STEP SAMPLE

echo "$(date +"%F %H:%M:%S")"
# Download taxonomy
if [[ -s "${kraken2_db_dir}/taxonomy/nucl_wgs.accession2taxid" ]]; then
	echo "Taxonomy exists";
else
	echo "Downloading taxonomy"
	kraken2-build --download-taxonomy --db "${kraken2_db_dir}" --use-ftp
fi

# Download databases
DB_NAMES=("bacteria" "archaea" "viral" "human" "UniVec_Core" "fungi" "plant" "plasmid" "protozoa")
for DB_NAME in  "${DB_NAMES[@]}"
do
	if [[ -s "${kraken2_db_dir}/library/${DB_NAME}/library.fna" ]]; then
		echo "${DB_NAME} database exists";	
	else
		echo "$(date +"%F %H:%M:%S")"
		echo "Downloading ${DB_NAME}"
		kraken2-build --download-library ${DB_NAME} --db "${kraken2_db_dir}" --use-ftp
	fi
done

# echo "$(date +"%F %H:%M:%S")"
# echo "Adding NMR"
# kraken2-build --add-to-library ~/common_data/reference_genomes/Heter_glaber.v1.7_hic_pac_genomic_kraken2.fna  \
#   --db "${kraken2_db_dir}"
# ### nt is for fragments of novel organisms (e.g. 16S rRNA gene) that don't have full genomes
# #kraken2-build --download-library nt --db "${kraken2_db_dir}"
# echo "$(date +"%F %H:%M:%S")"

## Build the database: uses taxonomy and library
echo "Building the database"
kraken2-build --build --threads "${nthreads_kraken2_build}" \
	--db "${kraken2_db_dir}" \
	--kmer-len "${kmer_len}"  \
	--minimizer-len "${minimizer_len}" \
	--minimizer-spaces "${minimizer_spaces}"
echo "$(date +"%F %H:%M:%S")"
### 1. minimizer length must be no more than 31 for nucleotide databases, and 15 for protein databases
### 2. minimizer length must be no more than the k-mer length
### 3. if you change minimizer length, be sure to adjust minimizer spaces:
### A number s < l/4 can be chosen, and s positions in the minimizer will be
#  masked out during all comparisons.
echo "generating kmer distribution"
python ~/miniconda3/envs/kraken2-tools-2.1.6/bin/generate_kmer_distribution.py \
	-i "${kraken2_db_dir}"/database"${read_len}"mers.kraken \
	-o "${kraken2_db_dir}"/database"${read_len}"mers.kmer_distrib

echo "$(date +"%F %H:%M:%S")"
echo "Building the bracken index"
bracken-build -d "${kraken2_db_dir}" \
	-t "${nthreads_kraken2_build}" \
	-k "${kmer_len}" \
	-l "${read_len}"
echo "$(date +"%F %H:%M:%S")"
