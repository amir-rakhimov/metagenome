#!/usr/bin/env bash
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
######
date_var=$(date -I|sed 's/-//g')
# kmer length for kraken2-build and braken-build
kmer_len=35 
minimizer_len=31
minimizer_spaces=7
# read length for kraken2-build and bracken
read_len=150
nthreads=1
kraken2_db_dir=data/kraken2_db/k2_large_20240307
#kraken2_db_dir=data/kraken2_db/k2_large_${date_var} 
source ~/miniconda3/etc/profile.d/conda.sh
conda activate kraken2-tools-2.1.3
## Download taxonomy
kraken2-build --download-taxonomy --db ${kraken2_db_dir}
## Download databases
kraken2-build --download-library bacteria --db ${kraken2_db_dir}
kraken2-build --download-library archaea --db ${kraken2_db_dir}
kraken2-build --download-library viral --db ${kraken2_db_dir}
kraken2-build --download-library human --db ${kraken2_db_dir}
kraken2-build --download-library UniVec_Core --db ${kraken2_db_dir}
kraken2-build --download-library fungi --db ${kraken2_db_dir}
kraken2-build --download-library plant --db ${kraken2_db_dir}
kraken2-build --download-library protozoa --db ${kraken2_db_dir}
kraken2-build --download-library plasmid --db ${kraken2_db_dir}
kraken2-build --add-to-library ~/common_data/reference_genomes/Heter_glaber.v1.7_hic_pac_genomic_kraken2.fna  \
  --db ${kraken2_db_dir}
### nt is for fragments of novel organisms (e.g. 16S rRNA gene) that don't have full genomes
#kraken2-build --download-library nt --db ${kraken2_db_dir}
## Build the database: uses taxonomy and library
kraken2-build --build --threads ${nthreads} \
 --db ${kraken2_db_dir} \
 --kmer-len ${kmer_len}  \
 --minimizer-len ${minimizer_len} \
 --minimizer-spaces ${minimizer_spaces}
bracken-build -d ${kraken2_db_dir} \
 -t ${nthreads} \
 -k ${kmer_len} \
 -l ${read_len}
