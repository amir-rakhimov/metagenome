#!/usr/bin/env bash
date_var=$(date -I|sed 's/-//g')
# kmer length for kraken2-build and braken-build
kmer_len=25 
minimizer_len=20
minimizer_spaces=5
# read length for kraken2-build and bracken
read_len=150
nthreads=1
kraken2_db_dir=data/kraken2_db/k2_large_20240111
# kraken2_db_dir=data/kraken2_db/k2_large_${date_var} 
source ~/miniconda3/etc/profile.d/conda.sh
conda activate kraken2-tools-2.1.3
## Download taxonomy
# kraken2-build --download-taxonomy --db ${kraken2_db_dir}
## Download databases
# kraken2-build --download-library bacteria --db ${kraken2_db_dir}
# kraken2-build --download-library archaea --db ${kraken2_db_dir}
# kraken2-build --download-library viral --db ${kraken2_db_dir}
# kraken2-build --download-library human --db ${kraken2_db_dir}
# kraken2-build --download-library UniVec_Core --db ${kraken2_db_dir}
# kraken2-build --download-library fungi --db ${kraken2_db_dir}
# kraken2-build --download-library plant --db ${kraken2_db_dir}
# kraken2-build --download-library protozoa --db ${kraken2_db_dir}
# kraken2-build --download-library plasmid --db ${kraken2_db_dir}
### nt is for fragments of novel organisms (e.g. 16S rRNA gene) that don't have full genomes
# kraken2-build --download-library nt --db ${kraken2_db_dir}
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