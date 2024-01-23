#!/usr/bin/env bash
### Get the current date in yyyy-mm-dd format with date -I
### Remove the hyphen with sed (use global search with /g to remove all hyphens)
### Store as a variable date_var
# date_var=$(date -I|sed 's/-//g')
date_var=20240122
# kmer length for kraken2-build and braken-build
kmer_len=35 
# read length for kraken2-build and bracken
read_len=150
nthreads=10
source ~/miniconda3/etc/profile.d/conda.sh
# directory with database and taxonomy
# kraken2_db_dir=data/kraken2_db/k2_large_${date_var} 
kraken2_db_dir=data/kraken2_db/k2_large_20240111
kraken2_reports_dir=output/kraken2_reports
kraken2_output_dir=output/kraken2_output
bracken_reports_dir=output/bracken_reports
bracken_output_dir=output/bracken_output
bracken_krona_txt_dir=output/bracken_krona_txt
krona_html_dir=output/krona_html
# directory with decontaminated fastq files
bowtie2_host_removed_fastq_dir=data/bowtie2_host_removed_fastq
# Install all tools with conda
## Setup channels for conda installation
# conda config --add channels defaults
# conda config --add channels bioconda
# conda config --add channels conda-forge

## Create a conda environment
# conda create -yqn kraken2-tools-2.1.3 -c conda-forge -c bioconda kraken2 krakentools bracken krona bowtie2

conda activate kraken2-tools-2.1.3

# Build a kraken2 database
# mkdir data
# mkdir data/metadata
# mkdir data/fastq
# mkdir data/kraken2_db
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
# kraken2-build --build --threads ${nthreads} \
#  --db ${kraken2_db_dir} \
#  --kmer-len ${kmer_len}  \
#  --minimizer-len 31
# bracken-build -d ${kraken2_db_dir} -t ${nthreads} -k ${kmer_len} -l ${read_len}

# 1. Remove host DNA with bowtie2 (done)

# 2. Classify microbiome samples using Kraken
# mkdir output/kraken2_reports
# mkdir output/kraken2_output
# mkdir output/bracken_reports
# mkdir output/bracken_output
# mkdir output/bracken_krona_txt
# mkdir output/krona_html
# mkdir output/fastqc_output
# We need the database and FASTA file of sequences
# kraken2 --db $DBNAME seqs.fq

# for FILE in ${bowtie2_host_removed_fastq_dir}/*host_removed_R1.fastq.gz
# do 
#  SAMPLE=$(echo ${FILE} | sed "s/_host_removed_R1\.fastq\.gz//")
#  base_name=$(basename "$SAMPLE" )
#  kraken2 --paired \
# 	--db ${kraken2_db_dir} \
# 	--threads ${nthreads} \
# 	--minimum-hit-groups 3 \
# 	--gzip-compressed \
# 	${bowtie2_host_removed_fastq_dir}/${base_name}_host_removed_R1.fastq.gz \
# 	${bowtie2_host_removed_fastq_dir}/${base_name}_host_removed_R2.fastq.gz \
# 	--output ${kraken2_output_dir}/${date_var}_${base_name}.kraken2 \
# 	--classified-out ${kraken2_output_dir}/${date_var}_${base_name}_classified_#.fq \
# 	--report ${kraken2_reports_dir}/${date_var}_${base_name}.k2report \
# 	--report-minimizer-data; #\
# 	#--confidence ???;
# done
### --db
### --threads
### -minimum-hit groups
### --gzip-compressed \
### ${bowtie2_host_removed_fastq_dir}/${date_var}_${base_name}_host_removed_R1.fastq.gz \
### ${bowtie2_host_removed_fastq_dir}/${date_var}_${base_name}_host_removed_R2.fastq.gz \
### --output ${kraken2_output_dir}/${base_name}.kraken2 \
### --classified-out ${kraken2_output_dir}/${base_name}_classified_#.fq \
	
### --classified-out will print classified sequences to filename
### --report-minimizer-data adds two columns to the report file (#4 and #5) , which
### represent the number of minimizers found to be associated with a taxon 
### in the read sequences (#4), and the estimate of the number of distinct 
### minimizers associated with a taxon in the read sequence data (#5)

# 3. Run bracken for abundance estimation of microbiome samples
for FILE in ${kraken2_reports_dir}/${date_var}_*.k2report
do 
  SAMPLE=$(echo ${FILE} | sed "s/\.k2report//"|sed "s/${date_var}_//")
  base_name=$(basename "$SAMPLE" )
  bracken -d ${kraken2_db_dir} \
	-i ${kraken2_reports_dir}/${date_var}_${base_name}.k2report \
	-r ${read_len} \
	-l S \
	-t 10 \
	-o ${bracken_output_dir}/${date_var}_${base_name}.bracken \
	-w ${bracken_reports_dir}/${date_var}_${base_name}.breport;
done
### -d
### -i
### -r is read length
### -l is classification level
### -t is threshold: specifies the minimum number of reads required for 
### a classification at the specified rank. Any classifications with less than 
### the specified threshold will not receive additional reads from higher taxonomy 
### levels when distributing reads for abundance estimation.
### -o
### -w

#FILE: testdir/20230104_f1.txt
#SAMPLE: testdir/f1
#base_name: f1

#FILE:testdir/20230104_f2.txt
#SAMPLE:testdir/f2
#base_name:f2

#FILE: testdir/20230104_file.txt
#SAMPLE:testdir/file
#base_name:file

# 4. Generate Krona plots
for FILE in ${bracken_reports_dir}/${date_var}_*.breport
do 
  SAMPLE=$(echo ${FILE} | sed "s/\.breport//"|sed "s/${date_var}_//")
  base_name=$(basename "$SAMPLE" )
  kreport2krona.py -r ${bracken_reports_dir}/${date_var}_${base_name}.breport \
    -o ${bracken_krona_txt_dir}/${date_var}_${base_name}.b.krona.txt \
	  --no-intermediate-ranks
  ktImportText ${bracken_krona_txt_dir}/${date_var}_${base_name}.b.krona.txt \
	-o ${krona_html_dir}/${date_var}_${base_name}.krona.html;
done
### For kreport2krona.py:
### -r
### -o
### --no-intermediate-ranks
### For ktImportText:
### input is
### -o
