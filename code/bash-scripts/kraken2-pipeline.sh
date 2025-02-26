#!/usr/bin/env bash
shopt -s nullglob
# When nullglob is enabled, if a glob pattern does not match any files,
# it expands to nothing (an empty string) instead of returning the pattern itself.
# So, if no matches are found, the script will skip the file

### Get the current date in yyyy-mm-dd format with date -I
### Remove the hyphen with sed (use global search with /g to remove all hyphens)
### Store as a variable date_var
kraken2_db_date=20241209
# time format is 16_28_52 
#time_var=$(date +%T |sed 's/:/_/g' )
#date_var=$(date -I|sed 's/-//g')
#date_time=${date_var}_${time_var}
date_time=20240409_17_32_40
### kmer length for kraken2-build and braken-build
kmer_len=30
minimizer_len=26
minimizer_spaces=6
### 1. minimizer length l must be no more than 31 for nucleotide databases, and 15 for protein databases
### 2. minimizer length l must be no more than the k-mer length
### 3. minimizer_space is s<l/4
### 4. default kmer=35, l=31, s=7
### confidence for kraken2 (default: 0):
confidence_kraken2=0.1
### threshold for bracken (default: 0)
bracken_threshold=10
### read length for kraken2-build and bracken:
read_len=150
nthreads=11
source ~/miniconda3/etc/profile.d/conda.sh
### directory with database and taxonomy
# kraken2_db_dir=data/kraken2_db/k2_large_${kraken2_db_date} 
# mkdir output/kraken2_pipeline/kraken2_reports
# mkdir output/kraken2_pipeline/kraken2_output
# mkdir output/kraken2_pipeline/bracken_reports
# mkdir output/kraken2_pipeline/bracken_output
# mkdir output/kraken2_pipeline/bracken_krona_txt
# mkdir output/kraken2_pipeline/krona_html
# mkdir output/fastqc_output
kraken2_db_dir=data/kraken2_db/k2_large_${kraken2_db_date}
kraken2_reports_dir=output/kraken2_pipeline/kraken2_reports
kraken2_output_dir=output/kraken2_pipeline/kraken2_output
bracken_reports_dir=output/kraken2_pipeline/bracken_reports
bracken_output_dir=output/kraken2_pipeline/bracken_output
bracken_krona_txt_dir=output/kraken2_pipeline/bracken_krona_txt
krona_html_dir=output/kraken2_pipeline/krona_html
### directory with decontaminated fastq files
bowtie2_decontam_fastq_dir=data/bowtie2_decontam_fastq
# Install all tools with conda
## Setup channels for conda installation
# conda config --add channels defaults
# conda config --add channels bioconda
# conda config --add channels conda-forge

## Create a conda environment
# conda create -yqn kraken2-tools-2.1.3 -c conda-forge -c bioconda kraken2 krakentools bracken krona bowtie2

conda activate kraken2-tools-2.1.3

# Build a kraken2 database
# mkdir -p data/metadata
# mkdir -p data/fastq
# mkdir -p data/kraken2_db
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
# kraken2-build --add-to-library ~/common_data/reference_genomes/Heter_glaber.v1.7_hic_pac_genomic_kraken2.fna  \
#   --db ${kraken2_db_dir}
### nt is for fragments of novel organisms (e.g. 16S rRNA gene) that don't have full genomes
# kraken2-build --download-library nt --db ${kraken2_db_dir}

## Build the database: uses taxonomy and library
#kraken2-build --build --threads ${nthreads} \
# --db ${kraken2_db_dir} \
# --kmer-len ${kmer_len}  \
# --minimizer-len ${minimizer_len} \
# --minimizer-spaces ${minimizer_spaces} 2>&1 |tee ${date_time}_kraken2_build_stdout.txt
#bracken-build -d ${kraken2_db_dir} \
# -t ${nthreads} \
# -k ${kmer_len} \
# -l ${read_len} 2>&1 |tee ${date_time}_bracken_build_stdout.txt
### 1. minimizer length must be no more than 31 for nucleotide databases, and 15 for protein databases
### 2. minimizer length must be no more than the k-mer length
### 3. if you change minimizer length, be sure to adjust minimizer spaces:
### A number s < l/4 can be chosen, and s positions in the minimizer will be
#  masked out during all comparisons.
# python ~/miniconda3/envs/kraken2-tools-2.1.3/bin/generate_kmer_distribution.py \
#   -i ${kraken2_db_dir}/database${read_len}mers.kraken \
#   -o ${kraken2_db_dir}/database${read_len}mers.kmer_distrib

# 1. Remove host DNA with bowtie2 (done)

# 2. Classify microbiome samples using Kraken

# We need the database and FASTA file of sequences
# kraken2 --db $DBNAME seqs.fq
# Run with no minimizer data
# If you want to run with minimizers, the output will be 
#  --output ${kraken2_output_dir}/${date_time}_${base_name}.kraken2 instead of 
#  --output ${kraken2_output_dir}/${date_time}_${base_name}_no_minimizer_data.kraken2
# And the --report will be 
#  --report ${kraken2_reports_dir}/${date_time}_${base_name}.k2report
# instead of --report ${kraken2_reports_dir}/${date_time}_${base_name}_no_minimizer_data.k2report
# And you add  --report-minimizer-data 
# to kraken2 command
# And gzip command will be gzip -9 --best ${kraken2_output_dir}/${date_time}_${base_name}.kraken2 
# instead of gzip -9 --best ${kraken2_output_dir}/${date_time}_${base_name}_no_minimizer_data.kraken2
# And tee command will be tee ${date_time}_kraken2_stdout.txt
# instead of tee ${date_time}_kraken2_stdout_no_minimizer_data.txt
for FILE in ${bowtie2_decontam_fastq_dir}/*decontam_R1.fastq.gz
do 
 SAMPLE=$(echo ${FILE} | sed "s/_trim_decontam_R1\.fastq\.gz//")
 base_name=$(basename "$SAMPLE" )
 kraken2 --paired \
	--db ${kraken2_db_dir} \
	--threads ${nthreads} \
	--minimum-hit-groups 3 \
	--gzip-compressed \
	${bowtie2_decontam_fastq_dir}/${base_name}_trim_decontam_R1.fastq.gz \
	${bowtie2_decontam_fastq_dir}/${base_name}_trim_decontam_R2.fastq.gz \
	--output ${kraken2_output_dir}/${date_time}_${base_name}_no_minimizer_data.kraken2 \
	--classified-out ${kraken2_output_dir}/${date_time}_${base_name}_classified_#.fq \
	--report ${kraken2_reports_dir}/${date_time}_${base_name}_no_minimizer_data.k2report \
	--confidence ${confidence_kraken2}
 gzip -9 --best ${kraken2_output_dir}/${date_time}_${base_name}_no_minimizer_data.kraken2
 gzip -9 --best ${kraken2_output_dir}/${date_time}_${base_name}_classified__1.fq
 gzip -9 --best ${kraken2_output_dir}/${date_time}_${base_name}_classified__2.fq;
done  2>&1 |tee ${date_time}_kraken2_stdout_no_minimizer_data.txt   
### --db
### --threads
### -minimum-hit groups
### --gzip-compressed \
### ${bowtie2_decontam_fastq_dir}/${date_time}_${base_name}_decontam_R1.fastq.gz \
### ${bowtie2_decontam_fastq_dir}/${date_time}_${base_name}_decontam_R2.fastq.gz \
### --output ${kraken2_output_dir}/${base_name}.kraken2 \
### --classified-out ${kraken2_output_dir}/${base_name}_classified_#.fq \

### --classified-out will print classified sequences to filename
### --report-minimizer-data adds two columns to the report file (#4 and #5) , which
### represent the number of minimizers found to be associated with a taxon
### in the read sequences (#4), and the estimate of the number of distinct
### minimizers associated with a taxon in the read sequence data (#5)

# 3. Run bracken for abundance estimation of microbiome samples
# Run with no minimizer data
# If you want to run with minimizers, the FOR command will be 
# for FILE in ${kraken2_reports_dir}/${date_time}_*.k2report
# instead of ${kraken2_reports_dir}/${date_time}_*_no_minimizer_data.k2report
# And sed will be SAMPLE=$(echo ${FILE} | sed "s/\.k2report//"|sed "s/${date_time}_//")
# instead of SAMPLE=$(echo ${FILE} | sed "s/\.k2report//"|sed "s/${date_time}_//"|sed "s/_no_minimizer_data//")
# And input will be -i ${kraken2_reports_dir}/${date_time}_${base_name}.k2report 
# instead of -i ${kraken2_reports_dir}/${date_time}_${base_name}_no_minimizer_data.k2report
# And output will be -o ${bracken_output_dir}/${date_time}_${base_name}.bracken 
# instead of -o ${bracken_output_dir}/${date_time}_${base_name}_no_minimizer_data.bracken 
# And the report will be -w ${bracken_reports_dir}/${date_time}_${base_name}.breport
# instead of -w ${bracken_reports_dir}/${date_time}_${base_name}_no_minimizer_data.breport
# And gzip will be gzip -9 --best ${kraken2_reports_dir}/${date_time}_${base_name}.k2report;
# instead of gzip -9 --best ${kraken2_reports_dir}/${date_time}_${base_name}_no_minimizer_data.k2report
# And tee command will be tee ${date_time}_bracken_stdout.txt
# instead of tee ${date_time}_bracken_stdout_no_minimizer_data.txt
for FILE in ${kraken2_reports_dir}/${date_time}_*_no_minimizer_data.k2report
do 
  SAMPLE=$(echo ${FILE} | sed "s/\.k2report//"|sed "s/${date_time}_//"|sed "s/_no_minimizer_data//")
  base_name=$(basename "$SAMPLE" )
  bracken -d ${kraken2_db_dir} \
	-i ${kraken2_reports_dir}/${date_time}_${base_name}_no_minimizer_data.k2report \
	-r ${read_len} \
	-l S \
	-t ${bracken_threshold} \
	-o ${bracken_output_dir}/${date_time}_${base_name}_no_minimizer_data.bracken \
	-w ${bracken_reports_dir}/${date_time}_${base_name}_no_minimizer_data.breport
  gzip -9 --best ${kraken2_reports_dir}/${date_time}_${base_name}_no_minimizer_data.k2report;
done 2>&1 |tee ${date_time}_bracken_stdout_no_minimizer_data.txt
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
for FILE in ${bracken_reports_dir}/${date_time}_*.breport
do 
  SAMPLE=$(echo ${FILE} | sed "s/\.breport//"|sed "s/${date_time}_//")
  base_name=$(basename "$SAMPLE" )
  kreport2krona.py -r ${bracken_reports_dir}/${date_time}_${base_name}.breport \
    -o ${bracken_krona_txt_dir}/${date_time}_${base_name}.b.krona.txt \
	  --no-intermediate-ranks
  ktImportText ${bracken_krona_txt_dir}/${date_time}_${base_name}.b.krona.txt \
	-o ${krona_html_dir}/${date_time}_${base_name}.krona.html;
done
### For kreport2krona.py:
### -r
### -o
### --no-intermediate-ranks
### For ktImportText:
### input is
### -o


# Analyse specific sequences: Macadamia integrifolia 
# 1. Extract reads that were mapped to macadamia
zcat 20240409_17_32_40_H21_trim_classified__1.fq.gz |grep ^@ |grep "kraken:taxid|60698" > 20240409_17_32_40_H21_trim_classified__1_macadamia.txt
conda activate qc-tools
# Use seqtk subseq : doesn't work if fastq header begins with "@"
seqtk subseq 20240409_17_32_40_H21_trim_classified__1.fq.gz 20240409_17_32_40_H21_trim_classified__1_macadamia.txt > H21_r1_macadamia.txt
seqkit grep -f 20240409_17_32_40_H21_trim_classified__1_macadamia.txt 20240409_17_32_40_H21_trim_classified__1.fq.gz -o H21_r1_macadamia.fq.gz

# extract  sequences that contain some pattern in the header 
# -A1 mean print the matching line and 1 line AFTER the match
zcat 20240409_17_32_40_H21_trim_classified__1.fq.gz |grep --no-group-separator  -A1 "kraken:taxid|60698" >  macadamia_r1.fq
zcat 20240409_17_32_40_H21_trim_classified__2.fq.gz |grep --no-group-separator  -A1 "kraken:taxid|60698" > macadamia_r2.fq

sed '1~2s/@/>/g' macadamia_r1.fq > macadamia_r1_blast_input.fq
sed '1~2s/@/>/g' macadamia_r2.fq > macadamia_r2_blast_input.fq

head -n 2 macadamia_r1_blast_input.fq > try.fq
 
blastn -db nt -query try.fq -out results.out -remote
blastn -db nt -query macadamia_r1_blast_input.fq -out macadamia_r1_blast_results.out -remote
blastn -db nt -query macadamia_r2_blast_input.fq -out macadamia_r2_blast_results.out -remote