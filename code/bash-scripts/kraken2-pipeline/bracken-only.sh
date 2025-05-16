#!/usr/bin/env bash
### Get the current date in yyyy-mm-dd format with date -I
### Remove the hyphen with sed (use global search with /g to remove all hyphens)
### Store as a variable date_var
kraken2_db_date=20240307
date_var=20240401
# time format is 16_28_52 
time_var=$(date +%T |sed 's/:/_/g' )
bracken_date=$(date -I|sed 's/-//g')
bracken_date_time=${bracken_date}_${time_var}
# date_var=20240307
### kmer length for kraken2-build and braken-build
kmer_len=25
minimizer_len=20
minimizer_spaces=5
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
nthreads=10
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

conda activate kraken2-tools-2.1.3


# 3. Run bracken for abundance estimation of microbiome samples
# for FILE in ${kraken2_reports_dir}/${date_time}_*.k2report
# do 
#   SAMPLE=$(echo ${FILE} | sed "s/\.k2report//"|sed "s/${date_time}_//")
#   base_name=$(basename "$SAMPLE" )
#   bracken -d ${kraken2_db_dir} \
# 	-i ${kraken2_reports_dir}/${date_time}_${base_name}.k2report \
# 	-r ${read_len} \
# 	-l S \
# 	-t ${bracken_threshold} \
# 	-o ${bracken_output_dir}/${date_time}_${base_name}.bracken \
# 	-w ${bracken_reports_dir}/${date_time}_${base_name}.breport;
# done
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
# Run with no minimizer data
for FILE in ${kraken2_reports_dir}/${date_time}_*_no_minimizer_data.k2report
do 
  SAMPLE=$(echo ${FILE} | sed "s/\.k2report//"|sed "s/${date_time}_//"|sed "s/_no_minimizer_data//")
  base_name=$(basename "$SAMPLE" )
  bracken -d ${kraken2_db_dir} \
	-i ${kraken2_reports_dir}/${date_time}_${base_name}_no_minimizer_data.k2report \
	-r ${read_len} \
	-l S \
	-t ${bracken_threshold} \
	-o ${bracken_output_dir}/${bracken_date_time}_${base_name}_no_minimizer_data.bracken \
	-w ${bracken_reports_dir}/${bracken_date_time}_${base_name}_no_minimizer_data.breport;
done 2>&1 |tee ${bracken_date_time}_bracken_stdout_no_minimizer_data.txt
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
for FILE in ${bracken_reports_dir}/${bracken_date_time}_*.breport
do 
  SAMPLE=$(echo ${FILE} | sed "s/\.breport//"|sed "s/${bracken_date_time}_//")
  base_name=$(basename "$SAMPLE" )
  kreport2krona.py -r ${bracken_reports_dir}/${bracken_date_time}_${base_name}.breport \
    -o ${bracken_krona_txt_dir}/${bracken_date_time}_${base_name}.b.krona.txt \
	  --no-intermediate-ranks
  ktImportText ${bracken_krona_txt_dir}/${bracken_date_time}_${base_name}.b.krona.txt \
	-o ${krona_html_dir}/${bracken_date_time}_${base_name}.krona.html;
done
### For kreport2krona.py:
### -r
### -o
### --no-intermediate-ranks
### For ktImportText:
### input is
### -o
# for FILE in ${kraken2_reports_dir}/${date_time}_*.k2report
# do 
#   SAMPLE=$(echo ${FILE} | sed "s/\.k2report//"|sed "s/${date_time}_//")
#   base_name=$(basename "$SAMPLE" )
#   cat ${kraken2_reports_dir}/${date_time}_${base_name}.k2report | \
#   	awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$6"\t"$7"\t"$8}'\
#   	>${kraken2_reports_dir}/${date_time}_${base_name}_no_minimizer_data.k2report;
# done

# for FILE in ${kraken2_reports_dir}/${date_time}_*.k2report
# do 
#   SAMPLE=$(echo ${FILE} | sed "s/\.k2report//"|sed "s/${date_time}_//")
#   base_name=$(basename "$SAMPLE" )
#   cut -f1-3,6-8 ${kraken2_reports_dir}/${date_time}_${base_name}.k2report > \
#   	${kraken2_reports_dir}/${date_time}_${base_name}_no_minimizer_data.k2report;
# done


# for FILE in ${kraken2_reports_dir}/${date_time}_*_no_minimizer_data.k2report
# do 
#   SAMPLE=$(echo ${FILE} | sed "s/\.k2report//"|sed "s/${date_time}_//"|sed ""s/_no_minimizer_data//"")
#   base_name=$(basename "$SAMPLE" )
#   kreport2mpa.py -r ${kraken2_reports_dir}/${date_time}_${base_name}_no_minimizer_data.k2report \
#     -o ${kraken2_reports_dir}/${date_time}_${base_name}_no_minimizer_data.mpa.tsv \
# 	--display-header;
# done
# for FILE in ${bracken_reports_dir}/${date_time}_*_no_minimizer_data.breport
# do 
#   SAMPLE=$(echo ${FILE} | sed "s/\.breport//"|sed "s/${date_time}_//"|sed "s/_no_minimizer_data//")
#   base_name=$(basename "$SAMPLE" )
#   kreport2mpa.py -r ${bracken_reports_dir}/${date_time}_${base_name}_no_minimizer_data.breport \
#     -o ${bracken_reports_dir}/${date_time}_${base_name}_no_minimizer_data.breport.mpa.tsv \
# 	--display-header;
# done


# combine_mpa.py -i ${bracken_reports_dir}/${date_time}_*_no_minimizer_data.breport.mpa.tsv \
# 	-o ${bracken_reports_dir}/${date_time}_combined_mpa.tsv

