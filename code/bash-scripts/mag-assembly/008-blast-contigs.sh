#!/usr/bin/env bash
#SBATCH -t 360:00:00
#SBATCH -N 4
#SBATCH -n 20
#SBATCH --mem-per-cpu 8g
#SBATCH -J 20250630_14-34-blast-contigs
#SBATCH --output jobreports/20250630_14-34-blast-contigs-out.txt
#SBATCH --error jobreports/20250630_14-34-blast-contigs-out.txt
#I am requesting 4 nodes containing 20 CPUs, with 8 GB memory per CPU. Total: 160 GB
source ~/miniconda3/etc/profile.d/conda.sh
shopt -s nullglob
# When nullglob is enabled, if a glob pattern does not match any files,
# it expands to nothing (an empty string) instead of returning the pattern itself.
# So, if no matches are found, the script will skip the file

# This script uses BLASTN on SILVA SSU Ref NR99 database on contigs from MEGAHIT
# for a broad taxonomic classification.
date_var=$(date -I|sed 's/-//g')
time_var=$(date +%T |sed 's/:/_/g' )
date_time=${date_var}_${time_var}
start_date_time=$(date +"%F %H:%M:%S")
megahit_date_time=20250303_17_54_26
blastn_date_time=20250701_05_24_10
nthreads=20
nthreads_sort=10
mem_req=8G
mem_req_sort=4G
project_home_dir=~/projects/metagenome
output_dir=~/projects/metagenome/output/mag_assembly
megahit_output_dir=output/mag_assembly/megahit_output
silva_db_dir=data/silva_db
blastdb_path=data/blastdbs/SILVA_138_2_db
blastn_output_dir=output/mag_assembly/blastn_output

cd ${project_home_dir}
mkdir -p output/mag_assembly/blastn_output
mkdir -p data/silva_db
mkdir -p data/blastdbs
# 0. Show the current time for logging
echo "${start_date_time}"

# 1. Get the SILVA SSU Ref NR99 database
# echo "Downloading SILVA SSU Ref NR99 database"
# wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz \
#  -O "${silva_db_dir}"/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz
# gunzip "${silva_db_dir}"/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz
# intermediate_date_time=$(date +"%F %H:%M:%S")
# echo "${intermediate_date_time}"

# 1.1 Build a reference database
# echo "Building a BLAST reference database"
# makeblastdb -in "${silva_db_dir}"/SILVA_138.2_SSURef_NR99_tax_silva.fasta \
#  -parse_seqids \
#  -out "${blastdb_path}" \
#  -logfile "jobreports/"${blastn_date_time}"-SILVA_138_2_SSURef_NR99_tax_silva-makedblastdb-out.txt" \
#  -title "SILVA 138_2 SSURef NR99 db" \
#  -dbtype nucl

# 2. Run BLASTN for a broad taxonomic screening
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
mkdir -p "${blastn_output_dir}"/"${blastn_date_time}"
echo "Running BLASTN"
for FILE_DIR in "${megahit_output_dir}"/"${megahit_date_time}"_*.megahit_asm 
do 
 intermediate_date_time=$(date +"%F %H:%M:%S")
 echo "${intermediate_date_time}"
 SAMPLE=$(echo "${FILE_DIR}" | sed "s/\.megahit_asm//" |sed "s/${megahit_date_time}_//")
 base_name=$(basename "$SAMPLE" )
 echo "Running BLASTN on ${base_name} using SILVA 138_2 SSURef NR99 database"
 blastn -query "${megahit_output_dir}"/"${megahit_date_time}"_"${base_name}".megahit_asm/"${megahit_date_time}"_"${base_name}"_final.contigs.fa \
 -db "${blastdb_path}" \
 -out "${blastn_output_dir}"/"${blastn_date_time}"/"${blastn_date_time}"_blastn_"${base_name}"_SILVA_results.tsv \
 -outfmt "6 qacc sacc qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sframe" \
 -num_threads "${nthreads}" \
 2>&1 |tee  "${blastn_output_dir}"/"${blastn_date_time}"/"${blastn_date_time}"_blastn_"${base_name}"_SILVA.log;
done
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
# qacc means Query accession
# sacc means Subject accession
# qseqid means Query Seq-id
# sseqid means Subject Seq-id
# pident means Percentage of identical matches
# length means Alignment length
# mismatch means Number of mismatches
# gapopen means Number of gap openings
# qstart means Start of alignment in query
# qend means End of alignment in query
# sstart means Start of alignment in subject
# send means End of alignment in subject
# evalue means Expect value
# bitscore means Bit score
# sframe means Subject frame


# 2.1 Extract top BLAST hits with sequence identity > 90% and alignment length > 500 bp
# 1) Filter based on identity and length
# 2) Sort by query ID and descending bitscore
# 3) Pick the top hit per query
# Then, extract qseqid, sseqid, pident, length, bitscore
for FILE in "${blastn_output_dir}"/"${blastn_date_time}"/"${blastn_date_time}"_blastn_*_SILVA_results.tsv
do
 base_name=$(basename "$FILE" )
 base_name=$(echo "${base_name}" | sed "s/_SILVA_results\.tsv//" |sed "s/${blastn_date_time}_blastn_//")
 intermediate_date_time=$(date +"%F %H:%M:%S")
 echo "${intermediate_date_time}"
 # Count all contigs from BLAST
 echo "Counting total number of contigs that had BLAST hits in ${base_name}"
 cut -f 1 "${blastn_output_dir}"/"${blastn_date_time}"/"${blastn_date_time}"_blastn_"${base_name}"_SILVA_results.tsv | \
  uniq  | \
   awk -v sample_id="${base_name}" ' { count++ }
    END {print count " contigs had BLAST hits in sample " sample_id }
    ' 
 # Filter by pident and length (sequence identity > 90% and alignment length > 500 bp), 
 #  and sort by bitscore
 echo "Filtering BLAST hits in ${base_name}"
 awk -F'\t' -v OFS='\t' '$5>90 && $6>500' "${blastn_output_dir}"/"${blastn_date_time}"/"${blastn_date_time}"_blastn_"${base_name}"_SILVA_results.tsv | \
   sort -k1,1 -k14,14nr | awk '!seen[$1]++' > \
    "${blastn_output_dir}"/"${blastn_date_time}"/"${blastn_date_time}"_blastn_"${base_name}"_SILVA_results_filtered.tsv
 #  sort -k1,1 -k14,14nr:
 # sort by qacc (column 1),
 # then by bitscore (column 14), numerically (-n), descending (-r).

 awk -F'\t' -v sample_id="${base_name}" \
    ' { count++ }
    END {print count " contigs had BLAST hits with sequence identity > 90% and alignment length > 500 bp in sample " sample_id }
    ' "${blastn_output_dir}"/"${blastn_date_time}"/"${blastn_date_time}"_blastn_"${base_name}"_SILVA_results_filtered.tsv 
 # Then, extract qseqid, sseqid, pident, length, bitscore
 awk -F'\t' -v OFS='\t' '{print $3, $4, $5, $6, $14}' "${blastn_output_dir}"/"${blastn_date_time}"/"${blastn_date_time}"_blastn_"${base_name}"_SILVA_results_filtered.tsv > \
    "${blastn_output_dir}"/"${blastn_date_time}"/"${blastn_date_time}"_blastn_"${base_name}"_results_for_taxonomy.tsv;
done

# grep "^>" "${silva_db_dir}"/SILVA_138.2_SSURef_NR99_tax_silva.fasta |sed 's/>//' > "${silva_db_dir}"/SILVA_138.2_headers.txt
# Substitute the space between ID and taxonomy with tab
grep "^>" "${silva_db_dir}"/SILVA_138.2_SSURef_NR99_tax_silva.fasta |sed -E 's/^>(\S+)\s+(.*)/\1\t\2/'  > "${silva_db_dir}"/SILVA_138.2_headers.txt

# 2.2 Match BLAST results with SILVA headers
# Final table: 
# contig_id | silva_ref_id | perc_identity | alignment_length | bitscore | taxonomy 
# qseqid, sseqid, pident, length, bitscore
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
for FILE in "${blastn_output_dir}"/"${blastn_date_time}"/"${blastn_date_time}"_blastn_*_results_for_taxonomy.tsv
do
 base_name=$(basename "$FILE" )
 base_name=$(echo "${base_name}" | sed "s/_results_for_taxonomy\.tsv//" |sed "s/${blastn_date_time}_blastn_//")
 intermediate_date_time=$(date +"%F %H:%M:%S")
 echo "${intermediate_date_time}"
 echo "Matching BLASTN results (sseqid) with SILVA database taxonomy on ${base_name}"
 awk -F'\t' -v OFS='\t' -v sample_name="${base_name}" ' 
    NR==FNR { silva[$1] = substr($0, index($0, $2)); next }
    $2 in silva { print sample_name "\t" $0 "\t" silva[$2] }
    ' "${silva_db_dir}"/SILVA_138.2_headers.txt \
        "${blastn_output_dir}"/"${blastn_date_time}"/"${blastn_date_time}"_blastn_"${base_name}"_results_for_taxonomy.tsv > \
        "${blastn_output_dir}"/"${blastn_date_time}"/"${blastn_date_time}"_blastn_"${base_name}"_blast_with_taxonomy.tsv;
done
 intermediate_date_time=$(date +"%F %H:%M:%S")
 echo "${intermediate_date_time}"

# Merge all BLAST taxonomies
cat "${blastn_output_dir}"/"${blastn_date_time}"/"${blastn_date_time}"_blastn_*_blast_with_taxonomy.tsv > \
   "${blastn_output_dir}"/"${blastn_date_time}"/"${blastn_date_time}"_blastn_all_samples.tsv