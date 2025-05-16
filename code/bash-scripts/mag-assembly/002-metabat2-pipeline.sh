#!/usr/bin/env bash
#SBATCH -t 360:00:00
#SBATCH -N 4
#SBATCH -n 40
#SBATCH --mem-per-cpu 8g
#SBATCH -J 20250515_11-06-mag_assembly
#SBATCH --output jobreports/20250515_11-06-mag_assembly-megahit-metabat2-pipeline-out.txt
#SBATCH --error jobreports/20250515_11-06-mag_assembly-megahit-metabat2-pipeline-out.txt
#I am requesting 4 nodes containing 40 CPUs, with 8 GB memory per CPU. Total: 320 GB. 280 GB for GTDB-TK
source ~/miniconda3/etc/profile.d/conda.sh
shopt -s nullglob
# When nullglob is enabled, if a glob pattern does not match any files,
# it expands to nothing (an empty string) instead of returning the pattern itself.
# So, if no matches are found, the script will skip the file

# Metaquast is a metagenome assembly evaluation tool
# MetaBAT2 is a binning software: group contigs into bins (genomes)
# Metaquast is also used to evaluate the quality of bins
# CheckM2 assesses the quality of bins
# GTDB-TK assigns taxonomy to genomes
# Prodigal is a gene prediction software
# Prokka annotates contigs to find bacterial genes
# Metaeuk annotates contigs to find eukaryotic genes
date_var=$(date -I|sed 's/-//g')
time_var=$(date +%T |sed 's/:/_/g' )
date_time=${date_var}_${time_var}
start_date_time=$(date +"%F %H:%M:%S")
# For binning with MetaBAT2
megahit_date_time=20250303_17_54_26
metaquast_megahit_date_time=20250303_17_54_26
prodigal_date_time=20250303_17_54_26
bbwrap_date_time=20250326_08_18_26
metabat2_date_time=20250417_23_21_29
metaquast_metabat2_date_time=20250417_23_21_29
busco_date_time=20250417_23_21_29
checkm2_date_time=20250501_07_53_55
gtdbtk_date_time=20250514_17_57_12
prokka_date_time="${date_time}"
nthreads=40
nthreads_sort=35
mem_req=8G
mem_req_sort=4G
project_home_dir=~/projects/metagenome
output_dir=~/projects/metagenome/output/mag_assembly
bowtie2_decontam_fastq_dir=data/bowtie2_decontam_fastq
megahit_output_dir=output/mag_assembly/megahit_output
megahit_aligned_reads_dir=output/mag_assembly/megahit_output/alignedreads
contig_coverages_dir=output/mag_assembly/contig_coverages
prodigal_output_dir=output/mag_assembly/prodigal_output
metaquast_output_dir=output/mag_assembly/metaquast_output
metaquast_script_dir=~/quast-5.2.0
bam_contig_depths_dir=output/mag_assembly/bam_contig_depths
metabat2_output_dir=output/mag_assembly/metabat2_output
metabat2_reports_dir=output/mag_assembly/metabat2_reports
samtools_reports_dir=output/mag_assembly/samtools_reports
checkm2_output_dir=output/mag_assembly/checkm2_output
gtdbtk_output_dir=output/mag_assembly/gtdbtk_output
bbwrap_refs_dir=output/mag_assembly/bbwrap_refs
busco_output_dir=output/mag_assembly/busco_output
prokka_output_dir=output/mag_assembly/prokka_output

cd ${project_home_dir}
mkdir -p output/mag_assembly/megahit_output/alignedreads
mkdir -p output/mag_assembly/metabat2_reports
mkdir -p output/mag_assembly/samtools_reports
mkdir -p output/mag_assembly/checkm2_output
mkdir -p output/mag_assembly/gtdbtk_output
mkdir -p output/mag_assembly/bam_contig_depths
mkdir -p output/mag_assembly/bbwrap_refs
mkdir -p output/mag_assembly/metaquast_output
mkdir -p output/mag_assembly/busco_output
mkdir -p output/mag_assembly/prodigal_output
mkdir -p output/mag_assembly/prokka_output
mkdir -p output/mag_assembly/contig_coverages

# 0. Show the current time for logging
echo "${start_date_time}"

# 1. Activate conda environment
conda activate mag_assembly-tools



# 6. Assemble contigs into bins: MetaBAT2, Metaquast, BUSCO, and CheckM2
# Be careful to have bams sorted first!
# 6.1 Activate the environment again
conda deactivate
conda activate mag_assembly-tools

# 6.2 Create a depth file
# The columns represent Contig depth per sample, and variation of that depth for a 
# particular contig in a particular sample.
# The rows represent all the contigs in the assembly
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
# for FILE_DIR in "${megahit_output_dir}"/"${megahit_date_time}"_*.megahit_asm 
# # for FILE_DIR in "${megahit_output_dir}"/"${megahit_date_time}"_2D14.megahit_asm \
# #  "${megahit_output_dir}"/"${megahit_date_time}"_G14.megahit_asm \
# #  "${megahit_output_dir}"/"${megahit_date_time}"_H15.megahit_asm
# do 
#  SAMPLE=$(echo "${FILE_DIR}" | sed "s/\.megahit_asm//" |sed "s/${megahit_date_time}_//")
#  base_name=$(basename "$SAMPLE" )
#  echo "Running jgi_summarize_bam_contig_depths on ${base_name}"
#  jgi_summarize_bam_contig_depths \
#     --percentIdentity 97 \
#     --minContigLength 1000 \
#     --minContigDepth 1.0  \
#     --referenceFasta "${megahit_output_dir}"/"${megahit_date_time}"_"${base_name}".megahit_asm/"${megahit_date_time}"_"${base_name}"_final.contigs.fa \
#     "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_either_read_mapped_sorted.bam \
#     --outputDepth "${bam_contig_depths_dir}"/"${metabat2_date_time}"_"${base_name}".depth.txt \
#     2>&1 |tee  "${metabat2_reports_dir}"/"${metabat2_date_time}"_"${base_name}"_jgi_summarize_bam_contig_depths_report.txt
#  gzip -9 --best "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_either_read_mapped_sorted.bam;
# done
#  docker run --workdir $(pwd) --volume $(pwd):$(pwd) metabat:latest jgi_summarize_bam_contig_depths \

# runMetaBat.sh <options> assembly.fasta sample1.bam [sample2.bam ...]
# runMetaBat.sh ../megahit_output/2D10_trim_decontam.megahit_asm/2D10_trim_decontam_final.contigs.fa     \
#     ../megahit_output/alignedreads/2D10_aln_either_read_mapped_sorted.bam

# 6.3 MetaBAT2
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
# echo "Running MetaBAT2"
# for FILE_DIR in "${megahit_output_dir}"/"${megahit_date_time}"*megahit_asm
# # for FILE_DIR in "${megahit_output_dir}"/"${megahit_date_time}"_2D14.megahit_asm \
# #  "${megahit_output_dir}"/"${megahit_date_time}"_G14.megahit_asm \
# #  "${megahit_output_dir}"/"${megahit_date_time}"_H15.megahit_asm
# do 
#  SAMPLE=$(echo "${FILE_DIR}" | sed "s/\.megahit_asm//" |sed "s/${megahit_date_time}_//")
#  base_name=$(basename "$SAMPLE" )
#  echo "Running MetaBAT2 on ${base_name}"
#  metabat2 \
#     -i "${megahit_output_dir}"/"${megahit_date_time}"_"${base_name}".megahit_asm/"${megahit_date_time}"_"${base_name}"_final.contigs.fa \
#     -a "${bam_contig_depths_dir}"/"${metabat2_date_time}"_"${base_name}".depth.txt  \
#     -o "${metabat2_output_dir}"/"${metabat2_date_time}"_"${base_name}"_bins/"${metabat2_date_time}"_"${base_name}"_bin \
#     -v \
#     2>&1 |tee "${metabat2_reports_dir}"/"${metabat2_date_time}"_"${base_name}"_metabat2_report.txt;
# done

#  docker run --workdir $(pwd) --volume $(pwd):$(pwd) metabat:latest metabat2 \

# metabat2 -i "${megahit_output_dir}"/demo_assembly.fa.gz \
#     -a "${bam_contig_depths_dir}"/demo.depth.txt  \
#     -o "${metabat2_output_dir}"/"${date_time}"_demo/"${date_time}"_demo_bin \
#     -v
# -o means the output directory and start of the file name (here, bin.1.fa, bin.2.fa, etc)

# 6.4 QC MetaBAT2 bins with QUAST
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
# for FILE_DIR in "${metabat2_output_dir}"/"${metabat2_date_time}"*bins
# # for FILE_DIR in "${metabat2_output_dir}"/"${metabat2_date_time}"_2D14_bins \
# #  "${metabat2_output_dir}"/"${metabat2_date_time}"_G14_bins \
# #  "${metabat2_output_dir}"/"${metabat2_date_time}"_H15_bins 
# do 
#  SAMPLE=$(echo "${FILE_DIR}" | sed "s/${metabat2_date_time}_//" | sed "s/_bins//")
#  base_name=$(basename "$SAMPLE" )
#  echo "Running MetaQuast on ${base_name} bins" 
#  "${metaquast_script_dir}"/metaquast.py \
#     "${metabat2_output_dir}"/"${metabat2_date_time}"_"${base_name}"_bins/"${metabat2_date_time}"_"${base_name}"_bin* \
#     -o "${metaquast_output_dir}"/"${metaquast_metabat2_date_time}"_metaquast_metabat2_"${base_name}" \
#     --threads "${nthreads}" \
#     2>&1 |tee "${metaquast_output_dir}"/"${metaquast_metabat2_date_time}"_"${base_name}"_metaquast_metabat2_report.txt;
# done

# 6.6 CheckM2
conda deactivate
conda activate checkm2-env
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
# echo "Running CheckM2"
# for FILE_DIR in "${metabat2_output_dir}"/"${metabat2_date_time}"*bins
# do 
#  SAMPLE=$(echo "${FILE_DIR}" | sed "s/${metabat2_date_time}_//"| sed "s/_bins//")
#  base_name=$(basename "$SAMPLE" )
#  echo "Running CheckM2 on ${base_name}" 
#  checkm2 predict --threads "${nthreads_sort}" \
#    -x fa \
#    --database_path "${project_home_dir}"/data/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd \
#    --input "${metabat2_output_dir}"/"${metabat2_date_time}"_"${base_name}"_bins \
#    --output-directory "${checkm2_output_dir}"/"${checkm2_date_time}"_"${base_name}"_checkm2 \
#    2>&1 |tee "${checkm2_output_dir}"/"${checkm2_date_time}"_"${base_name}"_checkm2_report.txt;
# done
# Usage: checkm2 predict --threads 30 --input <folder_with_bins> --output-directory <output_folder> 


# drep to get unique MAGs