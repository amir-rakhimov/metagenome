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


# 7. Classify with GTDB-TK
# GTDB-Tk requires ~110G of external data that needs to be downloaded and unarchived
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
conda deactivate
conda activate mag_assembly-tools
# download-db.sh
# wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz \
#  -O ~/miniconda3/envs/mag_assembly-tools/share/gtdbtk-2.4.0/db/gtdbtk_data.tar.gz
# tar -xvf "${GTDBTK_DATA_PATH}"/*tar.gz -C "${GTDBTK_DATA_PATH}"
echo "Running GTDB-TK"
# Fix pplacer locale error
export LC_ALL=C
# for FILE_DIR in "${metabat2_output_dir}"/"${metabat2_date_time}"*bins
# do 
#  SAMPLE=$(echo "${FILE_DIR}" | sed "s/${metabat2_date_time}_//" |sed "s/_bins//")
#  base_name=$(basename "$SAMPLE" )
#  echo "Running GTDB-TK on ${base_name}"
#  gtdbtk classify_wf --genome_dir "${metabat2_output_dir}"/"${metabat2_date_time}"_"${base_name}"_bins/ \
#     --out_dir "${gtdbtk_output_dir}"/"${gtdbtk_date_time}"_"${base_name}"_GTDBtk \
#     -x fa \
#     --skip_ani_screen \
#     --cpus "${nthreads_sort}";
# done
# Run classify separately because pplacer didn't recognize gz files
for FILE_DIR in "${metabat2_output_dir}"/"${metabat2_date_time}"*bins
do 
 SAMPLE=$(echo "${FILE_DIR}" | sed "s/${metabat2_date_time}_//" |sed "s/_bins//")
 base_name=$(basename "$SAMPLE" )
 echo "Running GTDB-TK on ${base_name}"
 gunzip "${gtdbtk_output_dir}"/"${gtdbtk_date_time}"_"${base_name}"_GTDBtk/align/*.gz
 gtdbtk classify --genome_dir "${metabat2_output_dir}"/"${metabat2_date_time}"_"${base_name}"_bins/ \
    --out_dir "${gtdbtk_output_dir}"/"${gtdbtk_date_time}"_"${base_name}"_GTDBtk/classify \
    --align_dir "${gtdbtk_output_dir}"/"${gtdbtk_date_time}"_"${base_name}"_GTDBtk/align \
    -x fa \
    --skip_ani_screen \
    --cpus "${nthreads_sort}";
done
