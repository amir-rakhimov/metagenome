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

# MEGAHIT assembles reads into contigs
# Metaquast is a metagenome assembly evaluation tool
# MetaBAT2 is a binning software: group contigs into bins (genomes)
# Metaquast is also used to evaluate the quality of bins
date_var=$(date -I|sed 's/-//g')
time_var=$(date +%T |sed 's/:/_/g' )
date_time=${date_var}_${time_var}
start_date_time=$(date +"%F %H:%M:%S")
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


# 4. Prodigal
# Anonymous mode (for metagenomes)
# should I use -p meta or -p anon?
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
# for FILE_DIR in ${megahit_output_dir}/"${megahit_date_time}"*megahit_asm
# do 
#  SAMPLE=$(echo "${FILE_DIR}" | sed "s/\.megahit_asm//" |sed "s/${megahit_date_time}_//")
#  base_name=$(basename "$SAMPLE" )
#  mkdir "${prodigal_output_dir}"/"${prodigal_date_time}"_"${base_name}"_prodigal
#  prodigal \
#     -i "${megahit_output_dir}"/"${megahit_date_time}"_"${base_name}".megahit_asm/"${megahit_date_time}"_"${base_name}"_final.contigs.fa \
#     -o "${prodigal_output_dir}"/"${prodigal_date_time}"_"${base_name}"_prodigal/"${prodigal_date_time}"_"${base_name}"_prodigal_coords.gbk \
#     -a "${prodigal_output_dir}"/"${prodigal_date_time}"_"${base_name}"_prodigal/"${prodigal_date_time}"_"${base_name}"_prodigal_proteins.faa \
#     -p meta 2>&1 |tee "${prodigal_output_dir}"/"${prodigal_date_time}"_"${base_name}"_prodigal_report.txt;
# done


# Prokka
# Input is contigs
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
echo "Running Prokka"
for FILE_DIR in ${megahit_output_dir}/"${megahit_date_time}"*megahit_asm
do 
 SAMPLE=$(echo "${FILE_DIR}" | sed "s/\.megahit_asm//" |sed "s/${megahit_date_time}_//")
 base_name=$(basename "$SAMPLE" )
 prokka \
    --outdir "${prokka_output_dir}"/"${prokka_date_time}"_"${base_name}"_prokka \
    --prefix "${prokka_date_time}"_"${base_name}"_pred_prokka \
    "${megahit_output_dir}"/"${megahit_date_time}"_"${base_name}".megahit_asm/"${megahit_date_time}"_"${base_name}"_final.contigs.fa \
    --cpus "${nthreads_sort}" 2>&1 |tee "${prokka_output_dir}"/"${prokka_date_time}"_"${base_name}"_prokka_report.txt;
done
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"

# Metaeuk for eukaryotes
# Input is contigs
