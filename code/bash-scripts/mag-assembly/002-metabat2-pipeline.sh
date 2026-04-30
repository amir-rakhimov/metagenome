#!/usr/bin/env bash
#SBATCH -t 360:00:00
#SBATCH -N 4
#SBATCH -n 40
#SBATCH --mem-per-cpu 8g
#SBATCH -J 20250612_13-36-metabat2
#SBATCH --output jobreports/20250612_13-36-metabat2-out.txt
#SBATCH --error jobreports/20250612_13-36-metabat2-out.txt
#I am requesting 4 nodes containing 40 CPUs, with 8 GB memory per CPU. Total: 320 GB
source ~/miniconda3/etc/profile.d/conda.sh
shopt -s nullglob
# When nullglob is enabled, if a glob pattern does not match any files,
# it expands to nothing (an empty string) instead of returning the pattern itself.
# So, if no matches are found, the script will skip the file

# This script uses MetaBAT2 to bin contigs (group contigs into bins (genomes)).
# MetaBAT2 is run locally, not in conda. Be careful to have bams sorted first!
# Before running metabat2, we need to create depth files with jgi_summarize_bam_contig_depths.
# MetaBAT2 creates a separate directory for each sample which contains all the bins.
# MetaBAT2 is using random seeds, but it's not necessary to set the random seed (results are usually similar).
# Final output: directories with bins (multi-FASTA files).
date_var=$(date -I|sed 's/-//g')
time_var=$(date +%T |sed 's/:/_/g' )
date_time=${date_var}_${time_var}
start_date_time=$(date +"%F %H:%M:%S")


# 0. Show the current time for logging
echo "${start_date_time}"


# 1. Assemble contigs into bins with MetaBAT2
# Be careful to have bams sorted first!

# 1.1 Create a depth file
# The columns represent Contig depth per sample, and variation of that depth for a 
# particular contig in a particular sample.
# The rows represent all the contigs in the assembly
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
for FILE_DIR in "${megahit_output_dir}"/*.megahit_asm 
# # for FILE_DIR in "${megahit_output_dir}"/2D14.megahit_asm \
# #  "${megahit_output_dir}"/G14.megahit_asm \
# #  "${megahit_output_dir}"/H15.megahit_asm
do 
 SAMPLE=$(echo "${FILE_DIR}" | sed "s/\.megahit_asm//")
 base_name=$(basename "$SAMPLE" )
 echo "Running jgi_summarize_bam_contig_depths on ${base_name}"
 jgi_summarize_bam_contig_depths \
    --percentIdentity 97 \
    --minContigLength 1000 \
    --minContigDepth 1.0  \
    --referenceFasta "${megahit_output_dir}"/"${base_name}".megahit_asm/"${base_name}"_final.contigs.fa \
    "${megahit_aligned_reads_dir}"/"${base_name}"_either_read_mapped_sorted.bam \
    --outputDepth "${bam_contig_depths_dir}"/"${base_name}".depth.txt \
    2>&1 |tee  "${metabat2_logs_dir}"/"${base_name}"_jgi_summarize_bam_contig_depths.log
 gzip -9 --best "${megahit_aligned_reads_dir}"/"${base_name}"_either_read_mapped_sorted.bam;
done
#  docker run --workdir $(pwd) --volume $(pwd):$(pwd) metabat:latest jgi_summarize_bam_contig_depths \

# runMetaBat.sh <options> assembly.fasta sample1.bam [sample2.bam ...]
# runMetaBat.sh ../megahit_output/2D10_trim_decontam.megahit_asm/2D10_trim_decontam_final.contigs.fa     \
#     ../megahit_output/alignedreads/2D10_aln_either_read_mapped_sorted.bam

# 1.2 MetaBAT2
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
echo "Running MetaBAT2"
for FILE_DIR in "${megahit_output_dir}"/*megahit_asm
# for FILE_DIR in "${megahit_output_dir}"/2D14.megahit_asm \
#  "${megahit_output_dir}"/G14.megahit_asm \
#  "${megahit_output_dir}"/H15.megahit_asm
do 
 SAMPLE=$(echo "${FILE_DIR}" | sed "s/\.megahit_asm//")
 base_name=$(basename "$SAMPLE" )
 echo "Running MetaBAT2 on ${base_name}"
 metabat2 \
    -i "${megahit_output_dir}"/"${base_name}".megahit_asm/"${base_name}"_final.contigs.fa \
    -a "${bam_contig_depths_dir}"/"${base_name}".depth.txt  \
    -o "${metabat2_output_dir}"/"${base_name}"_bins/"${base_name}"_bin \
    -v \
    2>&1 |tee "${metabat2_logs_dir}"/"${base_name}"_metabat2.log;
done
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"

#  docker run --workdir $(pwd) --volume $(pwd):$(pwd) metabat:latest metabat2 \

# -o means the output directory and start of the file name (here, bin.1.fa, bin.2.fa, etc)

