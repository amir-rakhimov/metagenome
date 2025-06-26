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

# MetaBAT2 is a binning software: group contigs into bins (genomes)
date_var=$(date -I|sed 's/-//g')
time_var=$(date +%T |sed 's/:/_/g' )
date_time=${date_var}_${time_var}
start_date_time=$(date +"%F %H:%M:%S")
# For binning with MetaBAT2
megahit_date_time=20250303_17_54_26
bbwrap_date_time=20250326_08_18_26
# metabat2_bin_date_time=20250417_23_21_29
metabat2_depth_file_date_time=20250417_23_21_29
metabat2_bin_date_time="${date_time}"
nthreads=40
nthreads_sort=35
mem_req=8G
mem_req_sort=4G
project_home_dir=~/projects/metagenome
output_dir=~/projects/metagenome/output/mag_assembly
megahit_output_dir=output/mag_assembly/megahit_output
megahit_aligned_reads_dir=output/mag_assembly/megahit_output/alignedreads
bam_contig_depths_dir=output/mag_assembly/bam_contig_depths
metabat2_output_dir=output/mag_assembly/metabat2_output
metabat2_reports_dir=output/mag_assembly/metabat2_reports

cd ${project_home_dir}
mkdir -p output/mag_assembly/megahit_output/alignedreads
mkdir -p output/mag_assembly/metabat2_reports
mkdir -p output/mag_assembly/bam_contig_depths
mkdir -p output/mag_assembly/bbwrap_refs

# 0. Show the current time for logging
echo "${start_date_time}"


# 6. Assemble contigs into bins: MetaBAT2, Metaquast, BUSCO, and CheckM2
# Be careful to have bams sorted first!
# 6.1 Activate the environment 


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
#     --outputDepth "${bam_contig_depths_dir}"/"${metabat2_depth_file_date_time}"_"${base_name}".depth.txt \
#     2>&1 |tee  "${metabat2_reports_dir}"/"${metabat2_depth_file_date_time}"_"${base_name}"_jgi_summarize_bam_contig_depths_report.txt
#  gzip -9 --best "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_either_read_mapped_sorted.bam;
# done
#  docker run --workdir $(pwd) --volume $(pwd):$(pwd) metabat:latest jgi_summarize_bam_contig_depths \

# runMetaBat.sh <options> assembly.fasta sample1.bam [sample2.bam ...]
# runMetaBat.sh ../megahit_output/2D10_trim_decontam.megahit_asm/2D10_trim_decontam_final.contigs.fa     \
#     ../megahit_output/alignedreads/2D10_aln_either_read_mapped_sorted.bam

# 6.3 MetaBAT2
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
echo "Running MetaBAT2"
for FILE_DIR in "${megahit_output_dir}"/"${megahit_date_time}"*megahit_asm
# for FILE_DIR in "${megahit_output_dir}"/"${megahit_date_time}"_2D14.megahit_asm \
#  "${megahit_output_dir}"/"${megahit_date_time}"_G14.megahit_asm \
#  "${megahit_output_dir}"/"${megahit_date_time}"_H15.megahit_asm
do 
 SAMPLE=$(echo "${FILE_DIR}" | sed "s/\.megahit_asm//" |sed "s/${megahit_date_time}_//")
 base_name=$(basename "$SAMPLE" )
 echo "Running MetaBAT2 on ${base_name}"
 metabat2 \
    -i "${megahit_output_dir}"/"${megahit_date_time}"_"${base_name}".megahit_asm/"${megahit_date_time}"_"${base_name}"_final.contigs.fa \
    -a "${bam_contig_depths_dir}"/"${metabat2_depth_file_date_time}"_"${base_name}".depth.txt  \
    -o "${metabat2_output_dir}"/"${metabat2_bin_date_time}"_"${base_name}"_bins/"${metabat2_bin_date_time}"_"${base_name}"_bin \
    -v \
    2>&1 |tee "${metabat2_reports_dir}"/"${metabat2_bin_date_time}"_"${base_name}"_metabat2_report.txt;
done
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"

#  docker run --workdir $(pwd) --volume $(pwd):$(pwd) metabat:latest metabat2 \

# metabat2 -i "${megahit_output_dir}"/demo_assembly.fa.gz \
#     -a "${bam_contig_depths_dir}"/demo.depth.txt  \
#     -o "${metabat2_output_dir}"/"${date_time}"_demo/"${date_time}"_demo_bin \
#     -v
# -o means the output directory and start of the file name (here, bin.1.fa, bin.2.fa, etc)

