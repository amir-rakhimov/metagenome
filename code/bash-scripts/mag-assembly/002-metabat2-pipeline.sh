#!/usr/bin/env bash
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate mag_assembly-tools
set -euo pipefail
shopt -s nullglob
source config/bash/config.sh
# When nullglob is enabled, if a glob pattern does not match any files,
# it expands to nothing (an empty string) instead of returning the pattern itself.
# So, if no matches are found, the script will skip the file

# This script uses MetaBAT2 to bin contigs (group contigs into bins (genomes)).
# MetaBAT2 is run locally, not in conda. Be careful to have bams sorted first!
# Before running metabat2, we need to create depth files with jgi_summarize_bam_contig_depths.
# MetaBAT2 creates a separate directory for each sample which contains all the bins.
# MetaBAT2 is using random seeds, but it's not necessary to set the random seed (results are usually similar).
# Final output: directories with bins (multi-FASTA files).
echo "$(date +"%F %H:%M:%S")"

# 1. Assemble contigs into bins with MetaBAT2
# Be careful to have bams sorted first!

# 1.1 Create a depth file
# The columns represent Contig depth per sample, and variation of that depth for a 
# particular contig in a particular sample.
# The rows represent all the contigs in the assembly
echo "$(date +"%F %H:%M:%S")"
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
echo "$(date +"%F %H:%M:%S")"
echo "Running MetaBAT2"
for FILE_DIR in "${megahit_output_dir}"/*megahit_asm
# for FILE_DIR in "${megahit_output_dir}"/2D14.megahit_asm \
#  "${megahit_output_dir}"/G14.megahit_asm \
#  "${megahit_output_dir}"/H15.megahit_asm
do 
    SAMPLE=$(echo "${FILE_DIR}" | sed "s/\.megahit_asm//")
    base_name=$(basename "$SAMPLE" )
    echo "$(date +"%F %H:%M:%S")"
    echo "Running MetaBAT2 on ${base_name}"
    metabat2 \
        -i "${megahit_output_dir}"/"${base_name}".megahit_asm/"${base_name}"_final.contigs.fa \
        -a "${bam_contig_depths_dir}"/"${base_name}".depth.txt  \
        -o "${metabat2_output_dir}"/"${base_name}"_bins/"${base_name}"_bin \
        -v \
        2>&1 |tee "${metabat2_logs_dir}"/"${base_name}"_metabat2.log;
done
echo "$(date +"%F %H:%M:%S")"

#  docker run --workdir $(pwd) --volume $(pwd):$(pwd) metabat:latest metabat2 \

# -o means the output directory and start of the file name (here, bin.1.fa, bin.2.fa, etc)

