#!/usr/bin/env bash
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate mag_assembly-tools
set -euo pipefail
shopt -s nullglob
source config/bash/config.sh
SAMPLES_TSV="${1}"

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
# Function that checks if a required file exists or not
require() { [[ -s "$1" ]] || { echo "ERROR: required file missing: $1" >&2; exit 1; }; }
# Function that marks a sample as done (in a hidden directory)
is_done()   { [[ -e "${mark_dir}/${2}.${1}.done" ]]; }   # is_done STEP SAMPLE
mark_done() { touch "${mark_dir}/${2}.${1}.done"; }      # mark_done STEP SAMPLE

# 1. Read the input file list
[[ -s "${SAMPLES_TSV}" ]] || { echo "ERROR: metadata $SAMPLES_TSV missing or empty" >&2; exit 1; }
SAMPLES=()
while IFS=$'\t' read -r sample_id \
    || [[ -n "${sample_id}" ]]
do
    [[ "${sample_id}" == "sample_id" || -z "${sample_id}" ]] && continue
    SAMPLES+=("${sample_id}"); 
done <  "${SAMPLES_TSV}"

(( ${#SAMPLES[@]} > 0 )) || { echo "ERROR: no samples found" >&2; exit 1; } # Counts the number of elements in the final array of samples. If
                                                    # the array is empty, raises an error

# 1. Assemble contigs into bins with MetaBAT2
# Be careful to have bams sorted first!

echo "$(date +"%F %H:%M:%S")"
mark_dir="${metabat2_logs_dir}/.done"
mkdir -p "${mark_dir}"
for SAMPLE_ID in "${SAMPLES[@]}" 
do
    if is_done metabat2 "${SAMPLE_ID}"; then
        echo "Skipping MetaBAT2 for ${SAMPLE_ID}"; continue
    fi
    asm="${megahit_output_dir}/${SAMPLE_ID}.megahit_asm"
    
    # 1.1 Create a depth file
    # The columns represent Contig depth per sample, and variation of that depth for a 
    # particular contig in a particular sample.
    # The rows represent all the contigs in the assembly
    echo "$(date +"%F %H:%M:%S")"
    echo "Running jgi_summarize_bam_contig_depths on ${SAMPLE_ID}"
    jgi_summarize_bam_contig_depths \
        --percentIdentity 97 \
        --minContigLength 1000 \
        --minContigDepth 1.0  \
        --referenceFasta "${asm}/${SAMPLE_ID}_final.contigs.fa" \
        "${megahit_aligned_reads_dir}"/"${SAMPLE_ID}"_either_read_mapped_sorted.bam \
        --outputDepth "${bam_contig_depths_dir}"/"${SAMPLE_ID}".depth.txt \
        2>&1 |tee  "${metabat2_logs_dir}"/"${SAMPLE_ID}"_jgi_summarize_bam_contig_depths.log

    # 1.2 MetaBAT2
    echo "$(date +"%F %H:%M:%S")"
    echo "Running MetaBAT2 on ${SAMPLE_ID}"
    metabat2 \
        -i "${asm}/${SAMPLE_ID}_final.contigs.fa" \
        -a "${bam_contig_depths_dir}"/"${SAMPLE_ID}".depth.txt  \
        -o "${metabat2_output_dir}"/"${SAMPLE_ID}"_bins/"${SAMPLE_ID}"_bin \
        -t "${NTHREADS_METABAT2}" \
        -v \
        2>&1 |tee "${metabat2_logs_dir}"/"${SAMPLE_ID}"_metabat2.log;
    mark_done metabat2 "${SAMPLE_ID}";
done



#  docker run --workdir $(pwd) --volume $(pwd):$(pwd) metabat:latest jgi_summarize_bam_contig_depths \

# runMetaBat.sh <options> assembly.fasta sample1.bam [sample2.bam ...]
# runMetaBat.sh ../megahit_output/2D10_trim_decontam.megahit_asm/2D10_trim_decontam_final.contigs.fa     \
#     ../megahit_output/alignedreads/2D10_aln_either_read_mapped_sorted.bam

echo "$(date +"%F %H:%M:%S")"

#  docker run --workdir $(pwd) --volume $(pwd):$(pwd) metabat:latest metabat2 \

# -o means the output directory and start of the file name (here, bin.1.fa, bin.2.fa, etc)

