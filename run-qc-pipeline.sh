#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob
qc_RUN_ID=$(date +%Y%m%d_%H%M%S)
# Load config
source config/bash/config.sh

## data
mkdir -p "${fastq_dir}"
mkdir -p $HOME/common_data
mkdir -p $HOME/common_data/reference_genomes

# QC and decontamination
mkdir -p "${fastqc_raw_output_dir}"
mkdir -p "${multiqc_raw_output_dir}"
mkdir -p "${fastqc_trimmed_output_dir}"
mkdir -p "${multiqc_trimmed_output_dir}"
mkdir -p "${cutadapt_output_dir}"

mkdir -p "${bowtie2_sam_dir}"
mkdir -p "${bowtie2_bam_dir}"
mkdir -p "${bowtie2_filtered_bam_dir}"
mkdir -p "${bowtie2_sorted_bam_dir}"
mkdir -p "${fastqc_decontam_dir}"
mkdir -p "${multiqc_decontam_dir}"
mkdir -p "${bowtie2_decontam_fastq_dir}"

job_id1=$(sbatch --parsable --export=ALL code/slurm-scripts/qc-pipeline/qc-commands-with-bowtie2.slurm)
echo "Pipeline submitted"
echo " step 1: ${job_id1}"
