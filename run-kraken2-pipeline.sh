#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob
kraken2_RUN_ID=$(date +%Y%m%d_%H%M%S)
qc_RUN_ID=$(date +%Y%m%d_%H%M%S)
# Load config
source config/config.sh

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

# Kraken2
mkdir -p "${kraken2_db_dir}"
## output 
mkdir -p "${kraken2_LOGDIR}"
mkdir -p "${kraken2_output_dir}"
mkdir -p "${kraken2_reports_dir}"
mkdir -p "${kraken2_classified_reads_dir}"

mkdir -p "${bracken_output_dir}"
mkdir -p "${bracken_reports_dir}"
mkdir -p "${bracken_krona_txt_dir}"
mkdir -p "${krona_html_dir}"

# mkdir -p results/1-kraken2/2-maaslin2
job_id1=$(sbatch --parsable --export=ALL code/slurm-scripts/get-reference-data.slurm)
job_id2=$(sbatch --parsable --dependency=afterok:$job_id1 --export=ALL code/slurm-scripts/get-toga-data.slurm)
echo "Pipeline submitted"
echo " step 1: ${job_id1}"
echo " step 2: ${job_id2}"
