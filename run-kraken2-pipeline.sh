#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob
kraken2_RUN_ID=$(date +%Y%m%d_%H%M%S)
# Load config
source config/bash/config.sh


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
job_id1=$(sbatch --parsable --export=ALL code/slurm-scripts/kraken2-pipeline/001-kraken2-build.slurm)
job_id2=$(sbatch --parsable --dependency=afterok:$job_id1 --export=ALL code/slurm-scripts/kraken2-pipeline/002-kraken2-pipeline.slurm)
echo "Pipeline submitted"
echo " step 1: ${job_id1}"
echo " step 2: ${job_id2}"
