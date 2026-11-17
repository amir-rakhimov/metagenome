#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob
# Load config
source config/bash/config.sh
SAMPLES_TSV="${1}"
mkdir -p jobreports
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

# Create metadata 
#  paste <(basename -s _trim_decontam_R1.fastq.gz data/2-bowtie2/decontam_fastq/*R1*gz) \
#     <(ls $PWD/data/2-bowtie2/decontam_fastq/*R1.*) \
#     <(ls $PWD/data/2-bowtie2/decontam_fastq/*R2.*) > data/metadata/yasuda/sample-paths-decontam.tsv
# # add header
# sed -i '1s/^/sample_id\tfwd_absolute_path\trev_absolute_path\n/' data/metadata/yasuda/sample-paths-decontam.tsv

 
 # mkdir -p results/1-kraken2/2-maaslin2
job_id1=$(sbatch --parsable --export=ALL code/slurm-scripts/kraken2-pipeline/001-kraken2-build.slurm)
job_id2=$(sbatch --parsable --dependency=afterok:$job_id1 --export=ALL code/slurm-scripts/kraken2-pipeline/002-kraken2-pipeline.slurm "${SAMPLES_TSV}")
echo "Pipeline submitted"
echo " step 1: ${job_id1}"
echo " step 2: ${job_id2}"
