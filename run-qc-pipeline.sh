#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob
# Load config
source config/bash/config.sh
SAMPLES_PATHS_RAW="${1}"
mkdir -p jobreports
## data
mkdir -p "${fastq_dir}"
mkdir -p $HOME/common_data
mkdir -p $HOME/common_data/reference_genomes

# QC and decontamination
mkdir -p "${fastqc_raw_dir}"
mkdir -p "${multiqc_raw_dir}"
mkdir -p "${fastqc_trimmed_dir}"
mkdir -p "${multiqc_trimmed_dir}"
mkdir -p "${cutadapt_output_dir}"

mkdir -p "${bowtie2_sam_dir}"
mkdir -p "${bowtie2_bam_dir}"
mkdir -p "${bowtie2_filtered_bam_dir}"
mkdir -p "${bowtie2_sorted_bam_dir}"
mkdir -p "${fastqc_decontam_dir}"
mkdir -p "${multiqc_decontam_dir}"
mkdir -p "${bowtie2_decontam_fastq_dir}"

# Tanes
mkdir -p "${SUNBEAM_OUTDIR}"
mkdir -p "${TRIMMOMATIC_OUTDIR}"
mkdir -p "${KOMPLEXITY_OUTDIR}"
mkdir -p "${BWA_INDICES_DIR}"
mkdir -p "${BWA_SAM_DIR}"


# job_id1=$(sbatch --parsable --export=ALL code/slurm-scripts/qc-pipeline/qc-commands-with-bowtie2.slurm "${SAMPLES_PATHS_RAW}")
job_id1=$(sbatch --parsable --export=ALL code/slurm-scripts/qc-pipeline/qc-tanes.slurm "${SAMPLES_PATHS_RAW}")
echo "Pipeline submitted"
echo " step 1: ${job_id1}"
