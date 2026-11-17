#!/usr/bin/env bash
#SBATCH -t 360:00:00
#SBATCH -N 1
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -J qc-commands-with-bowtie2
#SBATCH --output jobreports/%x_%j-out.txt
#SBATCH --error jobreports/%x_%j-out.txt
# %x_%j job name and job ID
set -euo pipefail
shopt -s nullglob
source config/bash/config.sh
echo "${qc_RUN_ID}"
SAMPLES_PATHS_RAW="${1}"
GENOMES_METADATA="{2}"
bash code/bash-scripts/qc-pipeline/001-qc-commands-with-bowtie2.sh "${SAMPLES_PATHS_RAW}" "${GENOMES_METADATA}"

