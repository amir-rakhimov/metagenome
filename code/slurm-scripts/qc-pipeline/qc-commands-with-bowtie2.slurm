#!/usr/bin/env bash
#SBATCH -t 360:00:00
#SBATCH -N 4
#SBATCH -n 12
#SBATCH --mem-per-cpu 8g
#SBATCH -J qc-commands-with-bowtie2
#SBATCH --output jobreports/%x_%j-out.txt
#SBATCH --error jobreports/%x_%j-out.txt
# %x_%j job name and job ID
set -euo pipefail
shopt -s nullglob
source config/bash/config.sh
echo "${qc_RUN_ID}"
bash code/bash-scripts/qc-pipeline/001-qc-commands-with-bowtie2.sh

