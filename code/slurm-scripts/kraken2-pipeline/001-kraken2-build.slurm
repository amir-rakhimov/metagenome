#!/usr/bin/env bash
#SBATCH -t 360:00:00
#SBATCH -N 1
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=24
#SBATCH --mem=320G
#SBATCH -J kraken2-build
#SBATCH --output jobreports/%x_%j-out.txt
#SBATCH --error jobreports/%x_%j-out.txt
# %x_%j job name and job ID
set -euo pipefail
shopt -s nullglob
source config/bash/config.sh
echo "${kraken2_RUN_ID}"
bash code/bash-scripts/kraken2-pipeline/001-kraken2-build.sh

