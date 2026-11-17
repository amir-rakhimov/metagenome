#!/usr/bin/env bash
#SBATCH -t 360:00:00
#SBATCH -N 1
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=24
#SBATCH --mem=320G
#SBATCH -J kraken2-pipeline
#SBATCH --output jobreports/%x_%j-out.txt
#SBATCH --error jobreports/%x_%j-out.txt
#I am requesting 1 node, 1 task, 49 CPUs, with 8 GB memory per CPU. Total: 392 GBs
set -euo pipefail
shopt -s nullglob
source config/bash/config.sh
echo "${kraken2_RUN_ID}"
SAMPLES_TSV="${1}"
bash code/bash-scripts/kraken2-pipeline/002-kraken2-pipeline.sh "${SAMPLES_TSV}"

