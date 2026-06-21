#!/usr/bin/env bash
#SBATCH -t 360:00:00
#SBATCH -N 4
#SBATCH -n 49
#SBATCH --mem-per-cpu 8g
#SBATCH -J kraken2-pipeline
#SBATCH --output jobreports/%x_%j-out.txt
#SBATCH --error jobreports/%x_%j-out.txt
#I am requesting 4 nodes containing 49 CPUs, with 8 GB memory per CPU. Total: 392 GB
set -euo pipefail
shopt -s nullglob
source config/bash/config.sh
echo "${RUN_ID}"
bash code/bash-scripts/kraken2-pipeline/002-kraken2-pipeline.sh

