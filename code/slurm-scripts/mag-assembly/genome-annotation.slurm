#!/usr/bin/env bash
#SBATCH -t 360:00:00
#SBATCH -N 1
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -J genome-annotation
#SBATCH --output jobreports/%x_%j-out.txt
#SBATCH --error jobreports/%x_%j-out.txt
#I am requesting 1 node with 1 task containing 40 CPUs, with 8 GB memory per CPU. Total: 320 GB
set -euo pipefail
shopt -s nullglob
source config/bash/config.sh
echo "${mag_RUN_ID}"
SAMPLES_TSV="${1}"
bash code/bash-scripts/mag-assembly/005-genome-annotation.sh "${SAMPLES_TSV}"