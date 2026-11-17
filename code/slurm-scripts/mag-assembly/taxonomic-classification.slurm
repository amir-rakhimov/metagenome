#!/usr/bin/env bash
#SBATCH -t 360:00:00
#SBATCH -N 1
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -J taxonomic-classification
#SBATCH --output jobreports/%x_%j-out.txt
#SBATCH --error jobreports/%x_%j-out.txt
#I am requesting 4 nodes containing 40 CPUs, with 8 GB memory per CPU. Total: 320 GB
set -euo pipefail
shopt -s nullglob
source config/bash/config.sh
echo "${mag_RUN_ID}"
bash code/bash-scripts/mag-assembly/004-taxonomic-classification.sh
