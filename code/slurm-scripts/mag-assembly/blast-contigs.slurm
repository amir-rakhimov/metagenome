#!/usr/bin/env bash
#SBATCH -t 360:00:00
#SBATCH -N 4
#SBATCH -n 20
#SBATCH --mem-per-cpu 8g
#SBATCH -J blast-contigs
#SBATCH --output jobreports/%x_%j-out.txt
#SBATCH --error jobreports/%x_%j-out.txt
#I am requesting 4 nodes containing 20 CPUs, with 8 GB memory per CPU. Total: 160 GB
set -euo pipefail
shopt -s nullglob
source config/bash/config.sh
echo "${RUN_ID}"
bash code/bash-scripts/mag-assembly/008-blast-contigs.sh