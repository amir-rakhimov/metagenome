#!/usr/bin/env bash
#SBATCH -t 360:00:00
#SBATCH -J quantify-genes
#SBATCH --output jobreports/%x_%j-out.txt
#SBATCH --error jobreports/%x_%j-out.txt
set -euo pipefail
shopt -s nullglob
source config/bash/config.sh
echo "${RUN_ID}"
bash code/bash-scripts/mag-assembly/006-quantify-genes.sh