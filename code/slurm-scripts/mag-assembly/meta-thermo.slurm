#!/usr/bin/env bash
#SBATCH -t 360:00:00
#SBATCH -J meta-thermo
#SBATCH --output jobreports/%x_%j-out.txt
#SBATCH --error jobreports/%x_%j-out.txt
set -euo pipefail
shopt -s nullglob
source config/bash/config.sh
echo "${RUN_ID}"
bash code/bash-scripts/mag-assembly/010-meta-thermo-script.sh