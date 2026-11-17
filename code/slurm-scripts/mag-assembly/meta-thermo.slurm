#!/usr/bin/env bash
#SBATCH -t 10:00:00
#SBATCH -J meta-thermo
#SBATCH --output jobreports/%x_%j-out.txt
#SBATCH --error jobreports/%x_%j-out.txt
set -euo pipefail
shopt -s nullglob
source config/bash/config.sh
echo "${mag_RUN_ID}"
SAMPLES_TSV="${1}"
bash code/bash-scripts/mag-assembly/010-meta-thermo-script.sh "${SAMPLES_TSV}"