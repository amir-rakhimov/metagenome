#!/usr/bin/env bash
#SBATCH -t 360:00:00
#SBATCH -N 4
#SBATCH -n 49
#SBATCH --mem-per-cpu 8g
#SBATCH -J kraken2-build
#SBATCH --output jobreports/%x_%j-out.txt
#SBATCH --error jobreports/%x_%j-out.txt
set -euo pipefail
source config/config.sh
bash code/bash-scripts/kraken2-pipeline/001-kraken2-build.sh

