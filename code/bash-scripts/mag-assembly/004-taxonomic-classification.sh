#!/usr/bin/env bash
#SBATCH -t 360:00:00
#SBATCH -N 4
#SBATCH -n 40
#SBATCH --mem-per-cpu 8g
#SBATCH -J 20250707_11-27-taxonomic-classification
#SBATCH --output jobreports/20250707_11-27-taxonomic-classification-out.txt
#SBATCH --error jobreports/20250707_11-27-taxonomic-classification-out.txt
#I am requesting 4 nodes containing 40 CPUs, with 8 GB memory per CPU. Total: 320 GB
source ~/miniconda3/etc/profile.d/conda.sh
shopt -s nullglob
# When nullglob is enabled, if a glob pattern does not match any files,
# it expands to nothing (an empty string) instead of returning the pattern itself.
# So, if no matches are found, the script will skip the file

# This script assigns taxonomy to dereplicated genomes with GTDB-TK. GTDB-TK is run
# on 95% and 99% ANI bins. classify workflow is run separately because there's an issue with 
# pplacer. We also run CompareM on dereplicated genomes. 
# Main final output: taxonomic classification of bins as tab-separated files for bacteria and archaea.
date_var=$(date -I|sed 's/-//g')
time_var=$(date +%T |sed 's/:/_/g' )
date_time="${date_var}_${time_var}"
start_date_time=$(date +"%F %H:%M:%S")
# 0. Show the current time for logging
echo "${start_date_time}"

# 1. Activate conda environment
conda activate mag_assembly-tools

# 2. Classify high-quality dereplicated MAGs with GTDB-TK
# GTDB-Tk requires ~110G of external data that needs to be downloaded and unarchived
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
# 2.1 Download the reference database
download-db.sh
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz \
 -O ~/miniconda3/envs/mag_assembly-tools/share/gtdbtk-2.4.0/db/gtdbtk_data.tar.gz
tar -xvf "${GTDBTK_DATA_PATH}"/*tar.gz -C "${GTDBTK_DATA_PATH}"
# 2.2 Classify MAGs
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
echo "Running GTDB-TK"
# Fix pplacer locale error
export LC_ALL=C
echo "Running GTDB-TK classify_wf on dereplicated genomes (95% ANI threshold)"
gtdbtk classify_wf \
    --genome_dir "${drep_output_dir}"/sa_95perc/dereplicated_genomes/ \
    --out_dir "${gtdbtk_output_dir}"/sa_95perc_GTDBtk \
    -x fa \
    --skip_ani_screen \
    --cpus "${nthreads_gtdbtk}" \
    --pplacer_cpus "${nthreads_pplacer}"
    # --scratch_dir "${gtdbtk_output_dir}"/sa_95perc_GTDBtk \
 # Run classify separately because pplacer didn't recognize gz files
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"

echo "Running GTDB-TK classify on dereplicated genomes (95% ANI threshold)"
gunzip "${gtdbtk_output_dir}"/sa_95perc_GTDBtk/align/*.gz
gtdbtk classify --genome_dir "${drep_output_dir}"/sa_95perc/dereplicated_genomes/ \
    --out_dir "${gtdbtk_output_dir}"/sa_95perc_GTDBtk \
    --align_dir "${gtdbtk_output_dir}"/sa_95perc_GTDBtk \
    -x fa \
    --skip_ani_screen \
    --cpus "${nthreads_gtdbtk}" \
    --pplacer_cpus "${nthreads_pplacer}"
    # --scratch_dir "${gtdbtk_output_dir}"/sa_95perc_GTDBtk \
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"

echo "Running GTDB-TK classify_wf on dereplicated genomes (99% ANI threshold)"
gtdbtk classify_wf \
    --genome_dir "${drep_output_dir}"/sa_99perc/dereplicated_genomes/ \
    --out_dir "${gtdbtk_output_dir}"/sa_99perc_GTDBtk \
    -x fa \
    --skip_ani_screen \
    --cpus "${nthreads_gtdbtk}" \
    --pplacer_cpus "${nthreads_pplacer}"
    # --scratch_dir "${gtdbtk_output_dir}"/sa_99perc_GTDBtk \
 # Run classify separately because pplacer didn't recognize gz files
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
echo "Running GTDB-TK classify on dereplicated genomes (99% ANI threshold)"
gunzip "${gtdbtk_output_dir}"/sa_99perc_GTDBtk/align/*.gz
gtdbtk classify --genome_dir "${drep_output_dir}"/sa_99perc/dereplicated_genomes/ \
    --out_dir "${gtdbtk_output_dir}"/sa_99perc_GTDBtk \
    --align_dir "${gtdbtk_output_dir}"/sa_99perc_GTDBtk \
    -x fa \
    --skip_ani_screen \
    --cpus "${nthreads_gtdbtk}" \
    --pplacer_cpus "${nthreads_pplacer}"
    # --scratch_dir "${gtdbtk_output_dir}"/sa_99perc_GTDBtk \
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"

# In case you want to run on all MAGs
# for FILE_DIR in "${metabat2_output_dir}"/*bins
# do 
#  SAMPLE=$(echo "${FILE_DIR}" | sed "s/_bins//")
#  base_name=$(basename "$SAMPLE" )
#  echo "Running GTDB-TK classify_wf on ${base_name}"
#  gtdbtk classify_wf --genome_dir "${metabat2_output_dir}"/"${base_name}"_bins/ \
#     --out_dir "${gtdbtk_output_dir}"/"${base_name}"_GTDBtk \
#     -x fa \
#     --skip_ani_screen \
#     --cpus "${nthreads_gtdbtk}" \
#     --pplacer_cpus "${nthreads_pplacer}"
#  # Run classify separately because pplacer didn't recognize gz files
#  echo "Running GTDB-TK classify on ${base_name}"
#  gunzip "${gtdbtk_output_dir}"/"${base_name}"_GTDBtk/align/*.gz
#  gtdbtk classify --genome_dir "${metabat2_output_dir}"/"${base_name}"_bins/ \
#     --out_dir "${gtdbtk_output_dir}"/"${base_name}"_GTDBtk \
#     --align_dir "${gtdbtk_output_dir}"/"${base_name}"_GTDBtk \
#     -x fa \
#     --skip_ani_screen \
#     --cpus "${nthreads_gtdbtk}" \
#     --pplacer_cpus "${nthreads_pplacer}";
# done


# Merge GTDB-Tk data within each sample
for FILE_DIR in "${gtdbtk_output_dir}"/*_GTDBtk
do
 base_name=$(basename "${FILE_DIR}" |sed "s/_GTDBtk//")
 cp "${FILE_DIR}"/classify/gtdbtk.bac120.summary.tsv \
  "${gtdbtk_output_dir}"/"${base_name}"_classification_combined.tsv
 tail --lines=+2 "${FILE_DIR}"/gtdbtk.ar53.summary.tsv >> \
  "${gtdbtk_output_dir}"/"${base_name}"_classification_combined.tsv;
done


# CompareM is used to calculate the average amino acid identity (AAI) among the MAGs
# comparem aai_wf <input_files> <output_dir>
# comparem --cpus 32 aai_wf my_genomes aai_output
echo "Running CompareM on dereplicated genomes (95% ANI threshold)"
mkdir -p "${comparem_output_dir}"/aai_output_sa95perc
comparem  aai_wf \
 --cpus "${nthreads_comparem}" \
 "${drep_output_dir}"/sa_95perc/dereplicated_genomes \
 -x fa \
 "${comparem_output_dir}"/aai_output_sa95perc
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"

echo "Running CompareM on dereplicated genomes (99% ANI threshold)"
mkdir -p "${comparem_output_dir}"/aai_output_sa99perc
comparem  aai_wf \
 --cpus "${nthreads_comparem}" \
 "${drep_output_dir}"/sa_99perc/dereplicated_genomes \
 -x fa \
 "${comparem_output_dir}"/aai_output_sa99perc
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"

# Checking output/mag_assembly/gtdbtk_output/20250705_19_24_02_sa_95perc_classification_combined.tsv
# There are 128 MAGs in total.
# Checking output/mag_assembly/gtdbtk_output/20250705_19_24_02_sa_99perc_classification_combined.tsv
# There are 138 MAGs in total.
for FILE in "${gtdbtk_output_dir}"/*_classification_combined.tsv
do
 echo "Checking" "${FILE}"
 awk 'NR > 1 { count++ }
  END {print "There are " count " MAGs in total." }
  '  "${FILE}";
done