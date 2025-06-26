#!/usr/bin/env bash
#SBATCH -t 360:00:00
#SBATCH -N 4
#SBATCH -n 40
#SBATCH --mem-per-cpu 8g
#SBATCH -J 20250625_14-52-mag-qc
#SBATCH --output jobreports/20250625_14-52-mag-qc-out.txt
#SBATCH --error jobreports/20250625_14-52-mag-qc-out.txt
#I am requesting 4 nodes containing 40 CPUs, with 8 GB memory per CPU. Total: 320 GB
source ~/miniconda3/etc/profile.d/conda.sh
shopt -s nullglob
# When nullglob is enabled, if a glob pattern does not match any files,
# it expands to nothing (an empty string) instead of returning the pattern itself.
# So, if no matches are found, the script will skip the file

# Metaquast is a metagenome assembly evaluation tool
# Metaquast is also used to evaluate the quality of bins
# CheckM2 assesses the quality of bins
# dRep dereplicates MAGs to retain unique ones
date_var=$(date -I|sed 's/-//g')
time_var=$(date +%T |sed 's/:/_/g' )
date_time=${date_var}_${time_var}
start_date_time=$(date +"%F %H:%M:%S")
# For binning with MetaBAT2
# metabat2_date_time=20250417_23_21_29
metabat2_date_time=20250612_13_37_47
# metaquast_metabat2_date_time=20250417_23_21_29
metaquast_metabat2_date_time="${date_time}"
# checkm2_date_time=20250501_07_53_55
checkm2_date_time=20250619_05_47_09
drep_date_time="${date_time}"
nthreads=40
nthreads_sort=35
mem_req=8G
mem_req_sort=4G
project_home_dir=~/projects/metagenome
output_dir=~/projects/metagenome/output/mag_assembly
metaquast_output_dir=output/mag_assembly/metaquast_output
metaquast_script_dir=~/quast-5.2.0
metabat2_output_dir=output/mag_assembly/metabat2_output
checkm2_output_dir=output/mag_assembly/checkm2_output
drep_output_dir=output/mag_assembly/drep_output

cd ${project_home_dir}
mkdir -p output/mag_assembly/checkm2_output
mkdir -p output/mag_assembly/metaquast_output

# 0. Show the current time for logging
echo "${start_date_time}"

# 6.4 QC MetaBAT2 bins with QUAST
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
# for FILE_DIR in "${metabat2_output_dir}"/"${metabat2_date_time}"*bins
# # # for FILE_DIR in "${metabat2_output_dir}"/"${metabat2_date_time}"_2D14_bins \
# # #  "${metabat2_output_dir}"/"${metabat2_date_time}"_G14_bins \
# # #  "${metabat2_output_dir}"/"${metabat2_date_time}"_H15_bins 
# do 
#  SAMPLE=$(echo "${FILE_DIR}" | sed "s/${metabat2_date_time}_//" | sed "s/_bins//")
#  base_name=$(basename "$SAMPLE" )
#  echo "Running MetaQuast on ${base_name} bins" 
#  "${metaquast_script_dir}"/metaquast.py \
#     "${metabat2_output_dir}"/"${metabat2_date_time}"_"${base_name}"_bins/"${metabat2_date_time}"_"${base_name}"_bin* \
#     -o "${metaquast_output_dir}"/"${metaquast_metabat2_date_time}"_metaquast_metabat2_"${base_name}" \
#     --threads "${nthreads}" \
#     2>&1 |tee "${metaquast_output_dir}"/"${metaquast_metabat2_date_time}"_"${base_name}"_metaquast_metabat2_report.txt;
# done

# 6.6 CheckM2
# conda activate checkm2-env
# intermediate_date_time=$(date +"%F %H:%M:%S")
# echo "${intermediate_date_time}"
# echo "Running CheckM2"
# for FILE_DIR in "${metabat2_output_dir}"/"${metabat2_date_time}"*bins
# do 
#  SAMPLE=$(echo "${FILE_DIR}" | sed "s/${metabat2_date_time}_//"| sed "s/_bins//")
#  base_name=$(basename "$SAMPLE" )
#  echo "Running CheckM2 on ${base_name}" 
#  checkm2 predict --threads "${nthreads_sort}" \
#    -x fa \
#    --database_path "${project_home_dir}"/data/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd \
#    --input "${metabat2_output_dir}"/"${metabat2_date_time}"_"${base_name}"_bins \
#    --output-directory "${checkm2_output_dir}"/"${checkm2_date_time}"_"${base_name}"_checkm2 \
#    2>&1 |tee "${checkm2_output_dir}"/"${checkm2_date_time}"_"${base_name}"_checkm2_report.txt;
# done
# # Usage: checkm2 predict --threads 30 --input <folder_with_bins> --output-directory <output_folder> 
# conda deactivate

intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
conda activate mag_assembly-tools
# Merge CheckM2 quality reports
# first_file=true
# for file in "${checkm2_output_dir}"/"${checkm2_date_time}"*_checkm2/quality_report.tsv; do
#   if [ "${first_file}" = true ]; then
#     cat "${file}"          # print header + data
#     first_file=false
#   else
#     tail -n +2 "${file}"   # skip header, print only data
#   fi
# done > "${checkm2_output_dir}"/"${checkm2_date_time}"_quality_reports_merged.tsv

# # 1347 MAGs in total
# awk 'NR > 1 { count++ }
#   END {print "There are " count " MAGs in total." }
#   '  "${checkm2_output_dir}"/"${checkm2_date_time}"_quality_reports_merged.tsv 

# Select high-quality MAGs (>=90% completeness, <=5% contamination): 313
# Completeness is column #2, contamination is column #3
# awk ' NR > 1 && $2 >= 90 && $3 <= 5 { count++ }
#   END {
#     print "There are " count " MAGs with completeness >= 90% and contamination <= 5%."
#   }
#  ' "${checkm2_output_dir}"/"${checkm2_date_time}"_quality_reports_merged.tsv

# Subset the high-quality MAGs and write into a separate table
# awk 'NR==1 || ($2 >= 90 && $3 <= 5)' \
#  "${checkm2_output_dir}"/"${checkm2_date_time}"_quality_reports_merged.tsv > \
#  "${checkm2_output_dir}"/"${checkm2_date_time}"_high_quality_mags.tsv

# Write high-quality bin paths into a separate file
# awk 'NR>1 {print $1}' "${checkm2_output_dir}"/"${checkm2_date_time}"_high_quality_mags.tsv |\
#  while read -r name; do
#   bin_subdir="${name%_bin.*}_bins"
#   bin_path="${metabat2_output_dir}/${bin_subdir}/${name}.fa"
#   if [ -f "${bin_path}" ]; then
#     echo "${bin_path}"
#   else
#     echo " ${bin_path}"
#   fi
#  done >  "${checkm2_output_dir}"/"${checkm2_date_time}"_high_quality_mags_paths.txt \
#       2> "${checkm2_output_dir}"/"${checkm2_date_time}"_high_quality_mags_paths_missing.txt

# Run drep to get unique MAGs
# dRep dereplicate -g path/to/genomes/*.fasta -p "${nthreads}"  "${drep_output_dir}"
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
# https://www.nature.com/articles/s42003-021-02827-2
# All MAGs were dereplicated at 99% ANI (equivalent to the strain level) and 95% ANI
#  (equivalent to the species level) using dRep (v2.6.2).

# https://link.springer.com/article/10.1186/s40168-022-01448-z
# All of the final bins were aggregated, and dRep ver. 1.1.2 software was then used 
# with the parameters ‘-p 16, -comp 80, -con 10, -str 100, -strW 0’ [27] to 
# remove duplicate bins. Then, dRep was used with the secondary clustering at the 
# threshold of 99% ANI with at least 25% overlap between genomes.

mkdir -p output/mag_assembly/drep_output/"${drep_date_time}"_sa_95perc
mkdir -p output/mag_assembly/drep_output/"${drep_date_time}"_sa_99perc
# dRep at species level: secondary clustering threshold of 95% ANI
dRep dereplicate \
 -g "${checkm2_output_dir}"/"${checkm2_date_time}"_high_quality_mags_paths.txt \
 -p "${nthreads_sort}" \
 -pa 0.95 \
 -sa 0.95 \
 -comp 80 \
 -con 10 \
 -strW 0 \
 -nc 0.25 \
 -cm larger \
 -d \
 "${drep_output_dir}"/"${drep_date_time}"_sa_95perc \
 2>&1 |tee "${drep_output_dir}"/"${drep_date_time}"_drep_sa_95perc_report.txt
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"

# dRep at strain level: secondary clustering threshold of 99% ANI
dRep dereplicate \
 -g "${checkm2_output_dir}"/"${checkm2_date_time}"_high_quality_mags_paths.txt \
 -p "${nthreads_sort}" \
 -pa 0.95 \
 -sa 0.99 \
 -comp 80 \
 -con 10 \
 -strW 0 \
 -nc 0.25 \
 -cm larger \
 -d \
 "${drep_output_dir}"/"${drep_date_time}"_sa_99perc \
 2>&1 |tee "${drep_output_dir}"/"${drep_date_time}"_drep_sa_99perc_report.txt
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
# Default:
#  -comp 75 \
#  -con 25 \
#  -strW 1 \
#  -nc 0.1 \
# -pa P_ANI, --P_ani P_ANI
#                       ANI threshold to form primary (MASH) clusters
#                       (default: 0.9)
# -sa S_ANI, --S_ani S_ANI
#                       ANI threshold to form secondary clusters (default:
#                       0.95)
# -nc COV_THRESH, --cov_thresh COV_THRESH
#                           Minmum level of overlap between genomes when doing
#                           secondary comparisons (default: 0.1)
# -strW STRAIN_HETEROGENEITY_WEIGHT, --strain_heterogeneity_weight STRAIN_HETEROGENEITY_WEIGHT
#                       strain heterogeneity weight (default: 1)
# -cm {total,larger}, --coverage_method {total,larger}
#                       Method to calculate coverage of an alignment
#                       (for ANIn/ANImf only; gANI and fastANI can only do larger method)
#                       total   = 2*(aligned length) / (sum of total genome lengths)
#                       larger  = max((aligned length / genome 1), (aligned_length / genome2))
#                        (default: larger)

# Filter the dereplicated genomes to keep high-quality ones
# dereplicated_mags_dir="${drep_output_dir}"/"${drep_date_time}"_sa_95perc/dereplicated_genomes
# checkm2_high_qual_mags="${checkm2_output_dir}"/"${checkm2_date_time}"_high_quality_mags.tsv
# drep_selected_high_qual_mags="${checkm2_output_dir}"/"${checkm2_date_time}"_drep_selected_bins.txt
# drep_missing_high_qual_mags="${checkm2_output_dir}"/"${checkm2_date_time}"_drep_missing_files.log
# drep_ani_perc=95

# # Clear output files
# > "$drep_selected_high_qual_mags"
# > "$drep_missing_high_qual_mags"

# # Read and process each line
# while IFS=$'\t' read -r line; do
#     # Extract the first column (filename)
#     bin_id=$(echo "$line" | cut -f1)
#     file=$(echo "$bin_id".fa )
#     # Check for existence in the target directory
#     if [ -e "$dereplicated_mags_dir/$file" ]; then
#         echo -e "$bin_id" >> "$drep_selected_high_qual_mags"
#     else
#         echo "$bin_id" >> "$drep_missing_high_qual_mags"
#     fi
# done < "$checkm2_high_qual_mags"

# code/bash-scripts/mag-assembly/filter-drep-mags.sh \
#  <input_dereplicated_mags_dir> \
#  <input_checkm2_high_qual_mags> \
#  <output_drep_selected_high_qual_mags> \
#  <output_drep_missing_high_qual_mags>

# 126 high-quality derpelicated MAGs
# 187 MAGs removed (don't forget the header line)
./code/bash-scripts/mag-assembly/filter-drep-mags.sh \
 "${drep_output_dir}"/"${drep_date_time}"_sa_95perc/dereplicated_genomes \
 "${checkm2_output_dir}"/"${checkm2_date_time}"_high_quality_mags.tsv \
 "${checkm2_output_dir}"/"${checkm2_date_time}"_drep_selected_bins_sa_95perc.txt \
 "${checkm2_output_dir}"/"${checkm2_date_time}"_drep_missing_files_sa_95perc.log


./code/bash-scripts/mag-assembly/filter-drep-mags.sh \
 "${drep_output_dir}"/"${drep_date_time}"_sa_99perc/dereplicated_genomes \
 "${checkm2_output_dir}"/"${checkm2_date_time}"_high_quality_mags.tsv \
 "${checkm2_output_dir}"/"${checkm2_date_time}"_drep_selected_bins_sa_99perc.txt \
 "${checkm2_output_dir}"/"${checkm2_date_time}"_drep_missing_files_sa_99perc.log

