#!/usr/bin/env bash
#SBATCH -t 360:00:00
#SBATCH -N 4
#SBATCH -n 40
#SBATCH --mem-per-cpu 8g
#SBATCH -J 20250619_19-13-004-taxonomic-classification
#SBATCH --output jobreports/20250619_19-13-004-taxonomic-classification-out.txt
#SBATCH --error jobreports/20250619_19-13-004-taxonomic-classification-out.txt
#I am requesting 4 nodes containing 40 CPUs, with 8 GB memory per CPU. Total: 320 GB. 280 GB for GTDB-TK
source ~/miniconda3/etc/profile.d/conda.sh
shopt -s nullglob
# When nullglob is enabled, if a glob pattern does not match any files,
# it expands to nothing (an empty string) instead of returning the pattern itself.
# So, if no matches are found, the script will skip the file

# Metaquast is a metagenome assembly evaluation tool
# MetaBAT2 is a binning software: group contigs into bins (genomes)
# Metaquast is also used to evaluate the quality of bins
# CheckM2 assesses the quality of bins
# GTDB-TK assigns taxonomy to genomes
# Prodigal is a gene prediction software
# Prokka annotates contigs to find bacterial genes
# Metaeuk annotates contigs to find eukaryotic genes
date_var=$(date -I|sed 's/-//g')
time_var=$(date +%T |sed 's/:/_/g' )
date_time="${date_var}_${time_var}"
start_date_time=$(date +"%F %H:%M:%S")
# metabat2_date_time=20250417_23_21_29
metabat2_date_time=20250612_13_37_47
drep_date_time="${date_time}"
# gtdbtk_date_time=20250514_17_57_12
gtdbtk_date_time="${date_time}"
nthreads=40
nthreads_sort=35
mem_req=8G
mem_req_sort=4G
project_home_dir=~/projects/metagenome
output_dir=~/projects/metagenome/output/mag_assembly
metabat2_output_dir=output/mag_assembly/metabat2_output
drep_output_dir=output/mag_assembly/drep_output
gtdbtk_output_dir=output/mag_assembly/gtdbtk_output

cd ${project_home_dir}
mkdir -p output/mag_assembly/gtdbtk_output

# 0. Show the current time for logging
echo "${start_date_time}"

# 1. Activate conda environment
conda activate mag_assembly-tools

# 7. Classify with GTDB-TK
# GTDB-Tk requires ~110G of external data that needs to be downloaded and unarchived
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
# download-db.sh
# wget https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/auxillary_files/gtdbtk_package/full_package/gtdbtk_data.tar.gz \
#  -O ~/miniconda3/envs/mag_assembly-tools/share/gtdbtk-2.4.0/db/gtdbtk_data.tar.gz
# tar -xvf "${GTDBTK_DATA_PATH}"/*tar.gz -C "${GTDBTK_DATA_PATH}"
echo "Running GTDB-TK"
# Fix pplacer locale error
export LC_ALL=C
for FILE_DIR in "${metabat2_output_dir}"/"${metabat2_date_time}"*bins
do 
 SAMPLE=$(echo "${FILE_DIR}" | sed "s/${metabat2_date_time}_//" |sed "s/_bins//")
 base_name=$(basename "$SAMPLE" )
 echo "Running GTDB-TK classify_wf on ${base_name}"
 gtdbtk classify_wf --genome_dir "${metabat2_output_dir}"/"${metabat2_date_time}"_"${base_name}"_bins/ \
    --out_dir "${gtdbtk_output_dir}"/"${gtdbtk_date_time}"_"${base_name}"_GTDBtk \
    -x fa \
    --skip_ani_screen \
    --cpus "${nthreads_sort}" \
    --pplacer_cpus 1
 # Run classify separately because pplacer didn't recognize gz files
 echo "Running GTDB-TK classify on ${base_name}"
 #  gunzip "${gtdbtk_output_dir}"/"${gtdbtk_date_time}"_"${base_name}"_GTDBtk/align/*.gz
 gtdbtk classify --genome_dir "${metabat2_output_dir}"/"${metabat2_date_time}"_"${base_name}"_bins/ \
    --out_dir "${gtdbtk_output_dir}"/"${gtdbtk_date_time}"_"${base_name}"_GTDBtk \
    --align_dir "${gtdbtk_output_dir}"/"${gtdbtk_date_time}"_"${base_name}"_GTDBtk \
    -x fa \
    --skip_ani_screen \
    --cpus "${nthreads_sort}" \
    --pplacer_cpus 1;
done
# Run classify separately because pplacer didn't recognize gz files
for FILE_DIR in "${metabat2_output_dir}"/"${metabat2_date_time}"*bins
do 
 SAMPLE=$(echo "${FILE_DIR}" | sed "s/${metabat2_date_time}_//" |sed "s/_bins//")
 base_name=$(basename "$SAMPLE" )
 echo "Running GTDB-TK classify on ${base_name}"
#  gunzip "${gtdbtk_output_dir}"/"${gtdbtk_date_time}"_"${base_name}"_GTDBtk/align/*.gz
 gtdbtk classify --genome_dir "${metabat2_output_dir}"/"${metabat2_date_time}"_"${base_name}"_bins/ \
    --out_dir "${gtdbtk_output_dir}"/"${gtdbtk_date_time}"_"${base_name}"_GTDBtk \
    --align_dir "${gtdbtk_output_dir}"/"${gtdbtk_date_time}"_"${base_name}"_GTDBtk \
    -x fa \
    --skip_ani_screen \
    --cpus "${nthreads_sort}" \
    --pplacer_cpus 1;
done
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"



# # Merge GTDB-Tk data within each sample
# for FILE_DIR in "${gtdbtk_output_dir}"/"${gtdbtk_date_time}"*_GTDBtk/classify
# do
#  SAMPLE=$(echo "${FILE_DIR}" | sed "s/${gtdbtk_date_time}_//" |sed "s/_GTDBtk\/classify//")
#  base_name=$(basename "$SAMPLE" )

#  awk -F'\t' '{print $0}' "${FILE_DIR}"/gtdbtk.bac120.summary.tsv > \
#   "${gtdbtk_output_dir}"/"${gtdbtk_date_time}"_"${base_name}"_classification.tsv
#  awk -F'\t' 'NR>1{print $0}' "${FILE_DIR}"/gtdbtk.ar53.summary.tsv >> \
#   "${gtdbtk_output_dir}"/"${gtdbtk_date_time}"_"${base_name}"_classification.tsv;
# done

# # Merge all GTDB-Tk data from all samples into one table
# first_file=true
# for file in "${gtdbtk_output_dir}"/"${gtdbtk_date_time}"*_classification.tsv; do
#   if [ "${first_file}" = true ]; then
#     cat "${file}"         # print header + data
#     first_file=false
#   else
#     tail -n +2 "${file}"   # skip header, print only data
#   fi
# done > "${gtdbtk_output_dir}"/"${gtdbtk_date_time}"_classification_all.tsv


# # 1347 MAGs in total
# awk 'NR > 1 { count++ }
#   END {print "There are " count " MAGs in total." }
#   '  "${gtdbtk_output_dir}"/"${gtdbtk_date_time}"_classification_all.tsv

# # Keep only high-quality dereplicated MAGs (max 126)
# #  "${gtdbtk_output_dir}"/"${gtdbtk_date_time}"_drep_filtered_taxonomy.tsv is just a list!
# # checkm2_date_time=20250501_07_53_55 # old binning (gtdbtk_date_time=20250514_17_57_12)
# checkm2_date_time=20250619_05_47_09 # new binning
# checkm2_output_dir=output/mag_assembly/checkm2_output

# awk 'NR==FNR {a[$1]; next} FNR==1 || $1 in a' \
#  "${checkm2_output_dir}"/"${checkm2_date_time}"_drep_selected_bins.txt \
#  "${gtdbtk_output_dir}"/"${gtdbtk_date_time}"_classification_all.tsv \
#  > "${gtdbtk_output_dir}"/"${gtdbtk_date_time}"_drep_filtered_taxonomy.tsv
# # For a different column, replace $1 with $2, $3, etc.

# # 126 MAGs left
# awk 'NR > 1 { count++ }
#   END {print "There are " count " high-quality derpelicated MAGs." }
#   '  "${gtdbtk_output_dir}"/"${gtdbtk_date_time}"_drep_filtered_taxonomy.tsv
