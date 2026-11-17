#!/usr/bin/env bash
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate mag_assembly-tools
set -euo pipefail
shopt -s nullglob
source config/bash/config.sh
SAMPLES_TSV="${1}"
# When nullglob is enabled, if a glob pattern does not match any files,
# it expands to nothing (an empty string) instead of returning the pattern itself.
# So, if no matches are found, the script will skip the file

# This script evaluates quality of assembled bins from MetaBAT2 using MetaQUAST and CheckM2.
# First, we use MetaQUAST (but it's not used in the subsequent analysis). Then, we run 
# CheckM2 on bins. Only high-quality MAGs (>=90% completeness, <=5% contamination)
# are used for the subsequent analysis. The statistics from CheckM2 is used for bin
# dereplication with dRep to retain unique MAGs. dRep is run at 95% and 99% ANI.
# The final output: CheckM2 statistics and dereplicated MAGs (copies of fasta files from MetaBAT2).
echo "$(date +"%F %H:%M:%S")"

# Function that checks if a required file exists or not
require() { [[ -s "$1" ]] || { echo "ERROR: required file missing: $1" >&2; exit 1; }; }
# Function that marks a sample as done (in a hidden directory)
is_done()   { [[ -e "${mark_dir}/${2}.${1}.done" ]]; }   # is_done STEP SAMPLE
mark_done() { touch "${mark_dir}/${2}.${1}.done"; }      # mark_done STEP SAMPLE

# 1. Read the input file list
[[ -s "${SAMPLES_TSV}" ]] || { echo "ERROR: metadata $SAMPLES_TSV missing or empty" >&2; exit 1; }
SAMPLES=()
while IFS=$'\t' read -r sample_id \
    || [[ -n "${sample_id}" ]]
do
    [[ "${sample_id}" == "sample_id" || -z "${sample_id}" ]] && continue
    SAMPLES+=("${sample_id}"); 
done <  "${SAMPLES_TSV}"

(( ${#SAMPLES[@]} > 0 )) || { echo "ERROR: no samples found" >&2; exit 1; } # Counts the number of elements in the final array of samples. If
                                                    # the array is empty, raises an error
echo $(date +"%F %H:%M:%S")
# 1. QC MetaBAT2 bins with QUAST
echo "$(date +"%F %H:%M:%S")"
mark_dir="${mag_qc_logs_dir}/.done"
mkdir -p "${mark_dir}"
for SAMPLE_ID in "${SAMPLES[@]}" 
do
    if is_done metaquast_metabat2 "${SAMPLE_ID}"; then
        echo "Skipping MetaQuast for bins: ${SAMPLE_ID}"; continue
    fi
    echo "Running MetaQuast on ${SAMPLE_ID} bins" 
    echo "$(date +"%F %H:%M:%S")"
    "${metaquast_script_dir}"/metaquast.py \
        "${metabat2_output_dir}"/"${SAMPLE_ID}"_bins/"${SAMPLE_ID}"_bin* \
        -o "${metaquast_metabat2_output_dir}"/metaquast_metabat2_"${SAMPLE_ID}" \
        --threads "${nthreads_qc}" \
        2>&1 |tee "${mag_qc_logs_dir}"/"${SAMPLE_ID}"_metaquast_metabat2.log;

    mark_done metaquast_metabat2 "${SAMPLE_ID}";
done

# 2. CheckM2
conda activate checkm2-env
echo "$(date +"%F %H:%M:%S")"
mark_dir="${mag_qc_logs_dir}/.done"
echo "Running CheckM2"
for SAMPLE_ID in "${SAMPLES[@]}" 
do
    if is_done checkm2 "${SAMPLE_ID}"; then
        echo "Skipping CheckM2 for bins: ${SAMPLE_ID}"; continue
    fi
    echo "Running CheckM2 on ${SAMPLE_ID}" 
    checkm2 predict --threads "${nthreads_qc_sort}" \
        -x fa \
        --database_path "${project_home_dir}"/data/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd \
        --input "${metabat2_output_dir}"/"${SAMPLE_ID}"_bins \
        --output-directory "${checkm2_output_dir}"/"${SAMPLE_ID}"_checkm2 \
        2>&1 |tee "${mag_qc_logs_dir}"/"${SAMPLE_ID}"_checkm2.log;
    mark_done checkm2 "${SAMPLE_ID}";
done
# Usage: checkm2 predict --threads 30 --input <folder_with_bins> --output-directory <output_folder> 
conda deactivate

echo "$(date +"%F %H:%M:%S")"
conda activate mag_assembly-tools
# 2.1 Merge CheckM2 quality reports
first_file=true
for file in "${checkm2_output_dir}"/*_checkm2/quality_report.tsv; do
    if [ "${first_file}" = true ]; then
        cat "${file}"          # print header + data
        first_file=false
    else
        tail -n +2 "${file}"   # skip header, print only data
    fi
done > "${checkm2_output_dir}"/quality_reports_merged.tsv

# 2.2 1347 MAGs in total
awk 'NR > 1 { count++ }
    END {print "There are " count " MAGs in total." }
    '  "${checkm2_output_dir}"/quality_reports_merged.tsv 

# 2.3 Select high-quality MAGs (>=90% completeness, <=5% contamination): 319
# Completeness is column #2, contamination is column #3
awk ' NR > 1 && $2 >= 90 && $3 <= 5 { count++ }
    END {
        print "There are " count " MAGs with completeness >= 90% and contamination <= 5%."
    }
    ' "${checkm2_output_dir}"/quality_reports_merged.tsv

# Subset the high-quality MAGs and write into a separate table
awk 'NR==1 || ($2 >= 90 && $3 <= 5)' \
    "${checkm2_output_dir}"/quality_reports_merged.tsv > \
    "${checkm2_output_dir}"/high_quality_mags.tsv

# 2.4 Prepare the high-quality MAG information for dRep.
# The file must be in .csv format and have the columns "genome"(basename of
# .fasta file of that genome), "completeness", and "contamination".
#  Columns “completeness” and “contamination” should be 0-100, and “genome” is the filename of the genome.
awk -F"\t" -v OFS="," 'FNR==1  {printf  "%s,%s,%s", "genome", "completeness", "contamination"; 
        print "";} 
        FNR>1 {print $1".fa",$2,$3} ' \
    "${checkm2_output_dir}"/high_quality_mags.tsv > \
    "${checkm2_output_dir}"/high_quality_mags_drep_info.csv

# 2.5 Write high-quality bin paths into a separate file
awk 'NR>1 {print $1}' "${checkm2_output_dir}"/high_quality_mags.tsv |\
while read -r name; do
    bin_subdir="${name%_bin.*}_bins"
    bin_path="${metabat2_output_dir}/${bin_subdir}/${name}.fa"
    if [ -f "${bin_path}" ]; then
        echo "${bin_path}"
    else
        echo " ${bin_path}"
    fi
done >  "${checkm2_output_dir}"/high_quality_mags_paths.txt \
      2> "${checkm2_output_dir}"/high_quality_mags_paths_missing.txt

# 3. Run drep to get unique MAGs
# dRep dereplicate -g path/to/genomes/*.fasta -p "${nthreads_qc}"  "${drep_output_dir}"
echo "$(date +"%F %H:%M:%S")"
# https://www.nature.com/articles/s42003-021-02827-2
# All MAGs were dereplicated at 99% ANI (equivalent to the strain level) and 95% ANI
#  (equivalent to the species level) using dRep (v2.6.2).

# https://link.springer.com/article/10.1186/s40168-022-01448-z
# All of the final bins were aggregated, and dRep ver. 1.1.2 software was then used 
# with the parameters ‘-p 16, -comp 80, -con 10, -str 100, -strW 0’ [27] to 
# remove duplicate bins. Then, dRep was used with the secondary clustering at the 
# threshold of 99% ANI with at least 25% overlap between genomes.

# 3.1 dRep at species level: secondary clustering threshold of 95% ANI
dRep dereplicate \
    -g "${checkm2_output_dir}"/high_quality_mags_paths.txt \
    -p "${nthreads_qc}" \
    --genomeInfo "${checkm2_output_dir}"/high_quality_mags_drep_info.csv \
    -pa 0.95 \
    -sa 0.95 \
    -comp 80 \
    -con 10 \
    -strW 0 \
    -nc 0.25 \
    -cm larger \
    -d \
    "${drep_output_dir}"/sa_95perc \
    2>&1 |tee "${mag_qc_logs_dir}"/drep_sa_95perc.log
echo "$(date +"%F %H:%M:%S")"

# 3.2 dRep at strain level: secondary clustering threshold of 99% ANI
dRep dereplicate \
    -g "${checkm2_output_dir}"/high_quality_mags_paths.txt \
    -p "${nthreads_qc}" \
    --genomeInfo "${checkm2_output_dir}"/high_quality_mags_drep_info.csv \
    -pa 0.95 \
    -sa 0.99 \
    -comp 80 \
    -con 10 \
    -strW 0 \
    -nc 0.25 \
    -cm larger \
    -d \
    "${drep_output_dir}"/sa_99perc \
    2>&1 |tee "${mag_qc_logs_dir}"/drep_sa_99perc.log
echo "$(date +"%F %H:%M:%S")"

# 3.3 Filter the dereplicated genomes to keep high-quality ones
dereplicated_mags_dir="${drep_output_dir}"/sa_95perc/dereplicated_genomes
checkm2_high_qual_mags="${checkm2_output_dir}"/high_quality_mags.tsv
drep_selected_high_qual_mags="${checkm2_output_dir}"/drep_selected_bins.txt
drep_missing_high_qual_mags="${mag_qc_logs_dir}"/drep_missing_files.log
drep_ani_perc=95

# Clear output files
> "${drep_selected_high_qual_mags}"
> "${drep_missing_high_qual_mags}"

# Read and process each line
if [[ ! -s "${drep_selected_high_qual_mags}" ]]; then
    while IFS=$'\t' read -r line; do
        # Extract the first column (filename)
        bin_id=$(echo "$line" | cut -f1)
        file=$(echo "${bin_id}".fa )
        # Check for existence in the target directory
        if [ -e "${dereplicated_mags_dir}/$file" ]; then
            echo -e "${bin_id}" >> "${drep_selected_high_qual_mags}"
        else
            echo "${bin_id}" >> "${drep_missing_high_qual_mags}"
        fi
    done < "${checkm2_high_qual_mags}"
fi


# ${bash_utils_scripts_dir}/filter-drep-mags.sh \
#  <input_dereplicated_mags_dir> \
#  <input_checkm2_high_qual_mags> \
#  <output_drep_selected_high_qual_mags> \
#  <output_drep_missing_high_qual_mags>

# 128 high-quality derpelicated MAGs at >95% level
# 192 MAGs removed (don't forget the header line)
bash "${bash_utils_scripts_dir}/filter-drep-mags.sh" \
    "${drep_output_dir}"/sa_95perc/dereplicated_genomes \
    "${checkm2_output_dir}"/high_quality_mags.tsv \
    "${checkm2_output_dir}"/drep_selected_bins_sa_95perc.txt \
    "${mag_qc_logs_dir}"/drep_missing_files_sa_95perc.log

# 138 high-quality derpelicated MAGs at >95% level
# 182 MAGs removed (don't forget the header line)
bash "${bash_utils_scripts_dir}/filter-drep-mags.sh" \
    "${drep_output_dir}"/sa_99perc/dereplicated_genomes \
    "${checkm2_output_dir}"/high_quality_mags.tsv \
    "${checkm2_output_dir}"/drep_selected_bins_sa_99perc.txt \
    "${mag_qc_logs_dir}"/drep_missing_files_sa_99perc.log

