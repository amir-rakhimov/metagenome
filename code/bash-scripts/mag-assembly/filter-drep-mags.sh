#!/bin/bash

# Exit on error
set -euo pipefail

# === USAGE EXAMPLE ===
# ./filter-drep-mags.sh <input_dereplicated_mags_dir> \
#  <input_checkm2_high_qual_mags> <output_drep_selected_high_qual_mags> \
#  <output_drep_missing_high_qual_mags>

# === INPUT ARGUMENTS ===
# drep_output_dir="$1"
# drep_date_time="$2"
# checkm2_output_dir="$3"
# checkm2_date_time="$4"
# drep_ani_percentage="$5"

# === DERIVED PATHS ===
# input_dereplicated_mags_dir="${drep_output_dir}/${drep_date_time}/dereplicated_genomes"
# input_checkm2_high_qual_mags="${checkm2_output_dir}/${checkm2_date_time}_high_quality_mags.tsv"
# output_drep_selected_high_qual_mags="${checkm2_output_dir}/${checkm2_date_time}_drep_selected_bins.txt"
# output_drep_missing_high_qual_mags="${checkm2_output_dir}/${checkm2_date_time}_drep_missing_files.log"
# === INPUT ARGUMENTS ===
input_dereplicated_mags_dir="$1"
input_checkm2_high_qual_mags="$2"
output_drep_selected_high_qual_mags="$3"
output_drep_missing_high_qual_mags="$4"

# === OUTPUT HEADS-UP ===
echo "Filtering MAGs from: $input_checkm2_high_qual_mags"
echo "Checking files in:   $input_dereplicated_mags_dir"
echo "Writing to:"
echo "  - Existing MAGs: $output_drep_selected_high_qual_mags"
echo "  - Missing MAGs:  $output_drep_missing_high_qual_mags"
echo

# === PREPARE OUTPUT FILES ===
> "$output_drep_selected_high_qual_mags"
> "$output_drep_missing_high_qual_mags"

# === SKIP HEADER AND FILTER ===
# tail -n +2 "$input_checkm2_high_qual_mags" | while IFS=$'\n' read -r line; do
#     bin_id=$(echo "$line" | cut -f1)
#     file="${bin_id}.fa"

#     if [ -e "$input_dereplicated_mags_dir/$file" ]; then
#         echo "$bin_id" >> "$output_drep_selected_high_qual_mags"
#     else
#         echo "$bin_id" >> "$output_drep_missing_high_qual_mags"
#     fi
# done

while IFS=$'\t' read -r line; do
    # Extract the first column (filename)
    bin_id=$(echo "$line" | cut -f1)
    file=$(echo "$bin_id".fa )
    # Check for existence in the target directory
    if [ -e "$input_dereplicated_mags_dir/$file" ]; then
        echo -e "$bin_id" >> "$output_drep_selected_high_qual_mags"
    else
        echo "$bin_id" >> "$output_drep_missing_high_qual_mags"
    fi
done < "$input_checkm2_high_qual_mags"

echo "Filtering complete."

# How many MAGs remained (no header)
awk ' { count++ }
  END {print "There are " count " high-quality derpelicated MAGs." }
  '  "$output_drep_selected_high_qual_mags"

# 187 MAGs removed (don't forget the header line)
awk 'NR>1 { count++ }
  END {print "Removed " count " MAGs." }
  '  "$output_drep_missing_high_qual_mags"