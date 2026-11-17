#!/usr/bin/env bash
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate meta-thermo-env
set -euo pipefail
shopt -s nullglob
source config/bash/config.sh
echo "$(date +"%F %H:%M:%S")"
SAMPLES_TSV="${1}"
[[ -s "${SAMPLES_TSV}" ]] || { echo "ERROR: metadata $SAMPLES_TSV missing or empty" >&2; exit 1; }
declare -A R1 R2
SAMPLES=()
while IFS=$'\t' read -r sample_id \
    || [[ -n "${sample_id}" ]]
do
    [[ "${sample_id}" == "sample_id" || -z "${sample_id}" ]] && continue
    r1="${bowtie2_decontam_fastq_dir}/${sample_id}_trim_decontam_R1.fastq.gz"
    r2="${bowtie2_decontam_fastq_dir}/${sample_id}_trim_decontam_R2.fastq.gz"
    for f in "${r1}" "${r2}"; do
        [[ -s "${f}" ]] || { echo "ERROR: missing input for ${sample_id}: ${f}" >&2; exit 1; } # True if file exists and has a size greater than zero.
    done
    SAMPLES+=("${sample_id}"); R1[$sample_id]="${r1}"; R2[$sample_id]="${r2}"
done <  "${SAMPLES_TSV}"

(( ${#SAMPLES[@]} > 0 )) || { echo "ERROR: no samples found" >&2; exit 1; } # Counts the number of elements in the final array of samples. If
                                                    # the array is empty, raises an error

cd $metathermo_out_dir
for SAMPLE_ID in "${SAMPLES[@]}" 
do 
    echo "${SAMPLE_ID}"
    python "${meta_thermo_mpt_calculation_script}" \
        "${prokka_output_dir}"/"${SAMPLE_ID}"_prokka/"${SAMPLE_ID}"_pred_prokka.faa \
        "${SAMPLE_ID}"_MPT.csv "${SAMPLE_ID}"_AA.csv; 
done

cat *MPT.csv > all_samples_MPT.csv
conda activate mag_assembly-tools

Rscript "${meta_thermo_mpt_summary_script}" all_samples_MPT.csv\
    all_samples_MPT.csv

# Run on raw reads
for SAMPLE_ID in "${SAMPLES[@]}" 
do 
    echo "Running Meta-thermo on ${SAMPLE_ID}" unassembled reads
    # Read 1
    python "${meta_thermo_mpt_calculation_script}" \
        "${R1[$SAMPLE_ID]}" \
        "${SAMPLE_ID}"_raw_read_R1_MPT.csv "${SAMPLE_ID}"_raw_read_R1_AA.csv
    # Read 2
    python "${meta_thermo_mpt_calculation_script}" \
        "${R2[$SAMPLE_ID]}" \
        "${SAMPLE_ID}"_raw_read_R2_MPT.csv "${SAMPLE_ID}"_raw_read_R2_AA.csv;
done
