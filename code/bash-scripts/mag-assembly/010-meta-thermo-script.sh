#!/usr/bin/env bash
source ~/miniconda3/etc/profile.d/conda.sh
conda activate meta-thermo-env
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"

cd $metathermo_out_dir
for FILE_DIR in "${prokka_output_dir}"/*_prokka
do 
    SAMPLE=$(echo "${FILE_DIR}" | sed "s/_prokka//")
    base_name=$(basename "$SAMPLE" )
    echo "${base_name}"
    python ~/meta-thermo/scripts/MPTcalculation.py \
        "${FILE_DIR}"/"${base_name}"_pred_prokka.faa \
        "${base_name}"_MPT.csv "${base_name}"_AA.csv; 
done

cat *MPT.csv > all_samples_MPT.csv
conda activate mag_assembly-tools

Rscript ~/projects/metagenome/code/r-scripts/mag-assembly/mpt-summary.R all_samples_MPT.csv\
    all_samples_MPT.csv

# Run on raw reads
for FILE in "${bowtie2_decontam_fastq_dir}"/*decontam_R1.fastq.gz
do 
    SAMPLE=$(echo "${FILE}" | sed "s/_trim_decontam_R1\.fastq\.gz//")
    base_name=$(basename "$SAMPLE" )
    echo "Running Meta-thermo on ${base_name}" unassembled reads
    # Read 1
    python ~/meta-thermo/scripts/MPTcalculation.py \
        "${bowtie2_decontam_fastq_dir}"/"${base_name}"_trim_decontam_R1.fastq.gz \
        "${base_name}"_raw_read_R1_MPT.csv "${base_name}"_raw_read_R1_AA.csv
    # Read 2
    python ~/meta-thermo/scripts/MPTcalculation.py \
        "${bowtie2_decontam_fastq_dir}"/"${base_name}"_trim_decontam_R2.fastq.gz \
        "${base_name}"_raw_read_R2_MPT.csv "${base_name}"_raw_read_R2_AA.csv;
done
