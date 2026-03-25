#!/usr/bin/env bash
prokka_date_time=20250626_22_11_43
project_home_dir=~/projects/metagenome
prokka_output_dir="${project_home_dir}"/output/mag_assembly/prokka_output
bowtie2_decontam_fastq_dir="${project_home_dir}"/data/bowtie2_decontam_fastq
metathermo_out_dir="${project_home_dir}"/output/mag_assembly/meta-thermo_output
conda activate meta-thermo-env

cd $metathermo_out_dir
for FILE_DIR in "${prokka_output_dir}"/"${prokka_date_time}"*_prokka
do 
    SAMPLE=$(echo "${FILE_DIR}" | sed "s/_prokka//" |sed "s/${prokka_date_time}_//")
    base_name=$(basename "$SAMPLE" )
    echo "${base_name}"
    python ~/meta-thermo/scripts/MPTcalculation.py \
        "${FILE_DIR}"/"${prokka_date_time}"_"${base_name}"_pred_prokka.faa \
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
