#!/usr/bin/env bash
# Metaquast is a metagenome assembly evaluation tools
date_var=$(date -I|sed 's/-//g')
time_var=$(date +%T |sed 's/:/_/g' )
date_time=${date_var}_${time_var}
megahit_date_time=20240607_12_52_52
nthreads=10
project_home_dir=~/projects/metagenome
mag_assembly_dir=output/mag_assembly
megahit_output_dir=output/mag_assembly/megahit_output
megahit_aligned_reads_dir=output/mag_assembly/megahit_output/alignedreads
metaquast_output_dir=output/mag_assembly/metaquast_output

mkdir output/mag_assembly/metaquast_output
# 1. Install QUAST
wget https://github.com/ablab/quast/releases/download/quast_5.2.0/quast-5.2.0.tar.gz -O ~/quast-5.2.0.tar.gz
tar -xzf ~/quast-5.2.0.tar.gz -C ~/

for FILE_DIR in ${megahit_output_dir}/*megahit_asm
do 
 SAMPLE=$(echo ${FILE_DIR} | sed "s/\.megahit_asm//" |sed "s/${megahit_date_time}_//")
 base_name=$(basename "$SAMPLE" )
 ~/quast-5.2.0/metaquast.py \
    ${megahit_output_dir}/${megahit_date_time}_${base_name}.megahit_asm/${megahit_date_time}_${base_name}_final.contigs.fa  \
    -o ${metaquast_output_dir}/${date_var}_metaquast_${base_name} \
    --threads ${nthreads} \
    2>&1 |tee ${metaquast_output_dir}/${date_var}_${base_name}_metaquast_report.txt;
done