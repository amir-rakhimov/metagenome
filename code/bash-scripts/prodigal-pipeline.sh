#!/usr/bin/env bash
# Prodigal is a gene prediction software
date_var=$(date -I|sed 's/-//g')
time_var=$(date +%T |sed 's/:/_/g' )
date_time=${date_var}_${time_var}
megahit_date_time=20240607_12_52_52
nthreads=10
project_home_dir=~/projects/metagenome
mag_assembly_dir=output/mag_assembly
megahit_output_dir=output/mag_assembly/megahit_output
megahit_aligned_reads_dir=output/mag_assembly/megahit_output/alignedreads
prodigal_output_dir=output/mag_assembly/prodigal_output

mkdir output/mag_assembly/prodigal_output
# Install Prodigal
# conda config --add channels defaults
# conda config --add channels bioconda
# conda config --add channels conda-forge

# conda create -yqn mag_assembly-tools -c conda-forge -c bioconda megahit \
#     prodigal metabat2 checkm-genome gtdbtk
conda activate mag_assembly-tools

# Anonymous mode (for metagenomes)
# should I use -p meta or -p anon?
for FILE_DIR in ${megahit_output_dir}/*megahit_asm
do 
 SAMPLE=$(echo ${FILE_DIR} | sed "s/\.megahit_asm//" |sed "s/${megahit_date_time}_//")
 base_name=$(basename "$SAMPLE" )
 mkdir ${prodigal_output_dir}/${date_time}_${base_name}_prodigal
 prodigal \
    -i ${megahit_output_dir}/${megahit_date_time}_${base_name}.megahit_asm/${megahit_date_time}_${base_name}_final.contigs.fa \
    -o ${prodigal_output_dir}/${date_time}_${base_name}_prodigal/${date_time}_${base_name}_prodigal_coords.gbk \
    -a ${prodigal_output_dir}/${date_time}_${base_name}_prodigal/${date_time}_${base_name}_prodigal_proteins.faa \
    -p meta 2>&1 |tee ${prodigal_output_dir}/${date_time}_${base_name}_prodigal_report.txt;
done