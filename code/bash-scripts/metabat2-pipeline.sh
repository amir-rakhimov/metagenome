#!/usr/bin/env bash
# MetaBAT2 is a binning software: contigs into genomes
date_var=$(date -I|sed 's/-//g')
time_var=$(date +%T |sed 's/:/_/g' )
date_time=${date_var}_${time_var}
megahit_date_time=20240607_12_52_52
nthreads=10
project_home_dir=~/projects/metagenome
bowtie2_decontam_fastq_dir=data/bowtie2_decontam_fastq
mag_assembly_dir=output/mag_assembly
megahit_output_dir=output/mag_assembly/megahit_output
megahit_aligned_reads_dir=output/mag_assembly/megahit_output/alignedreads
bam_contig_depths_dir=output/mag_assembly/bam_contig_depths
metabat2_output_dir=output/mag_assembly/metabat2_output
checkm_output_dir=output/mag_assembly/checkm_output
gtdbtk_output_dir=output/mag_assembly/gtdbtk_output
# Set conda channel priority
# conda config --add channels defaults
# conda config --add channels bioconda
# conda config --add channels conda-forge
# conda create -yqn mag_assembly-tools -c conda-forge -c bioconda megahit \
#     prodigal metabat2 checkm-genome gtdbtk
conda activate mag_assembly-tools

mkdir output/mag_assembly/metabat2_output
mkdir output/mag_assembly/checkm_output
mkdir output/mag_assembly/gtdbtk_output
mkdir output/mag_assembly/bam_contig_depths

# MetaBAT 2 USAGE: running on command line
# Be careful to have bams sorted first!
# 1. Create a depth file
for FILE_DIR in ${megahit_output_dir}/*megahit_asm
do 
 SAMPLE=$(echo ${FILE_DIR} | sed "s/\.megahit_asm//" |sed "s/${megahit_date_time}_//")
 base_name=$(basename "$SAMPLE" )
 jgi_summarize_bam_contig_depths --percentIdentity 97 \
    --minContigLength 1000 \
    --minContigDepth 1.0  \
    --referenceFasta ${megahit_output_dir}/${megahit_date_time}_${base_name}.megahit_asm/${megahit_date_time}_${base_name}_final.contigs.fa \
    ${megahit_aligned_reads_dir}/${megahit_date_time}_${base_name}_bothReadsMapped_sorted.bam \
    --outputDepth ${bam_contig_depths_dir}/${date_time}_${base_name}.depth.txt;
done
# runMetaBat.sh <options> assembly.fasta sample1.bam [sample2.bam ...]
# runMetaBat.sh ../megahit_output/2D10_trim_decontam.megahit_asm/2D10_trim_decontam_final.contigs.fa     \
#     ../megahit_output/alignedreads/2D10_aln_bothReadsMapped_sorted.bam

# 2. Run MetaBAT2 and CheckM
for FILE_DIR in ${megahit_output_dir}/*megahit_asm
do 
 SAMPLE=$(echo ${FILE_DIR} | sed "s/\.megahit_asm//" |sed "s/${megahit_date_time}_//")
 base_name=$(basename "$SAMPLE" )
 metabat2 -i ${megahit_output_dir}/${megahit_date_time}_${base_name}.megahit_asm/${megahit_date_time}_${base_name}_final.contigs.fa \
    -a ${bam_contig_depths_dir}/${date_time}_${base_name}.depth.txt  \
    -o ${metabat2_output_dir}/${date_time}_${base_name}/${date_time}_${base_name}_bin \
    -v;
done

# metabat2 -i ${megahit_output_dir}/demo_assembly.fa.gz \
#     -a ${bam_contig_depths_dir}/demo.depth.txt  \
#     -o ${metabat2_output_dir}/${date_time}_demo/${date_time}_demo_bin \
#     -v
# -o means the output directory and start of the file name (here, bin.1.fa, bin.2.fa, etc)

for FILE_DIR in ${metabat2_output_dir}/${date_time}*
do 
 SAMPLE=$(echo ${FILE_DIR} | sed "s/${date_time}_//")
 base_name=$(basename "$SAMPLE" )
 checkm lineage_wf -t ${nthreads} \
    -x fa \
    ${metabat2_output_dir}/${date_time}_${base_name} \
    ${checkm_output_dir}/${date_time}_${base_name}_SCG \
    -f ${checkm_output_dir}/${date_time}_${base_name}_CheckM.txt;
done
# -f, --file FILE       print results to file (default: stdout)
# -t, --threads THREADS
# -x, --extension EXTENSION: extension of bins (other files in directory are ignored) (default: fna)


export GTDBTK_DATA_PATH=~/miniconda3/envs/mag_assembly-tools/share/gtdbtk-1.7.0/db/

for FILE_DIR in ${metabat2_output_dir}/${date_time}*
do 
 SAMPLE=$(echo ${FILE_DIR} | sed "s/${date_time}_//")
 base_name=$(basename "$SAMPLE" )
 gtdbtk classify_wf --genome_dir ${metabat2_output_dir}/${date_time}_${base_name}/ \
    --out_dir ${gtdbtk_output_dir}/${date_time}_${base_name}_GTDBtk \
    -x fa \
    --cpus {config[cores][gtdbtk]};
done