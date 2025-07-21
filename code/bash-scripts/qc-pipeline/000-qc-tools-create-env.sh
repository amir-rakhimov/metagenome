#!/usr/bin/env bash
#SBATCH -t 10:00:00
#SBATCH -J 20250516_16-36-qc-tools-create-env
#SBATCH --output jobreports/20250516_16-36-qc-tools-create-env-out.txt
#SBATCH --error jobreports/20250516_16-36-qc-tools-create-env-out.txt
# This script creates a conda environment for QC procedures and decontamination. 

# Note 1: The environment will also be used in the MAG assembly because it contains samtools.
# Note 2: All FASTQ files are gzip-compressed
source ~/miniconda3/etc/profile.d/conda.sh
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
source ~/miniconda3/etc/profile.d/conda.sh
# conda create -y --name qc-tools -c bioconda bowtie2 cutadapt samtools trf fastqc \
#   multiqc seqtk bbmap blast csvtk htseq
# conda activate qc-tools

# fastqc -h

# Install Java (not necessary for supercomputer)
curl https://download.oracle.com/java/24/latest/jdk-24_linux-x64_bin.tar.gz \
  --output ~/jdk-24_linux-x64_bin.tar.gz
tar zxvf ~/jdk-24_linux-x64_bin.tar.gz -C ~/
rm ~/jdk-24_linux-x64_bin.tar.gz 

wget https://github.com/broadinstitute/picard/releases/download/3.4.0/picard.jar \
 -O ~/picard.jar

wget https://github.com/tao-bioinfo/gff3sort/archive/master.zip \
 -O ~/gff3sort.zip
unzip ~/gff3sort.zip -d ~/

curl -L https://github.com/lh3/minimap2/releases/download/v2.29/minimap2-2.29_x64-linux.tar.bz2 \
 --output ~/minimap2-2.29_x64-linux.tar.bz2 
tar -jxvf ~/minimap2-2.29_x64-linux.tar.bz2  -C ~/