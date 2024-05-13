#!/usr/bin/env bash
source ~/miniconda3/etc/profile.d/conda.sh
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels biobakery
source ~/miniconda3/etc/profile.d/conda.sh
#conda create -y --name qc-tools -c bioconda bowtie2 cutadapt samtools trf fastqc multiqc seqtk
conda activate qc-tools
conda install -y cutadapt
