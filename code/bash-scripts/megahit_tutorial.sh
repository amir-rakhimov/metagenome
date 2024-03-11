#!/usr/bin/env bash
cd ~
# 1. Download MEGAHIT
wget https://github.com/voutcn/megahit/releases/download/v1.2.9/MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz
tar zvxf MEGAHIT-1.2.9-Linux-x86_64-static.tar.gz
cd MEGAHIT-1.2.9-Linux-x86_64-static/bin/
./megahit --test  # run on a toy dataset
cd ../..

cd projects/megahit_protocol
mkdir data
mkdir data/fastq
mkdir data/fastq/raw
mkdir output
mkdir output/assembly
mkdir output/alignedreads

# 2. Download the data
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR341/SRR341725/SRR341725_[12].fastq.gz 
mv SRR341725_*.fastq.gz data/fastq/raw/.

# 3. Assemble
~/MEGAHIT-1.2.9-Linux-x86_64-static/bin/megahit \
    -1 data/fastq/raw/SRR341725_1.fastq.gz \
    -2 data/fastq/raw/SRR341725_2.fastq.gz \
    -o output/assembly/SRR341725.megahit_asm \
    -t ${nthreads}
# --k-min : minimum kmer size (<= 255), must be odd number [default 21] 
# --k-max : maximum kmer size (<= 255), must be odd number [default 141]
# --k-step : increment of kmer size of each iteration (<= 28), must be even number [default 12]

# Tip 1. for ultra complex metagenomics data such as soil, a larger kmin, say 27,
#  is recommended to reduce the complexity of the de Bruijn graph. Quality trimming is also recommended
# Tip 2. for high-depth generic data, large --k-min (25 to 31) is recommended
# Tip 3. smaller --k-step, say 10, is more friendly to low-coverage datasets


# 4. Calculate contig coverage and extract unassembled reads
# 4.1 Use BBmap and samtools
conda activate qc-tools

# 4.2 Align reads with bbwrap.sh:
bbwrap.sh ref=output/assembly/SRR341725.megahit_asm/final.contigs.fa \
    in=data/fastq/raw/SRR341725_1.fastq.gz \
    in2=data/fastq/raw/SRR341725_2.fastq.gz \
    out=output/alignedreads/aln.sam.gz \
    kfilter=22 \
    subfilter=15 \
    maxindel=80

# 4.3 Output per contig coverage to cov.txt with pileup.sh:
pileup.sh in=output/alignedreads/aln.sam.gz out=output/cov.txt

# 4.4 Extract unmapped reads (SE to unmapped.se.fq and PE to unmapped.pe.fq):
samtools view -u -f4 output/alignedreads/aln.sam.gz | \
    samtools bam2fq -s output/unmapped.se.fq - > output/unmapped.pe.fq
