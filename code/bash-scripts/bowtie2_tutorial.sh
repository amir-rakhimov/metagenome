#!/usr/bin/env bash
mkdir bowtie2_protocol
conda activate qc-tools
# Indexing reference genome
## Use lambda phage genome
mkdir reference
mkdir index
bowtie2-build reference/lambda_virus.fa index/lambda_virus

# 1. Aligning example reads
bowtie2 -x index/lambda_virus -U reads/reads_1.fq -S eg1.sam

# 2. Paired-end example: keep both aligned and unaligned reads
bowtie2 -x index/lambda_virus -1 reads/reads_1.fq -2 reads/reads_2.fq -S eg2.sam

bowtie2 -p 8 -x index/lambda_virus \
  -1 reads/reads_1.fq \
  -2 reads/reads_2.fq \
  --very-sensitive-local \
  --un-conc-gz \
  SAMPLE_host_removed \
  > SAMPLE_mapped_and_unmapped.sam
  
# SAMPLE_host_removed.1 and SAMPLE_host_removed.2 are gz files without gz ending
mv SAMPLE_host_removed.1 SAMPLE_host_removed_R1.fastq.gz
mv SAMPLE_host_removed.2 SAMPLE_host_removed_R2.fastq.gz

mkdir output

# Option --un-conc shows results like samtools options -F 2 
# (excluding reads "mapped in proper pair")
# 3. Bowtie with SAMtools
# 3.1 Bowtie2 mapping against host sequence
bowtie2 -x index/lambda_virus \
  -1 reads/reads_1.fq \
  -2 reads/reads_2.fq \
  -S output/SAMPLE_mapped_and_unmapped.sam

bowtie2 -x index/lambda_virus \
  -1 reads/reads_1.fq.gz \
  -2 reads/reads_2.fq.gz \
  -S output/SAMPLE_mapped_and_unmapped.sam

# 3.2 Convert file .sam to .bam
samtools view -bS output/SAMPLE_mapped_and_unmapped.sam > output/SAMPLE_mapped_and_unmapped.bam

# 3.3 Filter unmapped reads (unmapped to host genome)
## SAMtools SAM-flag filter: get unmapped pairs (both reads R1 and R2 unmapped)
samtools view -b -f 12 -F 256 \
   output/SAMPLE_mapped_and_unmapped.bam \
   > output/SAMPLE_bothReadsUnmapped.bam
###-f  12    # Extract only ( -f ) alignments with both reads unmapped: <read unmapped><mate unmapped>
###-F 256    # Do not (  -F  ) extract alignments which are: <not primary alignment>

# 3.4 Split paired-end reads into separated fastq files .._R1 .._R2
## sort bam file by read name ( -n ) to have paired reads next to each other 
## (2 parallel threads, each using up to 5G memory)
samtools sort -n -m 5G -@ 2 output/SAMPLE_bothReadsUnmapped.bam \
  -o output/SAMPLE_bothReadsUnmapped_sorted.bam
samtools fastq -@ 8 output/SAMPLE_bothReadsUnmapped_sorted.bam \
  -1 output/SAMPLE_host_removed_R1.fastq.gz \
  -2 output/SAMPLE_host_removed_R2.fastq.gz \
  -0 /dev/null -s /dev/null -n

