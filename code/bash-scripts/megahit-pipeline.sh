#!/usr/bin/env bash
# MEGAHIT assembles reads into contigs
date_var=$(date -I|sed 's/-//g')
time_var=$(date +%T |sed 's/:/_/g' )
date_time=${date_var}_${time_var}
nthreads=10
nthreads_sort=8
mem_req=8G
mem_req_sort=4G
project_home_dir=~/projects/metagenome
bowtie2_decontam_fastq_dir=data/bowtie2_decontam_fastq
mag_assembly_dir=output/mag_assembly
megahit_output_dir=output/mag_assembly/megahit_output
megahit_aligned_reads_dir=output/mag_assembly/megahit_output/alignedreads
# 1. Install conda environment
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda create -yqn mag_assembly-tools -c conda-forge -c bioconda megahit \
    prodigal metabat2 checkm-genome gtdbtk
conda activate mag_assembly-tools

megahit --test  # run on a toy dataset
# cd ../..

cd ~/projects/metagenome
mkdir output/mag_assembly
mkdir output/mag_assembly/megahit_output
mkdir output/mag_assembly/megahit_output/alignedreads

# 2. Assemble
for FILE in ${bowtie2_decontam_fastq_dir}/*decontam_R1.fastq.gz
do 
 SAMPLE=$(echo ${FILE} | sed "s/_decontam_R1\.fastq\.gz//")
 base_name=$(basename "$SAMPLE" )
 megahit \
    -1 ${bowtie2_decontam_fastq_dir}/${base_name}_decontam_R1.fastq.gz \
    -2 ${bowtie2_decontam_fastq_dir}/${base_name}_decontam_R2.fastq.gz \
    -o ${megahit_output_dir}/${date_time}_${base_name}_decontam.megahit_asm \
    -t ${nthreads} 2>&1 |tee ${megahit_output_dir}/${date_time}_${base_name}_megahit_report.txt;
done

# This one actually has one contig
# megahit \
#     -1 ~/miniconda3/envs/mag_assembly-tools/share/megahit/test_data/r3_1.fa \
#     -2 ~/miniconda3/envs/mag_assembly-tools/share/megahit/test_data/r3_2.fa \
#     -o ${megahit_output_dir}/${date_time}_r3.megahit_asm \
#     -t ${nthreads} 2>&1 |tee ${megahit_output_dir}/${date_time}_r3_megahit_report.txt
# generates a directory 2D10_subset.megahit_asm 
# which contains checkpoints.txt, done, final.contigs.fa, intermediate_contig/, log, options.json

# --k-min : minimum kmer size (<= 255), must be odd number [default 21] 
# --k-max : maximum kmer size (<= 255), must be odd number [default 141]
# --k-step : increment of kmer size of each iteration (<= 28), must be even number [default 12]

# Tip 1. for ultra complex metagenomics data such as soil, a larger kmin, say 27,
#  is recommended to reduce the complexity of the de Bruijn graph. Quality trimming is also recommended
# Tip 2. for high-depth generic data, large --k-min (25 to 31) is recommended
# Tip 3. smaller --k-step, say 10, is more friendly to low-coverage datasets


# 3. Calculate contig coverage and extract unassembled reads
# 3.1 Use BBmap and samtools
conda activate qc-tools

# 3.2 Align reads with bbwrap.sh:
for FILE_DIR in ${megahit_output_dir}/*megahit_asm
do 
 SAMPLE=$(echo ${FILE_DIR} | sed "s/\.megahit_asm//" |sed "s/${date_time}_//")
 base_name=$(basename "$SAMPLE" )
 mv ${FILE_DIR}/final.contigs.fa \
  ${FILE_DIR}/${date_time}_${base_name}_final.contigs.fa
 bbwrap.sh ref=${FILE_DIR}/${date_time}_${base_name}_final.contigs.fa \
    in=${bowtie2_decontam_fastq_dir}/${base_name}_R1.fastq.gz \
    in2=${bowtie2_decontam_fastq_dir}/${base_name}_R2.fastq.gz \
    out=${megahit_aligned_reads_dir}/${date_time}_${base_name}_mapped_and_unmapped.sam \
    kfilter=22 \
    subfilter=15 \
    maxindel=80 2>&1 |tee ${megahit_output_dir}/${date_time}_${base_name}_bbwrap_report.txt
 gzip -9 --best ${megahit_aligned_reads_dir}/${date_time}_${base_name}_mapped_and_unmapped.sam;
done

## TEST commands: work
# mv ${megahit_output_dir}/${date_time}_r3.megahit_asm/final.contigs.fa \
#   ${megahit_output_dir}/${date_time}_r3.megahit_asm/${date_time}_r3_final.contigs.fa
# bbwrap.sh ref=${megahit_output_dir}/${date_time}_r3.megahit_asm/${date_time}_r3_final.contigs.fa \
#     in=~/miniconda3/envs/mag_assembly-tools/share/megahit/test_data/r3_1.fa \
#     in2=~/miniconda3/envs/mag_assembly-tools/share/megahit/test_data/r3_2.fa \
#     out=${megahit_aligned_reads_dir}/${date_time}_r3_mapped_and_unmapped.sam.gz \
#     kfilter=22 \
#     subfilter=15 \
#     maxindel=80 2>&1 |tee ${megahit_output_dir}/${date_time}_r3_bbwrap_report.txt
###

# 3.3 Output per contig coverage to cov.txt with pileup.sh:
# BBMap generates coverage information by internally using Pileup
for FILE in ${megahit_aligned_reads_dir}/*_mapped_and_unmapped.sam.gz
do 
 SAMPLE=$(echo ${FILE} | sed "s/_mapped_and_unmapped\.sam\.gz//" |sed "s/${date_time}_//")
 base_name=$(basename "$SAMPLE" )
 pileup.sh in=${megahit_aligned_reads_dir}/${date_time}_${base_name}_mapped_and_unmapped.sam.gz \
    out=${megahit_aligned_reads_dir}/${date_time}_${base_name}_cov.txt \
    2>&1 |tee ${megahit_output_dir}/${date_time}_${base_name}_pileup_report.txt;
done

# 3.4 Extract unmapped reads (SE to unmapped.se.fq and PE to unmapped.pe.fq):
for FILE in ${megahit_aligned_reads_dir}/*_mapped_and_unmapped.sam.gz
do 
 SAMPLE=$(echo ${FILE} | sed "s/_mapped_and_unmapped\.sam\.gz//" |sed "s/${date_time}_//")
 base_name=$(basename "$SAMPLE" )
 samtools view -u -f4 ${megahit_aligned_reads_dir}/${date_time}_${base_name}_mapped_and_unmapped.sam.gz | \
    samtools bam2fq -s ${megahit_aligned_reads_dir}/${date_time}_${base_name}_unmapped.se.fq - > \
    ${megahit_aligned_reads_dir}/${date_time}_${base_name}_unmapped.pe.fq 
 gzip -9 --best ${megahit_aligned_reads_dir}/${date_time}_${base_name}_unmapped.se.fq
 gzip -9 --best ${megahit_aligned_reads_dir}/${date_time}_${base_name}_unmapped.pe.fq;
done
# samtools view – views and converts SAM/BAM/CRAM files
# -u, --uncompressed : Output uncompressed data
# -f4 flag filters out reads that are unmapped.

# bam2fq converts a bam into FASTQ format
# -s FILE : Write singleton reads in FASTQ format to FILE instead of outputting them.

# The final output is then directed to the standard output (stdout) using the - symbol 
#   followed by the redirection operator >
# The - > combination redirects the output of the second command (samtools bam2fq)
#    to a file named output/unmapped.pe.fq. The hyphen (-) after the redirection operator
#    indicates that the output should be sent to stdout (standard output), 
#    which can then be captured or saved to a file as needed.


# 3.5 Convert file .sam to .bam
for FILE in ${megahit_aligned_reads_dir}/*_mapped_and_unmapped.sam.gz
do 
  SAMPLE=$(echo ${FILE} | sed "s/_mapped_and_unmapped\.sam\.gz//" |sed "s/${date_time}_//")
  base_name=$(basename "$SAMPLE" )
  samtools view -b ${megahit_aligned_reads_dir}/${date_time}_${base_name}_mapped_and_unmapped.sam.gz > \
    ${megahit_aligned_reads_dir}/${date_time}_${base_name}_mapped_and_unmapped.bam;
done
# samtools view – views and converts SAM/BAM/CRAM files
# -b, --bam :Output in the BAM format.

## 4. Keep mapped reads (both mapped to contigs)
# 4.1 To get only the mapped reads use the parameter F, which works like -v of grep and skips the alignments for a specific flag.
for FILE in ${megahit_aligned_reads_dir}/*_mapped_and_unmapped.sam.gz
do 
  SAMPLE=$(echo ${FILE} | sed "s/_mapped_and_unmapped\.sam\.gz//" |sed "s/${date_time}_//")
  base_name=$(basename "$SAMPLE" )
  samtools view -b -F 12  \
   ${megahit_aligned_reads_dir}/${date_time}_${base_name}_mapped_and_unmapped.bam \
   > ${megahit_aligned_reads_dir}/${date_time}_${base_name}_bothReadsMapped.bam;
done

### 4.2 Sort bam file by coordinate (no options used)
for FILE in ${megahit_aligned_reads_dir}/*_bothReadsMapped.bam
do 
  SAMPLE=$(echo ${FILE} | sed "s/_bothReadsMapped\.bam//" |sed "s/${date_time}_//")
  base_name=$(basename "$SAMPLE" )
  samtools sort -m ${mem_req_sort} -@ ${nthreads_sort} \
  ${megahit_aligned_reads_dir}/${date_time}_${base_name}_bothReadsMapped.bam \
    -o ${megahit_aligned_reads_dir}/${date_time}_${base_name}_bothReadsMapped_sorted.bam;
done

## 4.3 Filter unmapped reads (unmapped to host genome)
### SAMtools SAM-flag filter: get unmapped pairs (both reads R1 and R2 unmapped)
for FILE in ${megahit_aligned_reads_dir}/*mapped_and_unmapped.bam 
do 
  SAMPLE=$(echo ${FILE} | sed "s/_mapped_and_unmapped\.bam//" |sed "s/${date_time}_//")
  base_name=$(basename "$SAMPLE" )
  samtools view -b -f 12 -F 256 \
   ${megahit_aligned_reads_dir}/${date_time}_${base_name}_mapped_and_unmapped.bam \
   > ${megahit_aligned_reads_dir}/${date_time}_${base_name}_bothReadsUnmapped.bam;
done
###-f  12    # Extract only ( -f ) alignments with both reads unmapped: <read unmapped><mate unmapped>
###-F 256    # Do not (  -F  ) extract alignments which are: <not primary alignment>

## 4.4 Split paired-end reads into separated fastq files .._R1 .._R2
### sort bam file by read name ( -n ) to have paired reads next to each other 
### (${nthreads} parallel threads, each using up to 5G memory)
for FILE in ${megahit_aligned_reads_dir}/*bothReadsUnmapped.bam 
do 
  SAMPLE=$(echo ${FILE} | sed "s/_bothReadsUnmapped\.bam//" |sed "s/${date_time}_//")
  base_name=$(basename "$SAMPLE" )
  samtools sort -n -m ${mem_req_sort} -@ ${nthreads_sort} \
  ${megahit_aligned_reads_dir}/${date_time}_${base_name}_bothReadsUnmapped.bam \
    -o ${megahit_aligned_reads_dir}/${date_time}_${base_name}_bothReadsUnmapped_sorted.bam
  samtools fastq -@ ${nthreads_sort} \
    ${megahit_aligned_reads_dir}/${date_time}_${base_name}_bothReadsUnmapped_sorted.bam \
  	-1 ${megahit_aligned_reads_dir}/${date_time}_${base_name}_decontam_R1.fastq.gz \
  	-2 ${megahit_aligned_reads_dir}/${date_time}_${base_name}_decontam_R2.fastq.gz \
  	-0 /dev/null -s /dev/null -n;
done

# samtools sort: Sort alignments by leftmost coordinates, 
#  by read name when -n or -N are used,
#  by tag contents with -t,
#  or a minimiser-based collation order with -M.
# -m : Approximately the maximum required memory per thread
# samtools fastq: Converts a BAM or CRAM into either FASTQ or FASTA format
#  depending on the command invoked