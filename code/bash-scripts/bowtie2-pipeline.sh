#!/usr/bin/env bash
nthreads=12
mem_req=8G
mem_req_sort=4G
nthreads_sort=8
# 1. Remove host DNA with bowtie2
### Activate qc-tools environment that has samtools and bowtie2
conda activate qc-tools
## 1.1 Bowtie2 mapping against host sequence
### RUN WITHOUT TRIMMOMATIC OR TRF
### FILE is the entire path (data/fastq/yasuda-fastq/2D10_wms_L3_1.fq.gz)
### SAMPLE is path and file without extension (data/fastq/yasuda-fastq/2D10_wms)
### base_name is just the file without extension or path (2D10_wms)
#mkdir output
#mkdir output/bowtie2_output_sam
#mkdir output/bowtie2_output_bam
#mkdir output/bowtie2_filtered_bam
#mkdir output/bowtie2_sorted_bam
#mkdir data/bowtie2_host_removed_fastq
fastq_dir=data/fastq/yasuda-fastq
bowtie2_output_sam_dir=output/bowtie2_output_sam
bowtie2_output_bam_dir=output/bowtie2_output_bam
bowtie2_filtered_bam_dir=output/bowtie2_filtered_bam
bowtie2_sorted_bam_dir=output/bowtie2_sorted_bam
bowtie2_host_removed_fastq_dir=data/bowtie2_host_removed_fastq
bowtie2_indices_dir=~/common_data/bowtie2_indices
for FILE in ${fastq_dir}/*L3_1.fq.gz
do 
  SAMPLE=$(echo ${FILE} | sed "s/_L3_1\.fq\.gz//")
  base_name=$(basename "$SAMPLE" )
  bowtie2 -p ${nthreads} -x ${bowtie2_indices_dir}/all_hosts_reference \
  -1 ${fastq_dir}/${base_name}_L3_1.fq.gz \
  -2 ${fastq_dir}/${base_name}_L3_2.fq.gz \
  -S ${bowtie2_output_sam_dir}/${base_name}_mapped_and_unmapped.sam;
done

## 1.2 Convert file .sam to .bam
for FILE in ${bowtie2_output_sam_dir}/*mapped_and_unmapped.sam
do 
  SAMPLE=$(echo ${FILE} | sed "s/_mapped_and_unmapped\.sam//")
  base_name=$(basename "$SAMPLE" )
  samtools view -bS ${bowtie2_output_sam_dir}/${base_name}_mapped_and_unmapped.sam > \
    ${bowtie2_output_bam_dir}/${base_name}_mapped_and_unmapped.bam;
done

## 1.3 Filter unmapped reads (unmapped to host genome)
### SAMtools SAM-flag filter: get unmapped pairs (both reads R1 and R2 unmapped)
for FILE in ${bowtie2_output_bam_dir}/*mapped_and_unmapped.bam 
do 
  SAMPLE=$(echo ${FILE} | sed "s/_mapped_and_unmapped\.bam//")
  base_name=$(basename "$SAMPLE" )
  samtools view -b -f 12 -F 256 \
   ${bowtie2_output_bam_dir}/${base_name}_mapped_and_unmapped.bam \
   > ${bowtie2_filtered_bam_dir}/${base_name}_bothReadsUnmapped.bam;
done
###-f  12    # Extract only ( -f ) alignments with both reads unmapped: <read unmapped><mate unmapped>
###-F 256    # Do not (  -F  ) extract alignments which are: <not primary alignment>

## 1.4 Split paired-end reads into separated fastq files .._R1 .._R2
### sort bam file by read name ( -n ) to have paired reads next to each other 
### (${nthreads} parallel threads, each using up to 5G memory)
for FILE in ${bowtie2_filtered_bam_dir}/*bothReadsUnmapped.bam 
do 
  SAMPLE=$(echo ${FILE} | sed "s/_bothReadsUnmapped\.bam//")
  base_name=$(basename "$SAMPLE" )
  samtools sort -n -m ${mem_req_sort} -@ ${nthreads_sort} ${bowtie2_filtered_bam_dir}/${base_name}_bothReadsUnmapped.bam \
    -o ${bowtie2_sorted_bam_dir}/${base_name}_bothReadsUnmapped_sorted.bam
  samtools fastq -@ ${nthreads_sort} ${bowtie2_sorted_bam_dir}/${base_name}_bothReadsUnmapped_sorted.bam \
  	-1 ${bowtie2_host_removed_fastq_dir}/${base_name}_host_removed_R1.fastq.gz \
  	-2 ${bowtie2_host_removed_fastq_dir}/${base_name}_host_removed_R2.fastq.gz \
  	-0 /dev/null -s /dev/null -n;
done