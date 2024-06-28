#!/usr/bin/env bash
nthreads=12
mem_req=8G
mem_req_sort=4G
nthreads_sort=8
project_home_dir=~/projects/metagenome
fastq_dir=data/fastq/yasuda-fastq
reference_genomes_dir=~/common_data/reference_genomes
bowtie2_indices_dir=~/common_data/bowtie2_indices
fastqc_output_dir=output/qc_pipeline/fastqc_output
multiqc_output_dir=output/qc_pipeline/multiqc_output
cutadapt_output_dir=output/qc_pipeline/cutadapt_output
temp_fasta_dir=output/qc_pipeline/temp_fasta
trf_output_dir=output/qc_pipeline/trf_output
fastq_norepeats_dir=output/qc_pipeline/fastq_norepeats
FWD_ADAPTER=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
REV_ADAPTER=GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG

# mkdir ~/common_data
# mkdir ~/common_data/reference_genomes
#mkdir output
# mkdir output/qc_pipeline
#mkdir output/qc_pipeline/fastqc_output
#mkdir output/qc_pipeline/multiqc_output
#mkdir output/qc_pipeline/cutadapt_output
#mkdir output/qc_pipeline/temp_fasta
#mkdir output/qc_pipeline/trf_output
#mkdir output/qc_pipeline/fastq_norepeats
#mkdir output/qc_pipeline/fastqc_output_decontam
#mkdir output/bowtie2_pipeline
#mkdir output/bowtie2_pipeline/bowtie2_output_sam
#mkdir output/bowtie2_pipeline/bowtie2_output_bam
#mkdir output/bowtie2_pipeline/bowtie2_filtered_bam
#mkdir output/bowtie2_pipeline/bowtie2_sorted_bam
#mkdir data/bowtie2_decontam_fastq
bowtie2_output_sam_dir=output/bowtie2_pipeline/bowtie2_output_sam
bowtie2_output_bam_dir=output/bowtie2_pipeline/bowtie2_output_bam
bowtie2_filtered_bam_dir=output/bowtie2_pipeline/bowtie2_filtered_bam
bowtie2_sorted_bam_dir=output/bowtie2_pipeline/bowtie2_sorted_bam
bowtie2_decontam_fastq_dir=data/bowtie2_decontam_fastq
fastqc_output_decontam_dir=output/qc_pipeline/fastqc_output_decontam
source ~/miniconda3/etc/profile.d/conda.sh
# Set conda channel priority
# conda config --add channels defaults
# conda config --add channels bioconda
# conda config --add channels conda-forge
# conda config --add channels biobakery
# # Create conda environment
# conda create -y --name qc-tools -c bioconda bowtie2 cutadapt samtools trf fastqc multiqc seqtk
conda activate qc-tools 
# Rename files
# for FILE in nmrF*.fq.gz; do mv "$FILE" $(echo "$FILE" | sed 's/nmrF_//'); done
# for FILE in *DKDN230040401-1A_HC3LYDSX7*.fq.gz; 
#   do mv "$FILE" $(echo "$FILE" | sed 's/DKDN230040401-1A_HC3LYDSX7/wms/'); done
# for FILE in *DKDN230040401-1A_HC557DSX7*.fq.gz; 
#   do mv "$FILE" $(echo "$FILE" | sed 's/DKDN230040401-1A_HC557DSX7/wms/'); done
# for FILE in *DKDN230040402-1A_HC557DSX7*.fq.gz; 
#   do mv "$FILE" $(echo "$FILE" | sed 's/DKDN230040402-1A_HC557DSX7/wms/'); done
# for FILE in *DKDN230040410-1A_HC52YDSX7*.fq.gz; 
#   do mv "$FILE" $(echo "$FILE" | sed 's/DKDN230040410-1A_HC52YDSX7/wms/'); done
# for FILE in *DKDN230040409-1A_HC3LYDSX7*.fq.gz; 
#   do mv "$FILE" $(echo "$FILE" | sed 's/DKDN230040409-1A_HC3LYDSX7/wms/'); done
# for FILE in *DKDN230040409-1A_HC557DSX7*.fq.gz; 
#   do mv "$FILE" $(echo "$FILE" | sed 's/DKDN230040409-1A_HC557DSX7/wms/'); done
# for FILE in *DKDN230040411-1A_HC52YDSX7*.fq.gz; 
#   do mv "$FILE" $(echo "$FILE" | sed 's/DKDN230040411-1A_HC52YDSX7/wms/'); done
# for FILE in *DKDN230040407-1A_HC557DSX7*.fq.gz; 
#   do mv "$FILE" $(echo "$FILE" | sed 's/DKDN230040407-1A_HC557DSX7/wms/'); done
# for FILE in *DKDN230040408-1A_HC557DSX7*.fq.gz; 
#   do mv "$FILE" $(echo "$FILE" | sed 's/DKDN230040408-1A_HC557DSX7/wms/'); done
# for FILE in *DKDN230040406-1A_HC557DSX7*.fq.gz; 
#   do mv "$FILE" $(echo "$FILE" | sed 's/DKDN230040406-1A_HC557DSX7/wms/'); done
# for FILE in *DKDN230040403-1A_HC557DSX7*.fq.gz; 
#   do mv "$FILE" $(echo "$FILE" | sed 's/DKDN230040403-1A_HC557DSX7/wms/'); done
# for FILE in *DKDN230040404-1A_HC557DSX7*.fq.gz; 
#   do mv "$FILE" $(echo "$FILE" | sed 's/DKDN230040404-1A_HC557DSX7/wms/'); done
# for FILE in *DKDN230040405-1A_HC3LYDSX7*.fq.gz; 
#   do mv "$FILE" $(echo "$FILE" | sed 's/DKDN230040405-1A_HC3LYDSX7/wms/'); done
# for FILE in *DKDN230040405-1A_HC557DSX7*.fq.gz; 
#   do mv "$FILE" $(echo "$FILE" | sed 's/DKDN230040405-1A_HC557DSX7/wms/'); done

# Download reference genomes
### Heter_glaber.v1.7_hic_pac
# curl https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCA_014060925.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT \
#   --output ${reference_genomes_dir}/Heter_glaber.v1.7_hic_pac.zip
### There is also a new assembly that is used by Ensembl: Naked_mole-rat_maternal
# curl https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCA_944319715.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT \
#   --output ${reference_genomes_dir}/Naked_mole-rat_maternal.zip
### DMR_v1.0_HiC Damaraland mole-rat
# curl https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_012274545.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT \
#   --output ${reference_genomes_dir}/DMR_v1.0_HiC.zip
# ### Human T2T
# curl https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_009914755.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT \
#   --output ${reference_genomes_dir}/T2T.zip
# ### Mouse GRCm39
# curl https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000001635.27/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT \
#   --output ${reference_genomes_dir}/GRCm39.zip
# ### human decoy genome
# wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5ss.fa.gz
# mv hs37d5ss.fa human_decoy.fa
# unzip ${reference_genomes_dir}/*.zip
# mv ncbi_dataset/data/GCF_012274545.1/GCF_012274545.1_DMR_v1.0_HiC_genomic.fna DMR_v1.0_HiC_genomic.fna
# # Merge them all into one FASTA file
# cat ${reference_genomes_dir}/*genomic.fna > ${reference_genomes_dir}/all_hosts_reference.fasta 
# cat ${reference_genomes_dir}/phiX174.fasta >>${reference_genomes_dir}/all_hosts_reference.fasta
# # Build a bowtie2 large index
# bowtie2-build --large-index --threads ${nthreads} ${reference_genomes_dir}/all_hosts_reference.fasta \
#   ${bowtie2_indices_dir}/all_hosts_reference

# The QC workflow is this:
### 1. Run FastQC and trim overrepresented sequences
### 1.1 Run MultiQC
### 1.2 Trim overrepresented sequences using Fastqc reports
### 2. Run Cutadapt
### Three filtering steps:
### (1) Remove reads containing adapters.
### (2) Remove reads containing N > 10% (N represents the base cannot be determined).
### (3) Remove reads containing low quality (Qscore<= 5) base which is over 50% of the total base.
### 3. Remove repeats with TRF
### 3.1 Create temporary fasta files with seqtk seq
### 3.2 Run TRF in parallel
### 3.3 Remove repeats identified by TRF from FASTQ files in parallel
### 4 Bowtie2
### 5. Run FastQC and MultiQC on decontaminated data as a final check

# 1. Run FastQC and trim overrepresented sequences
# fastqc ${fastq_dir}/*.fq.gz --outdir ${fastqc_output_dir}
## 1.1 Run MultiQC
# multiqc ${fastqc_output_dir}/ --outdir ${multiqc_output_dir}

## 1.2 Trim overrepresented sequences using Fastqc reports
#############
# TODO
#############################################################
# 2. Run Cutadapt
### Three filtering steps:
### (1) Remove reads containing adapters.
### (2) Remove reads containing N > 10% (N represents the base cannot be determined).
### (3) Remove reads containing low quality (Qscore<= 5) base which is over 50% of the total base.

for FILE in ${fastq_dir}/*wms_L3_1.fq.gz
do 
 SAMPLE=$(echo ${FILE} | sed "s/_wms_L3_1\.fq\.gz//")
 base_name=$(basename "$SAMPLE" )
 cutadapt \
  --cores=${nthreads} \
  -g ${FWD_ADAPTER} -G ${REV_ADAPTER} \
  --max-n 0.1 \
  -q 5 \
  -O 5 \
  --minimum-length 75 \
  --discard-trimmed \
  -o ${cutadapt_output_dir}/${base_name}_R1.trimmed.fastq.gz \
  -p ${cutadapt_output_dir}/${base_name}_R2.trimmed.fastq.gz \
  ${fastq_dir}/${base_name}_wms_L3_1.fq.gz ${fastq_dir}/${base_name}_wms_L3_2.fq.gz > \
  ${cutadapt_output_dir}/${base_name}_report.txt;
done

### FILE: data/fastq/yasuda-fastq/2D10_wms_L3_1.fq.gz
### SAMPLE: data/fastq/yasuda-fastq/2D10
### base_name: 2D10 
# cutadapt \
#   --cores=${nthreads} \                         ## number of cores for parallelisation
#   -g ${FWD_ADAPTER} -G ${REV_ADAPTER} \         ## 5' adapter
#   --max-n 0.1 \                                 ## do not allow > 10% Ns
#   -q 5 \                                        ## remove bases with basequal < 5
#   --minimum-length (readlength/2) \             ## minimum length after trimming
#   --discard-trimmed \                           ## Discard reads in which an adapter was found
#   -o read1_trimmed.fq -p read2_trimmed.fq \     ## outputs
#   read1.fq read2.fq                             ## inputs

# 3. Remove repeats with TRF
## 3.1 Create temporary fasta files with seqtk seq
# for FILE in ${cutadapt_output_dir}/*R1.trimmed.fastq.gz
# do 
#   SAMPLE=$(echo ${FILE} | sed "s/_R1\.trimmed\.fastq\.gz//")
#   base_name=$(basename "$SAMPLE" )
#   seqtk seq -a ${cutadapt_output_dir}/${base_name}_R1.trimmed.fastq.gz \
#     > ${temp_fasta_dir}/${base_name}_R1.trimmed.fasta
#   seqtk seq -a ${cutadapt_output_dir}/${base_name}_R2.trimmed.fastq.gz \
#     > ${temp_fasta_dir}/${base_name}_R2.trimmed.fasta;
# done

## 3.2 Run TRF in parallel
# cd ${trf_output_dir}
# # Find all files in the temp_fasta directory
# temp_fasta_files=$(find ${project_home_dir}/${temp_fasta_dir} -type f \( -name "*R1.trimmed.fasta"\
#  -o -name "*R2.trimmed.fasta" \))
# # run TRF parallel
# python3 ${project_home_dir}/code/python-scripts/parallel_trf.py \
#   --input-file-dir ${project_home_dir}/${temp_fasta_dir} \
#   --input-file-suffix "trimmed.fasta" \
#   --trf-output-dir . \
#   --trf-output-file-suffix "trf_out.dat" \
#   --nthreads ${nthreads} \
#   --files ${temp_fasta_files} 

# parallel_trf.py usage: 
#           python3 parallel_trf.py <input_file_dir> <input_file_suffix> \n"
#               "\t <trf_output_dir> <trf_output_file_suffix> <num_processes> \n"
#               "\t<file1> <file2> ...
# TRF command: 
# trf ${project_home_dir}/${temp_fasta_dir}/${base_name}_R1.trimmed.fasta \
    # 2 7 7 80 10 50 500 -f -d -m -h -ngs \
    # > ${base_name}_R1.trf_out.dat

### -f: flanking sequence around each repeat is recorded in the alignment file
### -d: A data file (.dat) is produced. This file is a text file which contains 
### the same information, in the same order, as the summary table file, 
### plus consensus pattern and repeat sequences.
### -m: instructs the program to generate a masked sequence file. 
### The masked sequence file is a FASTA format file containing a copy of the sequence 
### with every location that occurred in a tandem repeat changed to the letter N. 
### The word <<masked>> is added to the sequence description line just after the > character.
### -h: suppress HTML output 
### -ngs: More compact .dat output on multisequence files, returns 0 on success

## 3.3 Remove repeats identified by TRF from FASTQ files in parallel
# cd ${project_home_dir}
# ### Find all files in the cutadapt_output_dir directory
# trimmed_files=$(find ${project_home_dir}/${cutadapt_output_dir} -type f \( -name "*R1.trimmed.fastq.gz"\
#  -o -name "*R2.trimmed.fastq.gz" \))
# python3 code/python-scripts/parallel_remove_repeats_from_fastq.py \
#   --remove-repeats-script code/python-scripts/remove_repeats_from_fastq.py \
#   --input-file-dir ${cutadapt_output_dir} \
#   --input-file-suffix trimmed.fastq.gz \
#   --trf-file-dir ${trf_output_dir} \
#   --trf-file-suffix trf_out.dat \
#   --output-dir ${fastq_norepeats_dir} \
#   --output-file-suffix clean.fastq.gz \
#   --nthreads ${nthreads} \
#   --files ${trimmed_files}
# parallel_remove_repeats_from_fastq.py usage: 
#   python3 parallel_remove_repeats_from_fastq.py
#   <remove_repeats.py> <cutadapt_output_dir> <trf_output_dir>
#   <fastq_norepeats_dir> <num_processes>

# remove_repeats_from_fastq.py usage: 
#   python3 code/python-scripts/remove_repeats_from_fastq.py \
#   ${cutadapt_output_dir}/${base_name}_R1.trimmed.fastq.gz \
#   ${trf_output_dir}/${base_name}_R1.trf_out.dat \
#   ${fastq_norepeats_dir}/${base_name}_R1.clean.fastq.gz 

# trimmed fastq/fastq.gz to be cleaned: {cutadapt_output_dir}/2D10_wms_L1_1.trimmed.fastq.gz 
# trf output dat file from TRF: {trf_output_dir}/2D10_wms_L1_1.trf_out.dat
# name of cleaned fastq/fastq.gz: ${fastq_norepeats_dir}/2D10_wms_L1_1.clean.fastq.gz

# 4 Bowtie2
# The Bowtie2 pipeline consists of these steps:
### 4.1 Align the reads to the reference genome (remove host DNA with bowtie2)
### 4.2 Convert file .sam to .bam
### 4.3 Filter unmapped reads (unmapped to host genome)
### 4.4 Split paired-end reads into separated fastq files .._R1 .._R2
### sort bam file by read name ( -n ) to have paired reads next to each other 
### 5.1 Run FastQC and MultiQC on decontaminated data as a final check
### 5.2 Run MultiQC on decontaminated data

### RUN WITHOUT TRIMMOMATIC OR TRF
### FILE is the entire path (data/fastq/yasuda-fastq/2D10_wms_L3_4.fq.gz)
### SAMPLE is path and file without extension (data/fastq/yasuda-fastq/2D10_wms)
### base_name is just the file without extension or path (2D10_wms)

# 4.1 Align the reads to the reference genome (remove host DNA with bowtie2)
### On raw data (don't do it!)
# for FILE in ${fastq_dir}/*L3_1.fq.gz
# do 
#   SAMPLE=$(echo ${FILE} | sed "s/_L3_1\.fq\.gz//")
#   base_name=$(basename "$SAMPLE" )
#   bowtie2 -p ${nthreads} -x ${bowtie2_indices_dir}/all_hosts_reference \
#   -1 ${fastq_dir}/${base_name}_L3_1.fq.gz \
#   -2 ${fastq_dir}/${base_name}_L3_2.fq.gz \
#   -S ${bowtie2_output_sam_dir}/${base_name}_mapped_and_unmapped.sam;
# done

### On trimmed data (from cutadapt)
for FILE in ${cutadapt_output_dir}/*R1.trimmed.fastq.gz
do 
  SAMPLE=$(echo ${FILE} | sed "s/_R1\.trimmed\.fastq\.gz//")
  base_name=$(basename "$SAMPLE" )
  bowtie2 -p ${nthreads} -x ${bowtie2_indices_dir}/all_hosts_reference \
  -1 ${cutadapt_output_dir}/${base_name}_R1.trimmed.fastq.gz \
  -2 ${cutadapt_output_dir}/${base_name}_R2.trimmed.fastq.gz \
  -S ${bowtie2_output_sam_dir}/${base_name}_trim_mapped_and_unmapped.sam;
done

### On trimmed data where tandem repeats were removed (TRF)
# for FILE in ${fastq_norepeats_dir}/*R1.clean.fastq.gz
# do 
#   SAMPLE=$(echo ${FILE} | sed "s/_R1\.clean\.fastq\.gz//")
#   base_name=$(basename "$SAMPLE" )
#   bowtie2 -p ${nthreads} -x ${bowtie2_indices_dir}/all_hosts_reference \
#   -1 ${fastq_norepeats_dir}/${base_name}_R1.clean.fastq.gz \
#   -2 ${fastq_norepeats_dir}/${base_name}_R2.clean.fastq.gz \
#   -S ${bowtie2_output_sam_dir}/${base_name}_mapped_and_unmapped.sam;
# done


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
### (${nthreads_sort} parallel threads, each using up to 5G memory)
for FILE in ${bowtie2_filtered_bam_dir}/*bothReadsUnmapped.bam 
do 
  SAMPLE=$(echo ${FILE} | sed "s/_bothReadsUnmapped\.bam//")
  base_name=$(basename "$SAMPLE" )
  samtools sort -n -m ${mem_req_sort} -@ ${nthreads_sort} \
  ${bowtie2_filtered_bam_dir}/${base_name}_bothReadsUnmapped.bam \
    -o ${bowtie2_sorted_bam_dir}/${base_name}_bothReadsUnmapped_sorted.bam
  samtools fastq -@ ${nthreads_sort} ${bowtie2_sorted_bam_dir}/${base_name}_bothReadsUnmapped_sorted.bam \
  	-1 ${bowtie2_decontam_fastq_dir}/${base_name}_decontam_R1.fastq.gz \
  	-2 ${bowtie2_decontam_fastq_dir}/${base_name}_decontam_R2.fastq.gz \
  	-0 /dev/null -s /dev/null -n;
done

# 5.1 Run FastQC and MultiQC on decontaminated data as a final check
# fastqc ${bowtie2_decontam_fastq_dir}/2D10_wms_decontam_R1.fastq.gz --outdir fastqc_output
fastqc ${bowtie2_decontam_fastq_dir}/*.fastq.gz --outdir ${fastqc_output_decontam_dir}
## 5.2 Run MultiQC on decontaminated data
multiqc ${fastqc_output_decontam_dir} --outdir ${multiqc_output_dir}

