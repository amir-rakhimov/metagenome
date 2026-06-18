#!/usr/bin/env bash
source $HOME/miniconda3/etc/profile.d/conda.sh

######
# This script performs the following:
# 1. Renames FASTQ files to make them shorter because the sequencing company gave the files long 
# names
# 2. Downloads reference genomes for decontamination. We are using genomes of
# naked mole-rat (Heter_glaber.v1.7_hic_pac), Damaraland mole-rat (DMR_v1.0_HiC), human (T2T),
# mouse (GRCm39), and human decoy genome (as specified in 
# https://www.cureffi.org/2013/02/01/the-decoy-genome/). 

# The reference genomes are merged into one FASTA file and saved in 
# the `$HOME/common_data/reference_genomes` directory 
# as `all_hosts_reference.fasta` file. BowTie2 also builds an index using `all_hosts_reference.fasta` 
# and saves it in the `$HOME/common_data/bowtie2_indices` directory.

# 3. Runs FASTQC and MultiQC on raw reads.

# Input: FASTQ files in the `data/fastq/yasuda-fastq` directory.

# Output: reports saved in the `output/qc_pipeline/fastqc_output`
# and `output/qc_pipeline/multiqc_output` directories.
# 4. Trims raw reads with Cutadapt. The parameters are `--max-n` = 0.1 (do not allow > 10% Ns),
# `-q` = 5 (remove bases with basequal < 5), `-O` = 5 (Require MINLENGTH overlap between
# read and adapter for an adapter to be found), --discard-trimmed (discard reads in
# which an adapter was found), `--minimum-length` = 75 (minimum length after trimming). 

# Input: FASTQ files in the `data/fastq/yasuda-fastq` directory.

# Output: FASTQ files saved in the `output/qc_pipeline/cutadapt_output` 
# directory. The Cutadapt command also saves a report for each sample in the same directory.

# 5. Aligns trimmed reads to the reference genomes with BowTie2.

# Input: the BowTie2 index (`$HOME/common_data/bowtie2_indices/all_hosts_reference`), trimmed
# FASTQ files (`output/qc_pipeline/cutadapt_output`).

# Output: SAM and BAM files with reads that are mapped and unmapped to reference genomes. The SAM
# files for each sample are saved in the `output/bowtie2_pipeline/bowtie2_output_sam` directory as
# `${sample_name}_trim_mapped_and_unmapped.sam`. For example, the SAM file for 2D10 sample would be 
# `2D10_trim_mapped_and_unmapped.sam`. The BAM files are saved in the
# `output/bowtie2_pipeline/bowtie2_output_bam` directory in the same format as SAM files.

# 6. Filters out read pairs that did not map to reference genomes. These are reads free of 
# contamination.

# Input: BAM files with mapped and unmapped reads from the previous step
# (`output/bowtie2_pipeline/bowtie2_output_bam`).

# Output: BAM files where both reads did not map to reference genomes. They are saved in the
# `output/bowtie2_pipeline/bowtie2_filtered_bam` directory and are named like 
# `${sample_name}_bothReadsUnmapped.bam`.

# 7. Sorts the BAM file by read name to have paired reads next to each other. Then, it splits
# paired-end reads into separated fastq files: .._R1 .._R2. 

# Input: BAM files where both reads did not map to reference genomes 
# (`output/bowtie2_pipeline/bowtie2_filtered_bam/${sample_name}_bothReadsUnmapped.bam`).

# Output: Sorted BAM file (`${base_name}_bothReadsUnmapped_sorted.bam`) in the 
# `output/bowtie2_pipeline/bowtie2_sorted_bam` directory; decontaminated reads 
# (`${sample_name}_decontam_R1.fastq.gz` and `${sample_name}_decontam_R2.fastq.gz`) in the 
# `data/bowtie2_decontam_fastq` directory.

# 8. Runs FastQC and MultiQC on decontaminated data as a final check.

# Input: decontaminated reads in the `data/bowtie2_decontam_fastq` directory. 

# Output reports are saved in the 
# `output/qc_pipeline/fastqc_output` and `output/qc_pipeline/multiqc_output` directories.
######

# Main Input:
# FASTQ files (gzip-compressed) with raw paired-end demultiplexed 
# reads (with or without adapters). If adapters are not removed, please specify
# adapter sequences in the `bash-scripts/qc-pipeline/001-qc-commands-with-bowtie2.sh` script.

# Final output: 
# * Quality check reports from FASTQC and MultiQC on raw data
# * FASTQ files with trimmed reads (gzip-compressed)
# * SAM and BAM files with reads that are mapped and unmapped to reference genomes 
# * BAM file with read pairs that are unmapped (both R1 and R2 are unmapped) to the reference
# genomes
# * Sorted BAM file with read pairs that are unmapped (sorted by read name)
# * Final FASTQ files with decontaminated reads: only read pairs that were unmapped to 
# reference genomes are retained in the final files
# * Final quality check reports from FASTQC and MultiQC
######


# Activate conda environment
conda activate qc-tools 
# Rename files: done in the rename script
cd "${project_home_dir}"

# Download reference genomes
## Heter_glaber.v1.7_hic_pac
curl https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCA_014060925.1/download?include_annotation_type=GENOME_FASTA \
  --output "${reference_genomes_dir}"/Heter_glaber.v1.7_hic_pac.zip
## There is also a new assembly that is used by Ensembl: Naked_mole-rat_maternal
curl https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCA_944319715.1/download?include_annotation_type=GENOME_FASTA \
  --output "${reference_genomes_dir}"/Naked_mole-rat_maternal.zip
## DMR_v1.0_HiC Damaraland mole-rat
curl https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_012274545.1/download?include_annotation_type=GENOME_FASTA \
  --output "${reference_genomes_dir}"/DMR_v1.0_HiC.zip
### Human T2T
curl https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_009914755.1/download?include_annotation_type=GENOME_FASTA \
  --output "${reference_genomes_dir}"/T2T.zip
### Mouse GRCm39
curl https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000001635.27/download?include_annotation_type=GENOME_FASTA \
  --output "${reference_genomes_dir}"/GRCm39.zip
### human decoy genome
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5ss.fa.gz
mv hs37d5ss.fa human_decoy.fa
# TODO: check and rewrite if needed
# for file in ncbi_dataset/data/*/protein.faa
# do
# directory_name=$(dirname $file)
# accession=$(basename $directory_name)
# mv "${file}" "${directory_name}/${accession}_$(basename $file)"
# done

mv ncbi_dataset/data/GCF_012274545.1/GCF_012274545.1_DMR_v1.0_HiC_genomic.fna DMR_v1.0_HiC_genomic.fna
# Prepare NMR genome for Kraken2 database
sed 's/RPGA[0-9]\+\.1 Heterocephalus glaber isolate MEF-2018 //' GCA_014060925.1_Heter_glaber.v1.7_hic_pac_genomic.fna > \
 sed1_Heter_glaber.v1.7_hic_pac_genomic.fna
sed 's/, whole genome shotgun sequence/_Heter_glaber.v1.7_hic_pac|kraken:taxid|10181/' sed1_Heter_glaber.v1.7_hic_pac_genomic.fna > \
 Heter_glaber.v1.7_hic_pac_genomic_kraken2.fna

# # Merge them all into one FASTA file
cat "${reference_genomes_dir}"/*genomic.fna > "${reference_genomes_dir}"/all_hosts_reference.fasta 
cat "${reference_genomes_dir}"/phiX174.fasta >>"${reference_genomes_dir}"/all_hosts_reference.fasta
# # Build a bowtie2 large index
bowtie2-build --large-index --threads "${nthreads_qc}" "${reference_genomes_dir}"/all_hosts_reference.fasta \
  "${bowtie2_indices_dir}"/all_hosts_reference

# The QC workflow is this:
### 1. Run FastQC and trim overrepresented sequences
### 1.1 Run MultiQC
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
fastqc "${fastq_dir}"/*.fq.gz --outdir "${fastqc_output_dir}"
## 1.1 Run MultiQC
multiqc "${fastqc_output_dir}"/ --outdir ${multiqc_output_dir}

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
  --cores="${nthreads_qc}" \
  -g "${FWD_ADAPTER}" -G "${REV_ADAPTER}" \
  --max-n 0.1 \
  -q 5 \
  -O 5 \
  --minimum-length 75 \
  --discard-trimmed \
  -o "${cutadapt_output_dir}"/"${base_name}"_R1.trimmed.fastq.gz \
  -p "${cutadapt_output_dir}"/"${base_name}"_R2.trimmed.fastq.gz \
  "${fastq_dir}"/"${base_name}"_wms_L3_1.fq.gz "${fastq_dir}"/"${base_name}"_wms_L3_2.fq.gz > \
  "${cutadapt_output_dir}"/"${base_name}"_report.txt;
done

### FILE: data/fastq/yasuda-fastq/2D10_wms_L3_1.fq.gz
### SAMPLE: data/fastq/yasuda-fastq/2D10
### base_name: 2D10 
# cutadapt \
#   --cores=${nthreads_qc} \                         ## number of cores for parallelisation
#   -g ${FWD_ADAPTER} -G ${REV_ADAPTER} \         ## 5' adapter
#   --max-n 0.1 \                                 ## do not allow > 10% Ns
#   -q 5 \                                        ## remove bases with basequal < 5
#   -O 5  \                ## Require MINLENGTH overlap between read and adapter for an adapter to be found.
#   --minimum-length (readlength/2) \             ## minimum length after trimming
#   --discard-trimmed \                           ## Discard reads in which an adapter was found
#   -o read1_trimmed.fq -p read2_trimmed.fq \     ## outputs
#   read1.fq read2.fq                             ## inputs


# 3. Bowtie2
# The Bowtie2 pipeline consists of these steps:
### 3.1 Align the reads to the reference genome (remove host DNA with bowtie2)
### 3.2 Convert file .sam to .bam
### 3.3 Filter unmapped reads (unmapped to host genome)
### 3.4 Split paired-end reads into separated fastq files .._R1 .._R2
### sort bam file by read name ( -n ) to have paired reads next to each other 
### 4.1 Run FastQC and MultiQC on decontaminated data as a final check
### 4.2 Run MultiQC on decontaminated data

### FILE is the entire path (data/fastq/yasuda-fastq/2D10_wms_L3_4.fq.gz)
### SAMPLE is path and file without extension (data/fastq/yasuda-fastq/2D10_wms)
### base_name is just the file without extension or path (2D10_wms)

# 3.1 Align the reads to the reference genome (remove host DNA with bowtie2)
### On raw data (don't do it!)
# for FILE in ${fastq_dir}/*L3_1.fq.gz
# do 
#   SAMPLE=$(echo ${FILE} | sed "s/_L3_1\.fq\.gz//")
#   base_name=$(basename "$SAMPLE" )
#   bowtie2 -p ${nthreads_qc} -x ${bowtie2_indices_dir}/all_hosts_reference \
#   -1 ${fastq_dir}/${base_name}_L3_1.fq.gz \
#   -2 ${fastq_dir}/${base_name}_L3_2.fq.gz \
#   -S ${bowtie2_output_sam_dir}/${base_name}_mapped_and_unmapped.sam;
# done

### On trimmed data (from cutadapt)
for FILE in ${cutadapt_output_dir}/*R1.trimmed.fastq.gz
do 
  SAMPLE=$(echo ${FILE} | sed "s/_R1\.trimmed\.fastq\.gz//")
  base_name=$(basename "$SAMPLE" )
  bowtie2 -p "${nthreads_qc}" -x ${bowtie2_indices_dir}/all_hosts_reference \
  -1 "${cutadapt_output_dir}"/"${base_name}"_R1.trimmed.fastq.gz \
  -2 "${cutadapt_output_dir}"/"${base_name}"_R2.trimmed.fastq.gz \
  -S "${bowtie2_output_sam_dir}"/"${base_name}"_trim_mapped_and_unmapped.sam;
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


## 3.2 Convert file .sam to .bam
for FILE in ${bowtie2_output_sam_dir}/*mapped_and_unmapped.sam
do 
  SAMPLE=$(echo ${FILE} | sed "s/_mapped_and_unmapped\.sam//")
  base_name=$(basename "$SAMPLE" )
  samtools view -bS ${bowtie2_output_sam_dir}/${base_name}_mapped_and_unmapped.sam > \
    ${bowtie2_output_bam_dir}/${base_name}_mapped_and_unmapped.bam;
done

## 3.3 Filter unmapped reads (unmapped to host genome)
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

## 3.4 Split paired-end reads into separated fastq files .._R1 .._R2
### First, sort bam file by read name ( -n ) to have paired reads next to each other 
### (${nthreads_sort} parallel threads, each using up to 5G memory)
### Then, filter and split
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

# 4.1 Run FastQC and MultiQC on decontaminated data as a final check
# fastqc ${bowtie2_decontam_fastq_dir}/2D10_wms_decontam_R1.fastq.gz --outdir fastqc_output
fastqc ${bowtie2_decontam_fastq_dir}/*.fastq.gz --outdir ${fastqc_output_decontam_dir}
## 4.2 Run MultiQC on decontaminated data
multiqc ${fastqc_output_decontam_dir} --outdir ${multiqc_output_dir}
