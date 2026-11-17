#!/usr/bin/env bash
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate qc-tools 
set -euo pipefail
shopt -s nullglob
source config/bash/config.sh
cd "${project_home_dir}"
SAMPLES_PATHS_RAW="${1}"
GENOMES_METADATA="${2}"
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

# Function that checks if a required file exists or not
require() { [[ -s "$1" ]] || { echo "ERROR: required file missing: $1" >&2; exit 1; }; }
# Function that marks a sample as done (in a hidden directory)
is_done()   { [[ -e "${mark_dir}/${2}.${1}.done" ]]; }   # is_done STEP SAMPLE
mark_done() { touch "${mark_dir}/${2}.${1}.done"; }      # mark_done STEP SAMPLE

# Read the input file list
[[ -s "${SAMPLES_PATHS_RAW}" ]] || { echo "ERROR: metadata $SAMPLES_PATHS_RAW missing or empty" >&2; exit 1; }

# Remove CR LF
sed -i 's/\r$//' "${SAMPLES_PATHS_RAW}"

# We need to know the abosute paths of input FASTQ files.
declare -A R1 R2
SAMPLES=()
while IFS=$'\t' read -r sample_id fwd_absolute_path rev_absolute_path \
    || [[ -n "${sample_id}" ]]
do
    [[ "${sample_id}" == "sample_id" || -z "${sample_id}" ]] && continue
    r1="${fwd_absolute_path}"
    r2="${rev_absolute_path}"
    for f in "${r1}" "${r2}"; do
        [[ -s "${f}" ]] || { echo "ERROR: missing input for ${sample_id}: ${f}" >&2; exit 1; } # True if file exists and has a size greater than zero.
    done
    SAMPLES+=("${sample_id}"); R1[$sample_id]="${r1}"; R2[$sample_id]="${r2}"
done <  "${SAMPLES_PATHS_RAW}"

(( ${#SAMPLES[@]} > 0 )) || { echo "ERROR: no samples found" >&2; exit 1; } # Counts the number of elements in the final array of samples. If
                                                    # the array is empty, raises an error


# Download reference genomes for decontamination
# Heter_glaber.v1.7_hic_pac
# DMR_v1.0_HiC Damaraland mole-rat
# Human T2T
### Mouse GRCm39
### human decoy genome (downloaded separately)
# species_name	assembly_name	ncbi_id
# Heterocephalus glaber	Heter_glaber.v1.7_hic_pac	GCA_014060925.1
# Fukomys damarensis	DMR_v1.0_HiC	GCF_012274545.1
# Homo sapiens	Human_T2T-CHM13v2.0	GCF_009914755.1
# Mus musculus	GRCm39	GCF_000001635.27

while IFS=$'\t' read -r species_name assembly_name ncbi_id \
    || [[ -n "${species_name}" ]]
do
    [[ "${species_name}" == "species_name"|| -z "${species_name}" ]] && continue
    # Remove CR LF
    ncbi_id=$(echo "${ncbi_id}" | sed  's/\r$//' )
	assembly_name=$(echo "${assembly_name}" | sed  's/ //g' )
    echo "$(date +"%F %H:%M:%S")"
    GENOME_ID="${ncbi_id}__${assembly_name}"
    GENOME_FNA="${reference_genomes_dir}/${GENOME_ID}"__genomic.fna
    if [[ ! -s "${GENOME_FNA}" ]]; then
        curl https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/"${ncbi_id}"/download?include_annotation_type=GENOME_FASTA \
			--output "${reference_genomes_dir}/${GENOME_ID}".zip
        unzip "${reference_genomes_dir}/${GENOME_ID}".zip -d "${reference_genomes_dir}"
        mv "${reference_genomes_dir}/ncbi_dataset/data/${ncbi_id}/${ncbi_id}_${assembly_name}_genomic.fna" \
            "${GENOME_FNA}"
        rm "${reference_genomes_dir}"/md5sum.txt
        rm "${reference_genomes_dir}"/README.md
        rm -rf "${reference_genomes_dir}"/ncbi_dataset
        rm "${reference_genomes_dir}/${GENOME_ID}".zip

    fi
done < "${GENOMES_METADATA}"

### human decoy genome
GENOME_FNA="${reference_genomes_dir}"/hs37d5ss__human_decoy.fna
if [[ ! -s "${GENOME_FNA}" ]]; then
	curl ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5ss.fa.gz \
		--output "${reference_genomes_dir}"/hs37d5ss.fa.gz
	gunzip "${reference_genomes_dir}"/hs37d5ss.fa.gz
	mv "${reference_genomes_dir}"/hs37d5ss.fa "${GENOME_FNA}"
fi

conda activate nfkappab_evo-tools # need a different env to avoid conflicts
# phiX174 genome
GENOME_FNA="${reference_genomes_dir}"/NC_001422.1__phiX174__genomic.fna
if [[ ! -s "${GENOME_FNA}" ]]; then
	datasets download virus genome accession NC_001422.1 --include genome \
		--filename "${reference_genomes_dir}"/NC_001422.1__phiX174.zip
	unzip "${reference_genomes_dir}"/NC_001422.1__phiX174.zip -d "${reference_genomes_dir}"
	mv "${reference_genomes_dir}"/ncbi_dataset/data/genomic.fna \
		"${GENOME_FNA}"
	rm "${reference_genomes_dir}"/md5sum.txt
	rm "${reference_genomes_dir}"/README.md
	rm -rf "${reference_genomes_dir}"/ncbi_dataset
	rm "${reference_genomes_dir}"/NC_001422.1__phiX174.zip
fi
conda deactivate 

# Add a prefix with assembly ID just in case there're duplicates
for GENOME_FILE in "${reference_genomes_dir}"/*__genomic.fna
do
	GENOME_ID=$( basename "${GENOME_FILE}" | sed "s/__genomic\.fna//" )	
	GENOME_FILE_WITH_ASM_ID="${reference_genomes_dir}/${GENOME_ID}__genomic_id_replaced.fna"
	if [[ ! -s "${GENOME_FILE_WITH_ASM_ID}" ]]; then
		seqkit replace -p ^ -r "${GENOME_ID}__" "${GENOME_FILE}" > "${GENOME_FILE_WITH_ASM_ID}"
		rm "${GENOME_FILE}"
	fi

done
# separately for decoy
if [[ ! -s "${reference_genomes_dir}/hs37d5ss__human_decoy__id_replaced.fna" ]]; then
	seqkit replace -p ^ -r "hs37d5ss__human_decoy__" "${reference_genomes_dir}"/hs37d5ss__human_decoy.fna >\
		"${reference_genomes_dir}/hs37d5ss__human_decoy__id_replaced.fna"
fi

# Prepare NMR genome for Kraken2 database (old file name)
if [[ ! -s "${reference_genomes_dir}/Heter_glaber.v1.7_hic_pac_genomic_kraken2.fna" ]]; then
	sed 's/RPGA[0-9]\+\.1 Heterocephalus glaber isolate MEF-2018 //' \
		"${reference_genomes_dir}/GCA_014060925.1__Heter_glaber.v1.7_hic_pac__genomic.fna" > \
		"${reference_genomes_dir}/sed1_Heter_glaber.v1.7_hic_pac_genomic.fna"
	sed 's/, whole genome shotgun sequence/_Heter_glaber.v1.7_hic_pac|kraken:taxid|10181/' \
		"${reference_genomes_dir}/sed1_Heter_glaber.v1.7_hic_pac_genomic.fna" > \
		"${reference_genomes_dir}/Heter_glaber.v1.7_hic_pac_genomic_kraken2.fna"
fi

# Merge them all into one FASTA file
if [[ ! -s "${reference_genomes_dir}/${bowtie2_index_name}".fasta  ]]; then
	cat "${reference_genomes_dir}"/*_id_replaced.fna > "${reference_genomes_dir}/${bowtie2_index_name}".fasta 
fi

# Build a bowtie2 large index
if [[ ! -s "${bowtie2_indices_dir}/${bowtie2_index_name}"  ]]; then #not sure if this is the right file name to check
	bowtie2-build --large-index --threads "${nthreads_qc}" "${reference_genomes_dir}/${bowtie2_index_name}".fasta \
		"${bowtie2_indices_dir}/${bowtie2_index_name}"
fi

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
for SAMPLE_ID in "${SAMPLES[@]}" 
do
	echo "Running FastQC on ${SAMPLE_ID}"
	fastqc "${R1[$SAMPLE_ID]}" "${R2[$SAMPLE_ID]}" --outdir "${fastqc_raw_dir}"
done
## 1.1 Run MultiQC
multiqc "${fastqc_raw_dir}"/ --outdir "${multiqc_raw_dir}"

# 2. Run Cutadapt
### Three filtering steps:
### (1) Remove reads containing adapters.
### (2) Remove reads containing N > 10% (N represents the base cannot be determined).
### (3) Remove reads containing low quality (Qscore<= 5) base which is over 50% of the total base.

for SAMPLE_ID in "${SAMPLES[@]}" 
do
	echo "Running Cutadapt on ${SAMPLE_ID}"
	cutadapt \
		--cores="${nthreads_qc}" \
		-g "${FWD_ADAPTER}" -G "${REV_ADAPTER}" \
		--max-n 0.1 \
		-q 5 \
		-O 5 \
		--minimum-length 75 \
		--discard-trimmed \
		-o "${cutadapt_output_dir}"/"${SAMPLE_ID}"_R1.trimmed.fastq.gz \
		-p "${cutadapt_output_dir}"/"${SAMPLE_ID}"_R2.trimmed.fastq.gz \
		"${R1[$SAMPLE_ID]}" "${R2[$SAMPLE_ID]}" \
		"${cutadapt_output_dir}"/"${SAMPLE_ID}"_report.txt;
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
### On trimmed data (from cutadapt)
for SAMPLE_ID in "${SAMPLES[@]}" 
do
	echo "Running Bowtie2 on ${SAMPLE_ID}"
	bowtie2 -p "${nthreads_qc}" \
		-x "${bowtie2_indices_dir}/${bowtie2_index_name}" \
		-1 "${cutadapt_output_dir}"/"${SAMPLE_ID}"_R1.trimmed.fastq.gz \
		-2 "${cutadapt_output_dir}"/"${SAMPLE_ID}"_R2.trimmed.fastq.gz \
		-S "${bowtie2_output_sam_dir}"/"${SAMPLE_ID}"_trim_mapped_and_unmapped.sam
	## 3.2 Convert file .sam to .bam
	samtools view -bS ${bowtie2_output_sam_dir}/${SAMPLE_ID}_mapped_and_unmapped.sam > \
		${bowtie2_output_bam_dir}/${SAMPLE_ID}_mapped_and_unmapped.bam
	# 3.3 Filter unmapped reads (unmapped to host genome)
	### SAMtools SAM-flag filter: get unmapped pairs (both reads R1 and R2 unmapped)
	###-f  12    # Extract only ( -f ) alignments with both reads unmapped: <read unmapped><mate unmapped>
	###-F 256    # Do not (  -F  ) extract alignments which are: <not primary alignment>
	samtools view -b -f 12 -F 256 \
		${bowtie2_output_bam_dir}/${SAMPLE_ID}_mapped_and_unmapped.bam \
		> ${bowtie2_filtered_bam_dir}/${SAMPLE_ID}_bothReadsUnmapped.bam
	## 3.4 Split paired-end reads into separated fastq files .._R1 .._R2
	### First, sort bam file by read name ( -n ) to have paired reads next to each other 
	### (${nthreads_sort} parallel threads, each using up to 5G memory)
	### Then, filter and split
	samtools sort -n -m ${mem_req_sort} -@ ${nthreads_sort} \
	${bowtie2_filtered_bam_dir}/${SAMPLE_ID}_bothReadsUnmapped.bam \
		-o ${bowtie2_sorted_bam_dir}/${SAMPLE_ID}_bothReadsUnmapped_sorted.bam
	samtools fastq -@ ${nthreads_sort} ${bowtie2_sorted_bam_dir}/${SAMPLE_ID}_bothReadsUnmapped_sorted.bam \
		-1 ${bowtie2_decontam_fastq_dir}/${SAMPLE_ID}_decontam_R1.fastq.gz \
		-2 ${bowtie2_decontam_fastq_dir}/${SAMPLE_ID}_decontam_R2.fastq.gz \
		-0 /dev/null -s /dev/null -n;
done

# 4.1 Run FastQC and MultiQC on decontaminated data as a final check
# fastqc ${bowtie2_decontam_fastq_dir}/2D10_wms_decontam_R1.fastq.gz --outdir fastqc_output
fastqc ${bowtie2_decontam_fastq_dir}/*.fastq.gz --outdir ${fastqc_decontam_dir}
## 4.2 Run MultiQC on decontaminated data
multiqc ${fastqc_decontam_dir} --outdir ${multiqc_decontam_dir}
