#!/usr/bin/env bash
fastq_dir=data/fastq/yasuda-fastq
nthreads=12
mem_req=8G
# Rename files
for FILE in nmrF*.fq.gz; do mv "$FILE" $(echo "$FILE" | sed 's/nmrF_//'); done
for FILE in *DKDN230040401-1A_HC3LYDSX7*.fq.gz; 
  do mv "$FILE" $(echo "$FILE" | sed 's/DKDN230040401-1A_HC3LYDSX7/wms/'); done
for FILE in *DKDN230040401-1A_HC557DSX7*.fq.gz; 
  do mv "$FILE" $(echo "$FILE" | sed 's/DKDN230040401-1A_HC557DSX7/wms/'); done
for FILE in *DKDN230040402-1A_HC557DSX7*.fq.gz; 
  do mv "$FILE" $(echo "$FILE" | sed 's/DKDN230040402-1A_HC557DSX7/wms/'); done
for FILE in *DKDN230040410-1A_HC52YDSX7*.fq.gz; 
  do mv "$FILE" $(echo "$FILE" | sed 's/DKDN230040410-1A_HC52YDSX7/wms/'); done
for FILE in *DKDN230040409-1A_HC3LYDSX7*.fq.gz; 
  do mv "$FILE" $(echo "$FILE" | sed 's/DKDN230040409-1A_HC3LYDSX7/wms/'); done
for FILE in *DKDN230040409-1A_HC557DSX7*.fq.gz; 
  do mv "$FILE" $(echo "$FILE" | sed 's/DKDN230040409-1A_HC557DSX7/wms/'); done
for FILE in *DKDN230040411-1A_HC52YDSX7*.fq.gz; 
  do mv "$FILE" $(echo "$FILE" | sed 's/DKDN230040411-1A_HC52YDSX7/wms/'); done
for FILE in *DKDN230040407-1A_HC557DSX7*.fq.gz; 
  do mv "$FILE" $(echo "$FILE" | sed 's/DKDN230040407-1A_HC557DSX7/wms/'); done
for FILE in *DKDN230040408-1A_HC557DSX7*.fq.gz; 
  do mv "$FILE" $(echo "$FILE" | sed 's/DKDN230040408-1A_HC557DSX7/wms/'); done
for FILE in *DKDN230040406-1A_HC557DSX7*.fq.gz; 
  do mv "$FILE" $(echo "$FILE" | sed 's/DKDN230040406-1A_HC557DSX7/wms/'); done
for FILE in *DKDN230040403-1A_HC557DSX7*.fq.gz; 
  do mv "$FILE" $(echo "$FILE" | sed 's/DKDN230040403-1A_HC557DSX7/wms/'); done
for FILE in *DKDN230040404-1A_HC557DSX7*.fq.gz; 
  do mv "$FILE" $(echo "$FILE" | sed 's/DKDN230040404-1A_HC557DSX7/wms/'); done
for FILE in *DKDN230040405-1A_HC3LYDSX7*.fq.gz; 
  do mv "$FILE" $(echo "$FILE" | sed 's/DKDN230040405-1A_HC3LYDSX7/wms/'); done
for FILE in *DKDN230040405-1A_HC557DSX7*.fq.gz; 
  do mv "$FILE" $(echo "$FILE" | sed 's/DKDN230040405-1A_HC557DSX7/wms/'); done

source ~/miniconda3/etc/profile.d/conda.sh
# Set conda channel priority
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels biobakery
mkdir ~/common_data/reference_genomes
reference_genomes_dir=~/common_data/reference_genomes
bowtie2_indices_dir=~/common_data/bowtie2_indices
# Download reference genomes
### Heter_glaber.v1.7_hic_pac
curl https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCA_014060925.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT \
  --output ${reference_genomes_dir}/Heter_glaber.v1.7_hic_pac.zip
### Human T2T
curl https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_009914755.1/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT \
  --output ${reference_genomes_dir}/T2T.zip
### Mouse GRCm39
curl https://api.ncbi.nlm.nih.gov/datasets/v2alpha/genome/accession/GCF_000001635.27/download?include_annotation_type=GENOME_FASTA,GENOME_GFF,RNA_FASTA,CDS_FASTA,PROT_FASTA,SEQUENCE_REPORT \
  --output ${reference_genomes_dir}/GRCm39.zip
### human decoy genome
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5ss.fa.gz

cat ${reference_genomes_dir}/*genomic.fna > ${reference_genomes_dir}/all_hosts_reference.fasta 
cat ${reference_genomes_dir}/phiX174.fasta >>${reference_genomes_dir}/all_hosts_reference.fasta

conda create --name qc-tools -c bioconda bowtie2 trimmomatic samtools trf fastqc multiqc seqtk
conda activate qc-tools 
bowtie2-build --large-index --threads 12 ${reference_genomes_dir}/all_hosts_reference.fasta \
  ${bowtie2_indices_dir}/all_hosts_reference

# Run FastQC
mkdir output/fastqc_output
fastqc ${fastq_dir}/2D10_wms_L1_1.fq.gz --outdir output/fastqc_output
fastqc ${fastq_dir}/*.fq.gz --outdir output/fastqc_output
## Run MultiQC
multiqc ./output/fastqc_output/

# Run Trimmomatic
mkdir output/trimmomatic_output
trimmomatic_output_dir=output/trimmomatic_output
java -jar ~/miniconda3/envs/qc-tools/share/trimmomatic/trimmomatic.jar \
	PE -threads 4 -trimlog ${trimmomatic_output_dir}/trimmomatic_seq.log \
	input/seq_R1.fastq input/seq_R2.fastq \
	${trimmomatic_output_dir}/seq_R1.trimmed.fastq ${trimmomatic_output_dir}/seq_R1un.trimmed.fastq \
	${trimmomatic_output_dir}/seq_R2.trimmed.fastq ${trimmomatic_output_dir}/seq_R2un.trimmed.fastq \
	ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:TRUE SLIDINGWINDOW:4:20

# on fastq.gz files
java -jar ~/miniconda3/envs/qc-tools/share/trimmomatic/trimmomatic.jar \
	PE -threads 4 -trimlog ${trimmomatic_output_dir}/trimmomatic_seq.gz.log \
	input/seq_R1.fastq.gz input/seq_R2.fastq.gz \
	${trimmomatic_output_dir}/seq_R1.trimmed.fastq.gz ${trimmomatic_output_dir}/seq_R1un.trimmed.fastq.gz \
	${trimmomatic_output_dir}/seq_R2.trimmed.fastq.gz ${trimmomatic_output_dir}/seq_R2un.trimmed.fastq.gz \
	ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:TRUE SLIDINGWINDOW:4:20

java -jar ~/miniconda3/envs/qc-tools/share/trimmomatic/trimmomatic.jar \
	PE -threads 4 -trimlog ${trimmomatic_output_dir}/trimmomatic_seq.gz.log \
	${fastq_dir}/2D10_wms_L1_1.fq.gz ${fastq_dir}/2D10_wms_L1_2.fq.gz  \
	${trimmomatic_output_dir}/2D10_wms_L1_1.trimmed.fastq.gz ${trimmomatic_output_dir}/2D10_wms_L1_1un.trimmed.fastq.gz \
	${trimmomatic_output_dir}/2D10_wms_L1_2.trimmed.fastq.gz ${trimmomatic_output_dir}/2D10_wms_L1_2un.trimmed.fastq.gz \
	ILLUMINACLIP:Adapters.fa:2:30:10:8:TRUE SLIDINGWINDOW:4:20



'
Here is the order of file names in the trimmomatic command
input_forward.fq.gz input_reverse.fq.gz 
output_forward_paired.fq.gz output_forward_unpaired.fq.gz 
output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz
ILLUMINACLIP:adapters SLIDINGWINDOW
'

# TRF
mkdir output/temp_fasta
temp_fasta_dir=output/temp_fasta
## Create temporary fasta files
### using fastq files
seqtk seq -a ${trimmomatic_output_dir}/seq_R1.trimmed.fastq > ${temp_fasta_dir}/seq_R1.trimmed.fasta
seqtk seq -a ${trimmomatic_output_dir}/seq_R2.trimmed.fastq > ${temp_fasta_dir}/seq_R2.trimmed.fasta
### using fastq.gz files
seqtk seq -a ${trimmomatic_output_dir}/seq_R1.trimmed.fastq.gz > ${temp_fasta_dir}/seq_R1.trimmed.fasta
seqtk seq -a ${trimmomatic_output_dir}/seq_R2.trimmed.fastq.gz > ${temp_fasta_dir}/seq_R2.trimmed.fasta

seqtk seq -a ${trimmomatic_output_dir}/2D10_wms_L1_1.trimmed.fastq.gz > ${temp_fasta_dir}/2D10_wms_L1_1.trimmed.fasta
seqtk seq -a ${trimmomatic_output_dir}/2D10_wms_L1_2.trimmed.fastq.gz > ${temp_fasta_dir}/2D10_wms_L1_2.trimmed.fasta
## Run TRF
'
-f: flanking sequence around each repeat is recorded in the alignment file
-d: A data file (.dat) is produced. This file is a text file which contains 
the same information, in the same order, as the summary table file, 
plus consensus pattern and repeat sequences.
-m: instructs the program to generate a masked sequence file. 
The masked sequence file is a FASTA format file containing a copy of the sequence 
with every location that occurred in a tandem repeat changed to the letter N. 
The word <<masked>> is added to the sequence description line just after the > character.
-h: suppress HTML output 
-ngs: More compact .dat output on multisequence files, returns 0 on success
'
mkdir output/trf_output
cd output/trf_output
trf ../../${temp_fasta_dir}/seq_R1.trimmed.fasta 2 7 7 80 10 50 500 -f -d -m -h -ngs > seq_R1.trf_out.dat
trf ../../${temp_fasta_dir}/seq_R2.trimmed.fasta 2 7 7 80 10 50 500 -f -d -m -h -ngs > seq_R2.trf_out.dat

# from fastq.gz
trf ../../${temp_fasta_dir}/seq_R1.gz.trimmed.fasta 2 7 7 80 10 50 500 -f -d -m -h -ngs > seq_R1.trf_out.dat
trf ../../${temp_fasta_dir}/seq_R2.gz.trimmed.fasta 2 7 7 80 10 50 500 -f -d -m -h -ngs > seq_R2.trf_out.dat

trf ../../${temp_fasta_dir}/2D10_wms_L1_1.trimmed.fasta 2 7 7 80 10 50 500 -f -d -m -h -ngs > 2D10_wms_L1_1.trf_out.dat
trf ../../${temp_fasta_dir}/2D10_wms_L1_2.trimmed.fasta 2 7 7 80 10 50 500 -f -d -m -h -ngs > 2D10_wms_L1_2.trf_out.dat


## Remove repeats identified by TRF from FASTQ files
cd ..
#awk '/^@/ {print substr($0, 1)}' trf_output/seq_R2.trf_out.dat > trf_output/seq_R2.trf_list.txt
### /^@/: This specifies a pattern to match lines that start with "@".
### {print substr($0, 1)}: This is the action to perform when a matching line is found. 
### It prints the substring of the line starting from the first character.

#grep "^@" trf_output/seq_R2.trf_out.dat > trf_output/seq_R2.trf_list.txt

mkdir output/fastq_norepeats
# fastq
python3 code/python-scripts/remove_repeats_from_fastq.py \
${trimmomatic_output_dir}/seq_R2.trimmed.fastq \
trf_output/seq_R2.trf_out.dat \
fastq_norepeats/seq_R2.clean.fastq 

# fastq.gz
python3 code/python-scripts/remove_repeats_from_fastq.py \
${trimmomatic_output_dir}/seq_R1.trimmed.fastq.gz \
trf_output/seq_R1.trf_out.dat \
fastq_norepeats/seq_R1.clean.fastq.gz
python3 code/python-scripts/remove_repeats_from_fastq.py \
${trimmomatic_output_dir}/seq_R2.trimmed.fastq.gz \
trf_output/seq_R2.trf_out.dat \
fastq_norepeats/seq_R2.clean.fastq.gz

# My data
python3 code/python-scripts/remove_repeats_from_fastq.py \
${trimmomatic_output_dir}/2D10_wms_L1_1.trimmed.fastq.g \
trf_output/2D10_wms_L1_1.trf_out.dat \
fastq_norepeats/2D10_wms_L1_1.clean.fastq.gz


# Bowtie2
# 1. Bowtie2 mapping against host sequence
# RUN WITHOUT TRIMMOMATIC OR TRF
# done in bowtie2-pipeline.sh
bowtie2_host_removed_fastq_dir=data/bowtie2_host_removed_fastq

# Run FastQC 
mkdir output/fastqc_output_decontam
fastqc ${bowtie2_host_removed_fastq_dir}/2D10_wms_host_removed_R1.fastq.gz --outdir fastqc_output
fastqc ${bowtie2_host_removed_fastq_dir}/*.fastq.gz --outdir output/fastqc_output_decontam
## Run MultiQC on decontaminated data
multiqc ./fastqc_output_decontam/