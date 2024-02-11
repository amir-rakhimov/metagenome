#!/usr/bin/env bash
# on fastq.gz files
java -jar ~/miniconda3/envs/qc-tools/share/trimmomatic/trimmomatic.jar \
	PE -threads ${nthreads} -trimlog ${trimmomatic_output_dir}/trimmomatic_seq.log \
	input/seq_R1.fastq input/seq_R2.fastq \
	${trimmomatic_output_dir}/seq_R1.trimmed.fastq ${trimmomatic_output_dir}/seq_R1un.trimmed.fastq \
	${trimmomatic_output_dir}/seq_R2.trimmed.fastq ${trimmomatic_output_dir}/seq_R2un.trimmed.fastq \
	ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:8:TRUE SLIDINGWINDOW:4:20
java -jar ~/miniconda3/envs/qc-tools/share/trimmomatic/trimmomatic.jar \
	PE -threads ${nthreads} \
  -trimlog ${trimmomatic_output_dir}/trimmomatic_seq.log \
	${fastq_dir}/${base_name}_wms_L3_1.fq.gz ${fastq_dir}/${base_name}_wms_L3_2.fq.gz \
	${trimmomatic_output_dir}/${base_name}_R1.trimmed.fastq.gz ${trimmomatic_output_dir}/${base_name}_R1un.trimmed.fastq.gz \
	${trimmomatic_output_dir}/${base_name}_R2.trimmed.fastq.gz ${trimmomatic_output_dir}/${base_name}_R2un.trimmed.fastq.gz \
	ILLUMINACLIP:${adapters_file}:2:30:10:8:TRUE SLIDINGWINDOW:4:20

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
### Here is the order of file names in the trimmomatic command
### input_forward.fq.gz input_reverse.fq.gz 
### output_forward_paired.fq.gz output_forward_unpaired.fq.gz 
### output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz
### ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip threshold>:
### <simple clip threshold>:<minAdapterLength>:<keepBothReads> 
### SLIDINGWINDOW:<windowSize>:<requiredQuality>


# TRF
## Create temporary fasta files
### using fastq files
seqtk seq -a ${trimmomatic_output_dir}/seq_R1.trimmed.fastq > ${temp_fasta_dir}/seq_R1.trimmed.fasta
seqtk seq -a ${trimmomatic_output_dir}/seq_R2.trimmed.fastq > ${temp_fasta_dir}/seq_R2.trimmed.fasta
### using fastq.gz files
seqtk seq -a ${trimmomatic_output_dir}/seq_R1.trimmed.fastq.gz > ${temp_fasta_dir}/seq_R1.trimmed.fasta
seqtk seq -a ${trimmomatic_output_dir}/seq_R2.trimmed.fastq.gz > ${temp_fasta_dir}/seq_R2.trimmed.fasta

seqtk seq -a ${trimmomatic_output_dir}/2D10_wms_L1_1.trimmed.fastq.gz > ${temp_fasta_dir}/2D10_wms_L1_1.trimmed.fasta
seqtk seq -a ${trimmomatic_output_dir}/2D10_wms_L1_2.trimmed.fastq.gz > ${temp_fasta_dir}/2D10_wms_L1_2.trimmed.fasta

# from fastq.gz
trf ~/projects/metagenome/${temp_fasta_dir}/seq_R1.gz.trimmed.fasta 2 7 7 80 10 50 500 -f -d -m -h -ngs \
  > seq_R1.trf_out.dat
trf ~/projects/metagenome/${temp_fasta_dir}/seq_R2.gz.trimmed.fasta 2 7 7 80 10 50 500 -f -d -m -h -ngs \
  > seq_R2.trf_out.dat

trf ~/projects/metagenome/${temp_fasta_dir}/2D10_wms_L1_1.trimmed.fasta 2 7 7 80 10 50 500 -f -d -m -h -ngs \
  > 2D10_wms_L1_1.trf_out.dat
trf ~/projects/metagenome/${temp_fasta_dir}/2D10_wms_L1_2.trimmed.fasta 2 7 7 80 10 50 500 -f -d -m -h -ngs \
  > 2D10_wms_L1_2.trf_out.dat

for FILE in ${project_home_dir}/${temp_fasta_dir}/*R1.trimmed.fasta
do 
  SAMPLE=$(echo ${FILE} | sed "s/_R1\.trimmed\.fasta//")
  base_name=$(basename "$SAMPLE" )
  trf ${project_home_dir}/${temp_fasta_dir}/${base_name}_R1.trimmed.fasta 2 7 7 80 10 50 500 -f -d -m -h -ngs \
    > ${base_name}_R1.trf_out.dat
  trf ${project_home_dir}/${temp_fasta_dir}/${base_name}_R2.trimmed.fasta 2 7 7 80 10 50 500 -f -d -m -h -ngs \
    > ${base_name}_R2.trf_out.dat;
done


## 3.2 Remove repeats identified by TRF from FASTQ files
cd ${project_home_dir}
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

for FILE in ${trf_output_dir}/*R1.trf_out.dat
do 
 SAMPLE=$(echo ${FILE} | sed "s/_R1\.trf_out\.dat//")
 base_name=$(basename "$SAMPLE" )
 python3 code/python-scripts/remove_repeats_from_fastq.py \
  ${cutadapt_output_dir}/${base_name}_R1.trimmed.fastq.gz \
  ${trf_output_dir}/${base_name}_R1.trf_out.dat \
  ${fastq_norepeats_dir}/${base_name}_R1.clean.fastq.gz 
 python3 code/python-scripts/remove_repeats_from_fastq.py \
  ${cutadapt_output_dir}/${base_name}_R2.trimmed.fastq.gz \
  ${trf_output_dir}/${base_name}_R2.trf_out.dat \
  ${fastq_norepeats_dir}/${base_name}_R2.clean.fastq.gz ;
done