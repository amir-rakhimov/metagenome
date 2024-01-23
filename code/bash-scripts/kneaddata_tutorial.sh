# 1. Prepare the environment
cd projects
mkdir kneaddata_analysis
mkdir kneaddata_analysis
conda create --name kneaddata-tools -c bioconda bowtie2 kneaddata trimmomatic samtools trf fastqc multiqc
conda activate kneaddata-tools 

# 2. Download demo input files
cd kneaddata_analysis
wget https://github.com/biobakery/kneaddata/files/4703820/input.zip
unzip input.zip

### Fix trimmomatic error with jar
mv ~/miniconda3/envs/kneaddata-tools/bin/trimmomatic ~/miniconda3/envs/kneaddata-tools/bin/cmdline_trimmomatic
ln -s ~/miniconda3/envs/kneaddata-tools/share/trimmomatic/trimmomatic.jar ~/miniconda3/envs/kneaddata-tools/bin/trimmomatic
ln -s ~/miniconda3/envs/kneaddata-tools/share/trimmomatic/trimmomatic.jar ~/miniconda3/envs/kneaddata-tools/bin/trimmomatic.jar

# 3.1 Single end reads
kneaddata --unpaired input/singleEnd.fastq \
--reference-db input/demo_db \
--output kneaddataOutputSingleEnd
### add --run-trim-repetitive to trim overrepresented/repetitive sequences
### add --fastqc FastQC to create reports


## Kneaddata count table
kneaddata_read_count_table \
--input kneaddataOutputSingleEnd \
--output kneaddata_read_count_table.tsv

### Open the table
column -t -s $'\t' kneaddata_read_count_table.tsv | less -S

### check how many sequences were identified as contaminants 
wc -l kneaddataOutputSingleEnd/singleEnd_kneaddata_demo_db_bowtie2_contam.fastq

## FastQC
kneaddata --unpaired input/SE_extra.fastq \
--reference-db input/demo_db \
--output kneaddataOutputFastQC \
--run-fastqc-start \
--run-fastqc-end

# 3.2 Paired end reads
kneaddata --input1 input/seq1.fastq \
--input2 input/seq2.fastq \
--reference-db input/demo_db \
--output kneaddataOutputPairedEnd 
### add --run-trim-repetitive to trim overrepresented/repetitive sequences
### add --fastqc FastQC to create reports

## Kneaddata count table
kneaddata_read_count_table \
--input kneaddataOutputPairedEnd \
--output kneaddata_read_count_table_paired.tsv

## FastQC
kneaddata --unpaired input/SE_extra.fastq \
--reference-db input/demo_db \
--output kneaddataOutputFastQC \
--run-fastqc-start \
--run-fastqc-end 

# 4. Changing parameters
## Trimmomatic
### Change MINLEN to 90
kneaddata --unpaired input/singleEnd.fastq \
--output kneaddataOutputTrimmomatic \
--reference-db input/demo_db \
--trimmomatic-options="MINLEN:90"

kneaddata_read_count_table --input kneaddataOutputTrimmomatic/ --output kneaddata_read_count_table_trimmomatic.tsv
column -t -s $'\t' kneaddata_read_count_table_trimmomatic.tsv | less -S

### Keep all Kneaddata parameters except MINVALUE
kneaddata --unpaired input/singleEnd.fastq \
--output kneaddataOutputTrimmomatic2 \
--reference-db input/demo_db \
--trimmomatic-options="ILLUMINACLIP:/TruSeq3-SE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:90"

kneaddata_read_count_table --input kneaddataOutputTrimmomatic2/ \
--output kneaddata_read_count_table3.tsv

column -t -s $'\t' kneaddata_read_count_table3.tsv | less -S

mkdir fastqc_output
fastqc data/fastq/yasuda-fastq/2D10_wms_L1_1.fq.gz --outdir fastqc_output
fastqc data/fastq/yasuda-fastq/*.fq.gz --outdir fastqc_output

2. Once I run Kneaddata, I plan to build a Kraken reference database for classification. I would like to include bacteria, archaea, viruses, UniVec_Core, fungi, plants, protozoa, plasmids, and nt. Do you think it's too much, or is something missing? Initially, I thought adding plants would help to identify plant DNA from the host diet. But maybe there won't be much plant DNA. 


seq_R1_kneaddata.log										Log file of the kneaddata run
seq_R1_kneaddata.repeats.removed.1.fastq 					
seq_R1_kneaddata.repeats.removed.2.fastq
seq_R1_kneaddata.repeats.removed.unmatched.1.fastq 
seq_R1_kneaddata.repeats.removed.unmatched.2.fastq
seq_R1_kneaddata.trimmed.1.fastq						This file has trimmed reads (Mate 1) as a output of Paired Ends run in Trimmomatic
seq_R1_kneaddata.trimmed.2.fastq						This file has trimmed reads (Mate 2) as a output of Paired Ends run in Trimmomatic
seq_R1_kneaddata.trimmed.single.1.fastq				This file has trimmed reads (only Mate 1 survived) after running Trimmomatic
seq_R1_kneaddata.trimmed.single.2.fastq				This file has trimmed reads (only Mate 2 survived) after running Trimmomatic
seq_R1_kneaddata_demo_db_bowtie2_paired_contam_1.fastq 	FASTQ file containing reads that were identified as contaminants from the database.
seq_R1_kneaddata_demo_db_bowtie2_paired_contam_2.fastq 		FASTQ file containing reads that were identified as contaminants from the database.
seq_R1_kneaddata_demo_db_bowtie2_unmatched_1_contam.fastq	This file includes reads (Mate 1) that were identified as contaminants from the database
seq_R1_kneaddata_demo_db_bowtie2_unmatched_2_contam.fastq		This file includes reads (Mate 2) that were identified as contaminants from the database
seq_R1_kneaddata_paired_1.fastq			Final output of KneadData after running Trimmomatic + Bowtie2 for seq1
seq_R1_kneaddata_paired_2.fastq			Final output of KneadData after running Trimmomatic + Bowtie2 for seq1
seq_R1_kneaddata_unmatched_1.fastq	Final output of KneadData after running Trimmomatic + Bowtie2 for seq1 (only Mate 1 survived)
seq_R1_kneaddata_unmatched_2.fastq	Final output of KneadData after running Trimmomatic + Bowtie2 for seq1 (only Mate 2 survived)