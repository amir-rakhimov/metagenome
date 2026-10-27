#!/usr/bin/env bash
############################################################################
# 1. Project root and run IDs
############################################################################
project_home_dir=$HOME/projects/metagenome

# --- QC ---
qc_RUN_ID=$(date +%Y%m%d_%H%M%S)

# --- Kraken2 ---
kraken2_RUN_ID=20240409_17_32_40
kraken2_RUN_ID=$(date +%Y%m%d_%H%M%S)
kraken2_OUTDIR="${project_home_dir}"/results/1-kraken2/"${kraken2_RUN_ID}"
kraken2_LOGDIR="${kraken2_OUTDIR}"/logs
kraken2_TMPDIR="${kraken2_OUTDIR}"/tmp

# --- MAG ---
mag_RUN_ID=20250303_17_54_26 #20250417_23_21_29
# mag_RUN_ID=$(date +%Y%m%d_%H%M%S)
mag_OUTDIR="${project_home_dir}"/results/2-mag/"${mag_RUN_ID}"
mag_LOGDIR="${mag_OUTDIR}"/logs
mag_TMPDIR="${mag_OUTDIR}"/tmp

############################################################################
# 2. Sample Metadata
############################################################################
METADATA_DIR="${project_home_dir}"/data/metadata/tanes #???
SAMPLES_IDS="${METADATA_DIR}"/sample_ids.tsv
SAMPLES_TSV="${METADATA_DIR}"/sample-paths-decontam.tsv
SAMPLES_PATHS_RAW="${METADATA_DIR}"/sample-paths-raw.tsv

############################################################################
# 3. External tools and reference data
############################################################################
sratoolkit_scripts_dir=$HOME/sratoolkit.3.3.0-ubuntu64/bin
metaquast_script_dir=$HOME/quast-5.2.0
java_path=$HOME/jdk-24.0.1/bin/java
picard_jar_path=$HOME/picard.jar
reference_genomes_dir=$HOME/common_data/reference_genomes
bowtie2_indices_dir=$HOME/common_data/bowtie2_indices
bowtie2_index_name=all_hosts_reference_my

############################################################################
# 4. Adapters
############################################################################
FWD_ADAPTER=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT
REV_ADAPTER=GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG

############################################################################
# 5. Compute resouces (threads, memory)
############################################################################
# --- QC ---
# nthreads_qc=40
nthreads_qc=12
nthreads_qc_sort=8
mem_req_qc=8G
mem_req_qc_sort=4G

# --- Kraken2 ---
nthreads_kraken2=49
nthreads_kraken2_build=38

# --- MAG assembly ---
nthreads_megahit=40
nthreads_megahit_sort=35
mem_req_megahit=8G
mem_req_megahit_sort=4G

# --- MAG annotation ---
nthreads_gtdbtk=15
nthreads_pplacer=30
nthreads_comparem=40
nthreads_annotation=40
nthreads_prokka=35
nthreads_coverm=40
nthreads_blast=20

############################################################################
# 6. Tool parameters
############################################################################
# --- Kraken2 / Bracken database build ---
# kmer length for kraken2-build and braken-build
# kmer_len=35 
# minimizer_len=31
# minimizer_spaces=7
### kmer length for kraken2-build and braken-build
kmer_len=30
minimizer_len=26
minimizer_spaces=6
# Rules:
    ### 1. minimizer length l must be no more than 31 for nucleotide databases, and 15 for protein databases
    ### 2. minimizer length l must be no more than the k-mer length
    ### 3. minimizer_space is s<l/4
    ### 4. default kmer=35, l=31, s=7
### read length for kraken2-build and bracken:
read_len=150
# read length for kraken2-build and bracken
read_len=150

# --- Kraken2 database location ---
# kraken2_db_date=20241209
kraken2_db_date=20250919
# kraken2_db_date=$(date +%Y%m%d)
# kraken2_db_dir="${project_home_dir}"/data/kraken2_db/k2_large_20240307
kraken2_db_dir="${project_home_dir}"/data/kraken2_db/k2_large_${kraken2_db_date} 

# --- Kraken2 / Bracken run parameters ---
### confidence for kraken2 (default: 0):
confidence_kraken2=0.1
### threshold for bracken (default: 0)
bracken_threshold=10

############################################################################
# 7. Directories: Raw data, QC, Decontamination
############################################################################
# --- Input ---
sra_files_dir="${project_home_dir}"/data/sra-files
fastq_dir="${project_home_dir}"/data/0-raw/fastq

# --- QC (FastQC, MultiQC, cutadapt) ---
qc_dir="${project_home_dir}"/data/1-qc
fastqc_raw_output_dir="${qc_dir}"/fastqc_raw/"${qc_RUN_ID}"
multiqc_raw_output_dir="${qc_dir}"/multiqc_raw/"${qc_RUN_ID}"
fastqc_trimmed_output_dir="${qc_dir}"/fastqc_trimmed/"${qc_RUN_ID}"
multiqc_trimmed_output_dir="${qc_dir}"/multiqc_trimmed/"${qc_RUN_ID}"
cutadapt_output_dir="${qc_dir}"/cutadapt/"${qc_RUN_ID}"

# --- Decontamination (Bowtie2) ---
bowtie2_dir="${project_home_dir}"/data/2-bowtie2
bowtie2_sam_dir="${bowtie2_dir}"/sam
bowtie2_bam_dir="${bowtie2_dir}"/bam
bowtie2_filtered_bam_dir="${bowtie2_dir}"/filtered
bowtie2_sorted_bam_dir="${bowtie2_dir}"/sorted
fastqc_decontam_dir="${bowtie2_dir}"/fastqc_decontam
multiqc_decontam_dir="${bowtie2_dir}"/multiqc_decontam
bowtie2_decontam_fastq_dir="${bowtie2_dir}"/decontam_fastq # decontaminated fastq files

############################################################################
# 8. Directories: Kraken2 / Bracken / Krona
############################################################################
kraken2_output_dir="${kraken2_OUTDIR}"/01-kraken2_output
kraken2_reports_dir="${kraken2_OUTDIR}"/02-kraken2_reports
kraken2_classified_reads_dir="${kraken2_OUTDIR}"/03-kraken2_classified_reads
bracken_output_dir="${kraken2_OUTDIR}"/04-bracken_output
bracken_reports_dir="${kraken2_OUTDIR}"/05-bracken_reports
bracken_krona_txt_dir="${kraken2_OUTDIR}"/06-bracken_krona_txt
krona_html_dir="${kraken2_OUTDIR}"/07-krona_html

############################################################################
# 9. MAG assembly
############################################################################
# --- 01 Assembly (MEGAHIT) ---
megahit_output_dir="${mag_OUTDIR}"/01-megahit/assembled_contigs
megahit_logs_dir="${mag_OUTDIR}"/01-megahit/logs
megahit_aligned_reads_dir="${mag_OUTDIR}"/01-megahit/aligned_reads
contig_coverages_dir="${mag_OUTDIR}"/01-megahit/contig_coverages
metaquast_megahit_output_dir="${mag_OUTDIR}"/01-megahit/metaquast
bbwrap_refs_dir="${mag_OUTDIR}"/01-megahit/bbwrap_refs

# --- 02 Binning (MetaBAT2) ---
bam_contig_depths_dir="${mag_OUTDIR}"/02-metabat2/bam_contig_depths
metabat2_output_dir="${mag_OUTDIR}"/02-metabat2/bins
metabat2_logs_dir="${mag_OUTDIR}"/02-metabat2/logs
metaquast_metabat2_output_dir="${mag_OUTDIR}"/02-metabat2/metaquast

# --- 03 Bin QC (CheckM2, dRep) ---
checkm2_output_dir="${mag_OUTDIR}"/03-qc/checkm2
drep_output_dir="${mag_OUTDIR}"/03-qc/drep
mag_qc_logs_dir="${mag_OUTDIR}"/03-qc/logs

# --- 04 Taxonomy (GTDB-tk, CompareM) ---
gtdbtk_output_dir="${mag_OUTDIR}"/04-gtdbtk
comparem_output_dir="${mag_OUTDIR}"/04-gtdbtk/comparem

# --- 05 Annotation ---
prodigal_output_dir="${mag_OUTDIR}"/05-annotation/prodigal
prokka_output_dir="${mag_OUTDIR}"/05-annotation/prokka
prokka_logs_dir="${mag_OUTDIR}"/05-annotation/prokka/logs
mmseqs_easy_cluster_output_dir="${mag_OUTDIR}"/05-annotation/mmseqs2/easy_cluster
dbcan_output_dir="${mag_OUTDIR}"/05-annotation/dbcan
dbcan_logs_dir="${mag_OUTDIR}"/05-annotation/dbcan/logs
kofam_scan_output_dir="${mag_OUTDIR}"/05-annotation/kofam_scan
annotation_fasta_dir="${mag_OUTDIR}"/05-annotation/fasta

# --- 06 Gene quantification (Picard, HTseq, TPM) ---
picard_output_dir="${mag_OUTDIR}"/06-gene_quantification/picard
picard_logs_dir="${mag_OUTDIR}"/06-gene_quantification/picard/logs
htseq_gene_count_dir="${mag_OUTDIR}"/06-gene_quantification/htseq
htseq_logs_dir="${mag_OUTDIR}"/06-gene_quantification/htseq/logs
gene_lengths_dir="${mag_OUTDIR}"/06-gene_quantification/gene_lengths
tpm_files_dir="${mag_OUTDIR}"/06-gene_quantification/tpm
prokka_gtf_dir="${mag_OUTDIR}"/05-annotation/prokka/gtf
prokka_contig_gene_maps_dir="${mag_OUTDIR}"/05-annotation/prokka/contig_gene_maps

# --- 07 MAG stats (CoverM, seqkit) ---
coverm_output_dir="${mag_OUTDIR}"/07-mag_stats/coverm
coverm_logs_dir="${mag_OUTDIR}"/07-mag_stats/coverm/logs
coverm_drep_95_dir="${mag_OUTDIR}"/07-mag_stats/coverm/drep_95
coverm_drep_99_dir="${mag_OUTDIR}"/07-mag_stats/coverm/drep_99
coverm_metabat2_dir="${mag_OUTDIR}"/07-mag_stats/coverm/metabat2
coverm_megahit_dir="${mag_OUTDIR}"/07-mag_stats/coverm/megahit

seqkit_output_dir="${mag_OUTDIR}"/07-mag_stats/seqkit
seqkit_megahit_dir="${mag_OUTDIR}"/07-mag_stats/seqkit/megahit
seqkit_drep_95_dir="${mag_OUTDIR}"/07-mag_stats/seqkit/drep_95
seqkit_metabat2_dir="${mag_OUTDIR}"/07-mag_stats/seqkit/metabat2

# --- 08 Contig BLAST (SILVA) ---
blastn_output_dir="${mag_OUTDIR}"/08-contig_blast

# --- 09 CAZyme BLAST ---
blastp_output_dir="${mag_OUTDIR}"/09-cazyme_blast/blastp_results
blastp_logs_dir="${mag_OUTDIR}"/09-cazyme_blast/logs

# --- 09 MetaThermo ---
metathermo_out_dir="${mag_OUTDIR}"/10-meta_thermo

############################################################################
# 10. Databases (SILVA, CAZy)
############################################################################
silva_db_dir="${project_home_dir}"/data/4-mag/silva_db

blastn_db_path="${project_home_dir}"/data/4-mag/blast_dbs/SILVA_138_2_db
blastp_db_path="${project_home_dir}"/data/4-mag/blast_dbs/cazy_pep_db

cazy_data_dir="${project_home_dir}"/data/4-mag/cazy_data
cazy_db_taxonomy_dir="${project_home_dir}"/data/4-mag/cazy_data/taxonomy
cazy_db_fasta_dir="${project_home_dir}"/data/4-mag/cazy_data/fasta