# Whole metagenome analysis of fecal microbiota with Kraken2 and MegaHIT-MetaBAT2 pipeline

This project solves three problems (each problems is a pipeline):
1. Quality check and decontamination of whole metagenome sequencing reads with BowTie2
2. Taxonomic analysis of decontaminated whole metagenome sequencing reads with Kraken2
and Bracken
3. Assembly of decontaminated whole metagenome sequencing reads into 
metagenome-assembled genomes (MAGs) with MegaHIT and MetaBAT2 
(there are other tools used in the pipeline, but these two are the main ones). 
MegaHIT assembles reads into contigs, and MetaBAT2 groups contigs into bins (i.e. MAGs).

* The required packages and their version are listed at the end of this README.md file
# TODO: * Time required to run each pipeline (11 samples):
1. X hours
2. Y hours
3. Z hours

# Pipeline script, input, output, and example files
## 1. Quality check and decontamination with BowTie2
### 1.1 Script
The script of the quality check pipeline (pipeline #1): `qc-commands-with-bowtie2.sh` 
# TODO: (check the ___.md file for step-by-step explanantion)
### 1.2 Input
FASTQ files (gzip-compressed) with raw paired-end demultiplexed 
reads (with or without adapters). If adapters are not removed, please specify
 adapter sequences in the `qc-commands-with-bowtie2.sh` script
### 1.3 Output 
* Quality check reports from FASTQC and MultiQC on raw data
* FASTQ files with trimmed reads (gzip-compressed)
* SAM and BAM files with reads that are mapped and unmapped to reference genomes. 
* BAM file with read pairs that are unmapped (both R1 and R2 are unmapped) to the reference
genomes
* Sorted BAM file with read pairs that are unmapped (sorted by read name)
* Final FASTQ files with decontaminated reads: only read pairs that were unmapped to 
reference genomes are retained in the final files
* Final quality check reports from FASTQC and MultiQC

# TODO ### 1.4 Example files 
List them

## 2. Taxonomic analysis with Kraken2 and Bracken
### 2.1 Script
The script of the quality check pipeline (pipeline #1): `kraken2-pipeline.sh` 
### 2.2 Input
FASTQ files with decontaminated reads (either from pipeline #1 or supplied by user). Please 
adjust Kraken2 and Bracken parameters like kmer length, minimizer length, etc. according to your needs
### 2.3 Output
* Kraken2 output files (.kraken2) without minimizers. However, it's possible to add minimizers if needed
* Kraken2 reports (.k2report) with taxonomic classification
* FASTQ files with classified sequences (read 1 and read 2 are in separate files), gzip-compressed
* Bracken output files (.bracken)
* Bracken reports (.breport)
* TXT files for Krona generated from Bracken reports (one TXT file per sample)
* HTML files with Krona plots generated from TXT files (one HTML file per sample)
* FASTA files with reads that were mapped to Cryptomeria japonica (taxonomy id 3369),
Glycine soja (taxid 3848), Ipomoea triloba (taxid 35885), and Macadamia integrifolia (taxid 60698).
* TXT files with taxonomic classification from BLAST on selected species (using FASTA files)
### 2.4 Example files

## 3. MAG assembly with MegaHIT and MetaBAT2
### 3.1 Script
The script of the quality check pipeline (pipeline #1): `kraken2-pipeline.sh` 
### 3.2 Input
FASTQ files with decontaminated reads (either from pipeline #1 or supplied by user).
### 3.3 Output
SAM and BAM files are gzip-compressed to save disk space
* FASTA files with contigs from MEGAHIT
* Quality check reports on contigs from MetaQUAST
* FAA files with predicted genes from Prodigal
* GBK files from Prodigal
* SAM files with reads that are mapped and unmapped to contigs
* Coverage information for each contig from Pileup.sh (BBMap)
* FASTQ files with unmapped reads (single read unmapped)
* FASTQ files with unmapped read pairs (both reads unmapped)
* BAM files with reads that mapped and unmapped to contigs
* BAM files with read pairs mapped to contigs (both reads)
* BAM files with read pairs mapped to contigs and sorted by coordinate (both reads)
* BAM files with read pairs unmapped to contigs (both reads)
* BAM files with read pairs unmapped to contigs and sorted by read name
* FASTQ files with read pairs unmapped to contigs and sorted by read name (each pair is split 
into R1 and R2 files)
# TODO: * Contig depth files
* Bins from MetaBAT2 as FASTA files saved in the directory that corresponds to their sample.
* Quality check reports on bins from MetaQUAST
* TXT files with BUSCO reports on each bin
# TODO * CheckM2 reports
# TODO * Taxonomic annotation with GTDB-tk
### 3.4 Example files

# Scripts in detail
`bash-scripts/qc-tools-create-env.sh`:
Creates a conda environment for QC procedures and decontamination. 

Note: The environment will also be used in the MAG assembly because it contains samtools.

`bash-scripts/qc-commands-with-bowtie2.sh`:
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

1. Renames FASTQ files to make them shorter because the sequencing company gave the files long 
names
2. Downloads reference genomes for decontamination. We are using genomes of
naked mole-rat (Heter_glaber.v1.7_hic_pac), Damaraland mole-rat (DMR_v1.0_HiC), human (T2T),
mouse (GRCm39), and human decoy genome (as specified in https://www.cureffi.org/2013/02/01/the-decoy-genome/). 

The reference genomes are merged into one FASTA file and saved in 
the `~/common_data/reference_genomes` directory 
as `all_hosts_reference.fasta` file. BowTie2 also builds an index using `all_hosts_reference.fasta` 
and saves it in the `~/common_data/bowtie2_indices` directory.
3. Runs FASTQC and MultiQC on raw reads
4. Trims raw reads with Cutadapt. The parameters are `--max-n` = 0.1 (do not allow > 10% Ns),
`-q` = 5 (remove bases with basequal < 5), `-O` = 5 (Require MINLENGTH overlap between
read and adapter for an adapter to be found), --discard-trimmed (discard reads in
which an adapter was found), `--minimum-length` = 75 (minimum length after trimming). 

The output FASTQ files are saved in the `output/qc_pipeline/cutadapt_output` 
directory. The Cutadapt command also saves a report for each sample in the same directory.
5. 

Prepares a BowTie2 index
# TODO 







# TODO: Required packages
List them 