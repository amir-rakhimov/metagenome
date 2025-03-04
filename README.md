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
* SAM and BAM files with reads that are mapped and unmapped to reference genomes 
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

Note 1: The environment will also be used in the MAG assembly because it contains samtools.
Note 2: All FASTQ files are gzip-compressed
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
bowtie2_output_sam_dir=output/bowtie2_pipeline/bowtie2_output_sam
bowtie2_output_bam_dir=output/bowtie2_pipeline/bowtie2_output_bam
bowtie2_filtered_bam_dir=output/bowtie2_pipeline/bowtie2_filtered_bam
bowtie2_sorted_bam_dir=output/bowtie2_pipeline/bowtie2_sorted_bam
bowtie2_decontam_fastq_dir=data/bowtie2_decontam_fastq
fastqc_output_decontam_dir=output/qc_pipeline/fastqc_output_decontam


1. Renames FASTQ files to make them shorter because the sequencing company gave the files long 
names
2. Downloads reference genomes for decontamination. We are using genomes of
naked mole-rat (Heter_glaber.v1.7_hic_pac), Damaraland mole-rat (DMR_v1.0_HiC), human (T2T),
mouse (GRCm39), and human decoy genome (as specified in https://www.cureffi.org/2013/02/01/the-decoy-genome/). 

The reference genomes are merged into one FASTA file and saved in 
the `~/common_data/reference_genomes` directory 
as `all_hosts_reference.fasta` file. BowTie2 also builds an index using `all_hosts_reference.fasta` 
and saves it in the `~/common_data/bowtie2_indices` directory.

3. Runs FASTQC and MultiQC on raw reads.

Input: FASTQ files in the `data/fastq/yasuda-fastq` directory.

Output: reports saved in the `output/qc_pipeline/fastqc_output`
and `output/qc_pipeline/multiqc_output` directories.
4. Trims raw reads with Cutadapt. The parameters are `--max-n` = 0.1 (do not allow > 10% Ns),
`-q` = 5 (remove bases with basequal < 5), `-O` = 5 (Require MINLENGTH overlap between
read and adapter for an adapter to be found), --discard-trimmed (discard reads in
which an adapter was found), `--minimum-length` = 75 (minimum length after trimming). 

Input: FASTQ files in the `data/fastq/yasuda-fastq` directory.

Output: FASTQ files saved in the `output/qc_pipeline/cutadapt_output` 
directory. The Cutadapt command also saves a report for each sample in the same directory.

5. Aligns trimmed reads to the reference genomes with BowTie2.

Input: the BowTie2 index (`~/common_data/bowtie2_indices/all_hosts_reference`), trimmed
FASTQ files (`output/qc_pipeline/cutadapt_output`).

Output: SAM and BAM files with reads that are mapped and unmapped to reference genomes. The SAM
files for each sample are saved in the `output/bowtie2_pipeline/bowtie2_output_sam` directory as
`${sample_name}_trim_mapped_and_unmapped.sam`. For example, the SAM file for 2D10 sample would be 
`2D10_trim_mapped_and_unmapped.sam`. The BAM files are saved in the
`output/bowtie2_pipeline/bowtie2_output_bam` directory in the same format as SAM files.

6. Filters out read pairs that did not map to reference genomes. These are reads free of 
contamination.

Input: BAM files with mapped and unmapped reads from the previous step
(`output/bowtie2_pipeline/bowtie2_output_bam`).

Output: BAM files where both reads did not map to reference genomes. They are saved in the
`output/bowtie2_pipeline/bowtie2_filtered_bam` directory and are named like 
`${sample_name}_bothReadsUnmapped.bam`.

7. Sorts the BAM file by read name to have paired reads next to each other. Then, it splits
paired-end reads into separated fastq files: .._R1 .._R2. 

Input: BAM files where both reads did not map to reference genomes 
(`output/bowtie2_pipeline/bowtie2_filtered_bam/${sample_name}_bothReadsUnmapped.bam`).

Output: Sorted BAM file (`${base_name}_bothReadsUnmapped_sorted.bam`) in the 
`output/bowtie2_pipeline/bowtie2_sorted_bam` directory; decontaminated reads 
(`${sample_name}_decontam_R1.fastq.gz` and `${sample_name}_decontam_R2.fastq.gz`) in the 
`data/bowtie2_decontam_fastq` directory.

8. Runs FastQC and MultiQC on decontaminated data as a final check.

Input: decontaminated reads in the `data/bowtie2_decontam_fastq` directory. 

Output reports are saved in the 
`output/qc_pipeline/fastqc_output` and `output/qc_pipeline/multiqc_output` directories.


# Required packages for the conda environment
Name                    Version                   Build    Channel
_libgcc_mutex             0.1                 conda_forge    conda-forge
_openmp_mutex             4.5                       2_gnu    conda-forge
alsa-lib                  1.2.11               hd590300_1    conda-forge
annotated-types           0.7.0              pyhd8ed1ab_0    conda-forge
bbmap                     39.06                h92535d8_1    bioconda
blast                     2.15.0          pl5321h6f7f691_1    bioconda
bowtie2                   2.5.4                he20e202_0    bioconda
brotli-python             1.1.0           py311hb755f60_1    conda-forge
bzip2                     1.0.8                hd590300_5    conda-forge
c-ares                    1.28.1               hd590300_0    conda-forge
ca-certificates           2025.1.31            hbcca054_0    conda-forge
cairo                     1.18.0               h3faef2a_0    conda-forge
certifi                   2025.1.31          pyhd8ed1ab_0    conda-forge
cffi                      1.16.0          py311hb3a22ac_0    conda-forge
charset-normalizer        3.3.2              pyhd8ed1ab_0    conda-forge
click                     8.1.7           unix_pyh707e725_0    conda-forge
colorama                  0.4.6              pyhd8ed1ab_0    conda-forge
coloredlogs               15.0.1             pyhd8ed1ab_3    conda-forge
colormath                 3.0.0                      py_2    conda-forge
csvtk                     0.32.0               h9ec0b34_0    conda-forge
curl                      8.8.0                he654da7_0    conda-forge
cutadapt                  4.8             py311hdad781d_1    bioconda
dnaio                     1.2.0           py311hdad781d_2    bioconda
entrez-direct             22.1                 he881be0_0    bioconda
expat                     2.6.2                h59595ed_0    conda-forge
fastqc                    0.12.1               hdfd78af_0    bioconda
font-ttf-dejavu-sans-mono 2.37                 hab24e00_0    conda-forge
font-ttf-inconsolata      3.000                h77eed37_0    conda-forge
font-ttf-source-code-pro  2.038                h77eed37_0    conda-forge
font-ttf-ubuntu           0.83                 h77eed37_2    conda-forge
fontconfig                2.14.2               h14ed4e7_0    conda-forge
fonts-conda-ecosystem     1                             0    conda-forge
fonts-conda-forge         1                             0    conda-forge
freetype                  2.12.1               h267a509_2    conda-forge
gettext                   0.22.5               h59595ed_2    conda-forge
gettext-tools             0.22.5               h59595ed_2    conda-forge
giflib                    5.2.2                hd590300_0    conda-forge
graphite2                 1.3.13            h59595ed_1003    conda-forge
harfbuzz                  8.5.0                hfac3d4d_0    conda-forge
htslib                    1.20                 h5efdd21_1    bioconda
humanfriendly             10.0               pyhd8ed1ab_6    conda-forge
humanize                  4.9.0              pyhd8ed1ab_0    conda-forge
icu                       73.2                 h59595ed_0    conda-forge
idna                      3.7                pyhd8ed1ab_0    conda-forge
importlib-metadata        7.1.0              pyha770c72_0    conda-forge
importlib_metadata        7.1.0                hd8ed1ab_0    conda-forge
isa-l                     2.31.0               hd590300_1    conda-forge
jinja2                    3.1.4              pyhd8ed1ab_0    conda-forge
kaleido-core              0.2.1                h3644ca4_0    conda-forge
keyutils                  1.6.1                h166bdaf_0    conda-forge
krb5                      1.21.2               h659d440_0    conda-forge
lcms2                     2.16                 hb7c19ff_0    conda-forge
ld_impl_linux-64          2.40                 hf3520f5_2    conda-forge
lerc                      4.0.0                h27087fc_0    conda-forge
libasprintf               0.22.5               h661eb56_2    conda-forge
libasprintf-devel         0.22.5               h661eb56_2    conda-forge
libblas                   3.9.0           22_linux64_openblas    conda-forge
libcblas                  3.9.0           22_linux64_openblas    conda-forge
libcups                   2.3.3                h4637d8d_4    conda-forge
libcurl                   8.8.0                hca28451_0    conda-forge
libdeflate                1.20                 hd590300_0    conda-forge
libedit                   3.1.20191231         he28a2e2_2    conda-forge
libev                     4.33                 hd590300_2    conda-forge
libexpat                  2.6.2                h59595ed_0    conda-forge
libffi                    3.4.2                h7f98852_5    conda-forge
libgcc                    14.2.0               h77fa898_1    conda-forge
libgcc-ng                 14.2.0               h69a702a_1    conda-forge
libgettextpo              0.22.5               h59595ed_2    conda-forge
libgettextpo-devel        0.22.5               h59595ed_2    conda-forge
libgfortran-ng            13.2.0               h69a702a_7    conda-forge
libgfortran5              13.2.0               hca663fb_7    conda-forge
libglib                   2.80.2               hf974151_0    conda-forge
libgomp                   14.2.0               h77fa898_1    conda-forge
libhwloc                  2.10.0          default_h5622ce7_1001    conda-forge
libiconv                  1.17                 hd590300_2    conda-forge
libidn2                   2.3.7                hd590300_0    conda-forge
libjpeg-turbo             3.0.0                hd590300_1    conda-forge
liblapack                 3.9.0           22_linux64_openblas    conda-forge
libnghttp2                1.58.0               h47da74e_1    conda-forge
libnsl                    2.0.1                hd590300_0    conda-forge
libopenblas               0.3.27          pthreads_h413a1c8_0    conda-forge
libpng                    1.6.43               h2797004_0    conda-forge
libsqlite                 3.45.3               h2797004_0    conda-forge
libssh2                   1.11.0               h0841786_0    conda-forge
libstdcxx-ng              13.2.0               hc0a3c3a_7    conda-forge
libtiff                   4.6.0                h1dd3fc0_3    conda-forge
libunistring              0.9.10               h7f98852_0    conda-forge
libuuid                   2.38.1               h0b41bf4_0    conda-forge
libwebp-base              1.4.0                hd590300_0    conda-forge
libxcb                    1.15                 h0b41bf4_0    conda-forge
libxcrypt                 4.4.36               hd590300_1    conda-forge
libxml2                   2.12.7               hc051c1a_1    conda-forge
libzlib                   1.2.13               h4ab18f5_6    conda-forge
markdown                  3.6                pyhd8ed1ab_0    conda-forge
markdown-it-py            3.0.0              pyhd8ed1ab_0    conda-forge
markupsafe                2.1.5           py311h459d7ec_0    conda-forge
mathjax                   2.7.7                ha770c72_3    conda-forge
mdurl                     0.1.2              pyhd8ed1ab_0    conda-forge
multiqc                   1.22.2             pyhdfd78af_0    bioconda
ncbi-vdb                  3.1.1                h4ac6f70_0    bioconda
ncurses                   6.5                  h59595ed_0    conda-forge
networkx                  3.3                pyhd8ed1ab_1    conda-forge
nspr                      4.35                 h27087fc_0    conda-forge
nss                       3.100                hca3bf56_0    conda-forge
numpy                     1.26.4          py311h64a7726_0    conda-forge
openjdk                   22.0.1               hb622114_0    conda-forge
openjpeg                  2.5.2                h488ebb8_0    conda-forge
openssl                   3.4.1                h7b32b05_0    conda-forge
packaging                 24.0               pyhd8ed1ab_0    conda-forge
pbzip2                    1.1.13               h1fcc475_2    conda-forge
pcre                      8.45                 h9c3ff4c_0    conda-forge
pcre2                     10.43                hcad00b1_0    conda-forge
perl                      5.32.1          7_hd590300_perl5    conda-forge
perl-archive-tar          2.40            pl5321hdfd78af_0    bioconda
perl-carp                 1.50            pl5321hd8ed1ab_0    conda-forge
perl-common-sense         3.75            pl5321hd8ed1ab_0    conda-forge
perl-compress-raw-bzip2   2.201           pl5321h166bdaf_0    conda-forge
perl-compress-raw-zlib    2.202           pl5321h166bdaf_0    conda-forge
perl-encode               3.21            pl5321hd590300_0    conda-forge
perl-exporter             5.74            pl5321hd8ed1ab_0    conda-forge
perl-exporter-tiny        1.002002        pl5321hd8ed1ab_0    conda-forge
perl-extutils-makemaker   7.70            pl5321hd8ed1ab_0    conda-forge
perl-io-compress          2.201           pl5321hdbdd923_2    bioconda
perl-io-zlib              1.14            pl5321hdfd78af_0    bioconda
perl-json                 4.10            pl5321hdfd78af_0    bioconda
perl-json-xs              2.34            pl5321h4ac6f70_6    bioconda
perl-list-moreutils       0.430           pl5321hdfd78af_0    bioconda
perl-list-moreutils-xs    0.430           pl5321h031d066_2    bioconda
perl-parent               0.241           pl5321hd8ed1ab_0    conda-forge
perl-pathtools            3.75            pl5321h166bdaf_0    conda-forge
perl-scalar-list-utils    1.63            pl5321h166bdaf_0    conda-forge
perl-storable             3.15            pl5321h166bdaf_0    conda-forge
perl-types-serialiser     1.01            pl5321hdfd78af_0    bioconda
pigz                      2.8                  h2797004_0    conda-forge
pillow                    10.3.0          py311h18e6fac_0    conda-forge
pip                       24.0               pyhd8ed1ab_0    conda-forge
pixman                    0.43.2               h59595ed_0    conda-forge
plotly                    5.22.0             pyhd8ed1ab_0    conda-forge
pthread-stubs             0.4               h36c2ea0_1001    conda-forge
pyaml-env                 1.1.0              pyhd8ed1ab_0    conda-forge
pycparser                 2.22               pyhd8ed1ab_0    conda-forge
pydantic                  2.7.3              pyhd8ed1ab_0    conda-forge
pydantic-core             2.18.4          py311h5ecf98a_0    conda-forge
pygments                  2.18.0             pyhd8ed1ab_0    conda-forge
pysocks                   1.7.1              pyha2e5f31_6    conda-forge
python                    3.11.9          hb806964_0_cpython    conda-forge
python-isal               1.6.1           py311h459d7ec_0    conda-forge
python-kaleido            0.2.1              pyhd8ed1ab_0    conda-forge
python-zlib-ng            0.4.3           py311hb3d061c_0    conda-forge
python_abi                3.11                    4_cp311    conda-forge
pyyaml                    5.4.1           py311hd4cff14_4    conda-forge
readline                  8.2                  h8228510_1    conda-forge
requests                  2.32.3             pyhd8ed1ab_0    conda-forge
rich                      13.7.1             pyhd8ed1ab_0    conda-forge
rich-click                1.8.2              pyhd8ed1ab_0    conda-forge
samtools                  1.20                 h50ea8bc_0    bioconda
seqkit                    2.8.2                h9ee0642_0    bioconda
seqtk                     1.4                  he4a0461_2    bioconda
setuptools                70.0.0             pyhd8ed1ab_0    conda-forge
spectra                   0.0.11                     py_1    conda-forge
sqlite                    3.45.3               h2c6b66d_0    conda-forge
tbb                       2021.12.0            h297d8ca_1    conda-forge
tenacity                  8.3.0              pyhd8ed1ab_0    conda-forge
tk                        8.6.13          noxft_h4845f30_101    conda-forge
tqdm                      4.66.4             pyhd8ed1ab_0    conda-forge
trf                       4.09.1               h031d066_5    bioconda
typeguard                 4.3.0              pyhd8ed1ab_0    conda-forge
typing-extensions         4.12.1               hd8ed1ab_0    conda-forge
typing_extensions         4.12.1             pyha770c72_0    conda-forge
tzdata                    2024a                h0c530f3_0    conda-forge
urllib3                   2.2.1              pyhd8ed1ab_0    conda-forge
wget                      1.21.4               hda4d442_0    conda-forge
wheel                     0.43.0             pyhd8ed1ab_1    conda-forge
xopen                     2.0.1           py311h38be061_0    conda-forge
xorg-fixesproto           5.0               h7f98852_1002    conda-forge
xorg-inputproto           2.3.2             h7f98852_1002    conda-forge
xorg-kbproto              1.0.7             h7f98852_1002    conda-forge
xorg-libice               1.1.1                hd590300_0    conda-forge
xorg-libsm                1.2.4                h7391055_0    conda-forge
xorg-libx11               1.8.9                h8ee46fc_0    conda-forge
xorg-libxau               1.0.11               hd590300_0    conda-forge
xorg-libxdmcp             1.1.3                h7f98852_0    conda-forge
xorg-libxext              1.3.4                h0b41bf4_2    conda-forge
xorg-libxfixes            5.0.3             h7f98852_1004    conda-forge
xorg-libxi                1.7.10               h7f98852_0    conda-forge
xorg-libxrender           0.9.11               hd590300_0    conda-forge
xorg-libxt                1.3.0                hd590300_1    conda-forge
xorg-libxtst              1.2.3             h7f98852_1002    conda-forge
xorg-recordproto          1.14.2            h7f98852_1002    conda-forge
xorg-renderproto          0.11.1            h7f98852_1002    conda-forge
xorg-xextproto            7.3.0             h0b41bf4_1003    conda-forge
xorg-xproto               7.0.31            h7f98852_1007    conda-forge
xz                        5.2.6                h166bdaf_0    conda-forge
yaml                      0.2.5                h7f98852_2    conda-forge
zipp                      3.17.0             pyhd8ed1ab_0    conda-forge
zlib                      1.2.13               h4ab18f5_6    conda-forge
zlib-ng                   2.0.7                h0b41bf4_0    conda-forge
zstandard                 0.19.0          py311hd4cff14_0    conda-forge
zstd                      1.5.6                ha6fb4c9_0    conda-forge