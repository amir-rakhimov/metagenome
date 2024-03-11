#!/usr/bin/env bash
# Setup channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# Use conda to install all tools
conda create -yqn kraken2-tools -c conda-forge -c bioconda kraken2 krakentools bracken krona
conda install scipy
conda install bowtie2
source ~/miniconda3/etc/profile.d/conda.sh
conda activate kraken2-tools

# Alternatively, install software with binaries
# Create a directory to store all of the executable programs used in this protocol (if none already exists):
mkdir $HOME/bin
# Add the above directory to your PATH environment variable: 
export PATH=$HOME/bin:$PATH

# Install Kraken2
wget https://github.com/DerrickWood/kraken2/archive/refs/tags/v2.1.3.zip
unzip v2.1.3.zip
cd kraken2-2.1.3/
./install_kraken.sh .

# Install Bracken
wget https://github.com/jenniferlu717/Bracken/archive/refs/tags/v2.9.zip
unzip v2.9.zip
cd Bracken-2.9/
sh install_bracken.sh

# Install KrakenTools
wget https://github.com/jenniferlu717/KrakenTools/archive/refs/tags/v1.2.zip
unzip v1.2

# Bowtie2 
wget https://github.com/BenLangmead/bowtie2/releases/download/v2.5.2/bowtie2-2.5.2-linux-x86_64.zip
unzip bowtie2-2.5.2-linux-x86_64.zip
cp bowtie2-2.5.2-linux-x86_64/bowtie2* $HOME/bin

# Copy the kraken executables to a directory in your PATH, 
# and create symlinks for the Bracken and KrakenTools executables

cp kraken2-2.1.3/kraken2 $HOME/bin
cp kraken2-2.1.3/kraken2-build $HOME/bin
cp kraken2-2.1.3/kraken2-inspect $HOME/bin
ln -s ~/Bracken-2.9/bracken $HOME/bin/bracken
ln -s ~/Bracken-2.9/bracken-build $HOME/bin/bracken-build
ln -s ~/KrakenTools-1.2 $HOME/bin/KrakenTools-1.2

# Download Minikraken database with the cleaned eukaryotic pathogen genomes
wget https://genome-idx.s3.amazonaws.com/kraken/minikraken2_v1_8GB_201904.tgz
tar zxvf minikraken2_v1_8GB_201904.tgz
'
If the desired database is not available, we describe here how to create a
custom Kraken 2 database using the kraken2-build script options:
'
'
# First, download the NCBI taxonomy

kraken2-build --db krakendb --download-taxonomy

# Second, download one or more reference libraries. 
kraken2-build --db krakendb --download-library bacteria
kraken2-build --db krakendb --download-library archaea
kraken2-build --db krakendb --download-library viral
kraken2-build --db krakendb --download-library protozoa
kraken2-build --db krakendb --download-library UniVec_Core
'
'
Third, download additional genomes by adding multi-FASTA or single-FASTA files. 
The FASTA sequence headers must include either (1) NCBI accession numbers or (2) 
the text kraken:taxid followed by the taxonomy ID for the genome 
(e.g., >sequence100|kraken:taxid|9606|). 
If this requirement is met, the following commands will add the sequences to the database:
'

'
kraken2-build --db krakendb --add-to-library chr1.fa
kraken2-build --db krakendb --add-to-library chr2.fa

# Finally, build the Kraken 2 database and generate the Bracken database files
kraken2-build --db krakendb --build --threads 8
bracken-build -d krakendb -t 8 -k 35 -l 100
'

# Download the database from https://genome-idx.s3.amazonaws.com/kraken/k2_standard_eupath_20201202.tar.gz
wget https://genome-idx.s3.amazonaws.com/kraken/k2_standard_eupath_20201202.tar.gz

# Create folders for required data
mkdir k2protocol_db
mkdir m_samples
mkdir b_index

# Download SRA samples (not shown)
~/software/sratoolkit.3.0.7-ubuntu64/bin/prefetch \
	--option-file data/sra/kraken-sra.txt --output-directory m_samples/
 mv m_samples/*/*.sra m_samples/
~/software/sratoolkit.3.0.7-ubuntu64/bin/fastq-dump --outdir m_samples/ \
	--gzip --skip-technical  --readids --read-filter pass \
	--dumpbase --split-3 --clip m_samples/*sra

# Inside the b_index folder, download the k2protocol_bowtie2indices.tgz file 
# from http://ccb.jhu.edu/data/kraken2_protocol/ and unpack the files:
cd b_index
tar -xzvf k2protocol_bowtie2indices.tgz

# Unzip fastq files
cd m_samples
for i in *_pass_*.fastq.gz
do 
  gunzip ${i};
done

# 1. Remove host DNA
wget https://genome-idx.s3.amazonaws.com/bt/GRCh38_noalt_as.zip
unzip GRCh38_noalt_as.zip
mkdir bowtie_output
# The provided sample set has already undergone removal of host (human) DNA.
bowtie2 -p 8 -x GRCh38_noalt_as/GRCh38_noalt_as \
  -1 m_samples/SRR14092160_1.fastq \
  -2 m_samples/SRR14092160_2.fastq \
  --very-sensitive-local \
  --un-conc bowtie_output/nonhuman_reads.fastq \
  -S bowtie_output/human_reads.sam
# GRCh38_noalt_as/GRCh38_noalt_as means that we find the index files 
# in the GRCh38_noalt_as directory. The index file names start with GRCh38_noalt_as
# using --un-conc to get fastq output files;  8 processors
# use alignment mode "local"  especially in case raw-reads are not adapter/quality trimmed
'
36510750 reads; of these:
	36510750 (100.00%) were paired; of these:
		36510748 (100.00%) aligned concordantly 0 times
		1 (0.00%) aligned concordantly exactly 1 time
		1 (0.00%) aligned concordantly >1 times
		----
		36510748 pairs aligned concordantly 0 times; of these:
			1 (0.00%) aligned discordantly 1 time
		----
		36510747 pairs aligned 0 times concordantly or discordantly; of these:
			73021494 mates make up the pairs; of these:
				73019511 (100.00%) aligned 0 times
				1688 (0.00%) aligned exactly 1 time
				295 (0.00%) aligned >1 times
0.00% overall alignment rate
'
# 2. Classify microbiome samples using Kraken
mkdir kreports
mkdir kraken_outputs
mkdir breports
mkdir bracken_outputs
' We need the database and FASTA file of sequences
kraken2 --db $DBNAME seqs.fq
'

kraken2 --db minikraken2_v1_8GB \
  --threads 8 \
  --paired \
  --report kreports/SRR14092160.k2report \
  --minimum-hit-groups 3 \
  --report-minimizer-data \
  m_samples/SRR14092160_pass_1.fastq.gz \
  m_samples/SRR14092160_pass_2.fastq.gz \
  > kraken_outputs/SRR14092160.kraken2
kraken2 --db minikraken2_v1_8GB \
  --threads 8 \
  --paired \
  --report kreports/SRR14092310.k2report \
  --minimum-hit-groups 3 \
  --report-minimizer-data \
  m_samples/SRR14092310_pass_1.fastq.gz \
  m_samples/SRR14092310_pass_2.fastq.gz \
  > kraken_outputs/SRR14092310.kraken2
kraken2 --db minikraken2_v1_8GB \
  --threads 8 \
  --paired \
  --report kreports/SRR14143424.k2report \
  --minimum-hit-groups 3 \
  --report-minimizer-data \
  m_samples/SRR14143424_pass_1.fastq.gz \
  m_samples/SRR14143424_pass_2.fastq.gz \
  > kraken_outputs/SRR14143424.kraken2

# Loading database information...   done.
# 60431111 sequences (14461.89 Mbp) processed in 355.423s (10201.6 Kseq/m, 2441.36 Mbp/m).
# 22429192 sequences classified (37.12%)
# 38001919 sequences unclassified (62.88%)
## no minimizer data
kraken2 --db minikraken2_v1_8GB \
  --threads 8 \
  --paired \
  --report kreports/SRR14092160_no_minimizer_data.k2report \
  --minimum-hit-groups 3 \
  m_samples/SRR14092160_pass_1.fastq.gz \
  m_samples/SRR14092160_pass_2.fastq.gz \
  > kraken_outputs/SRR14092160_no_minimizer_data.kraken2
kraken2 --db minikraken2_v1_8GB \
  --threads 8 \
  --paired \
  --report kreports/SRR14092310_no_minimizer_data.k2report \
  --minimum-hit-groups 3 \
  m_samples/SRR14092310_pass_1.fastq.gz \
  m_samples/SRR14092310_pass_2.fastq.gz \
  > kraken_outputs/SRR14092310_no_minimizer_data.kraken2
kraken2 --db minikraken2_v1_8GB \
  --threads 8 \
  --paired \
  --report kreports/SRR14143424_no_minimizer_data.k2report \
  --minimum-hit-groups 3 \
  m_samples/SRR14143424_pass_1.fastq.gz \
  m_samples/SRR14143424_pass_2.fastq.gz \
  > kraken_outputs/SRR14143424_no_minimizer_data.kraken2




# 3. Run bracken for abundance estimation of microbiome samples
' Generic Bracken invocation
bracken -d kraken_database -i sample.k2report -r read_length \
-l taxonomic_level -t read_threshold -o sample.bracken -w sample.breport\
'
bracken -d minikraken2_v1_8GB \
  -i kreports/SRR14143424.k2report -r 100 -l S -t 10 \
  -o bracken_outputs/SRR14143424.bracken -w breports/SRR14143424.breport

bracken -d minikraken2_v1_8GB \
  -i kreports/SRR14092160.k2report -r 100 -l S -t 10 \
  -o bracken_outputs/SRR14092160.bracken -w breports/SRR14092160.breport

bracken -d minikraken2_v1_8GB \
  -i kreports/SRR14092310.k2report -r 100 -l S -t 10 \
  -o bracken_outputs/SRR14092310.bracken -w breports/SRR14092310.breport

bracken -d minikraken2_v1_8GB \
  -i kreports/SRR14143424_no_minimizer_data.k2report -r 100 -l S -t 10 \
  -o bracken_outputs/SRR14143424_no_minimizer_data.bracken \
  -w breports/SRR14143424_no_minimizer_data.breport


bracken -d minikraken2_v1_8GB \
  -i kreports/SRR14143424_no_minimizer_data.k2report -r 100 -l G -t 10 \
  -o bracken_outputs/SRR14143424_no_minimizer_data_G.bracken \
  -w breports/SRR14143424_no_minimizer_data_G.breport

bracken -d minikraken2_v1_8GB \
  -i kreports/SRR14092160_no_minimizer_data.k2report -r 100 -l S -t 10 \
  -o bracken_outputs/SRR14092160_no_minimizer_data.bracken \
  -w breports/SRR14092160_no_minimizer_data.breport

bracken -d minikraken2_v1_8GB \
  -i kreports/SRR14092310_no_minimizer_data.k2report -r 100 -l S -t 10 \
  -o bracken_outputs/SRR14092310_no_minimizer_data.bracken \
  -w breports/SRR14092310_no_minimizer_data.breport

# No difference between these reports
combine_kreports.py -r kreports/SRR14092160.k2report kreports/SRR14092310.k2report \
  kreports/SRR14143424.k2report \
  -o kraken_report_all.txt
combine_kreports.py -r kreports/*_no_minimizer_data.k2report \
  -o kraken_report_all_no_minimizer_data.txt

combine_kreports.py -r breports/SRR14092160.breport breports/SRR14092310.breport breports/SRR14143424.breport -o bracken_report_all.txt
combine_kreports.py -r breports/*_no_minimizer_data.breport \
  -o bracken_report_all_no_minimizer_data.txt

# perc............percentage of total reads rooted at this clade
# tot_all ........total reads rooted at this clade (including reads at more specific clades)
# tot_lvl.........total reads at this clade (not including reads at more specific clades)
# 1_all...........reads from Sample 1 rooted at this clade
# 1_lvl...........reads from Sample 1 at this clade
# 2_all...........""
# 2_lvl...........""
# etc..
# lvl_type........Clade level type (R, D, P, C, O, F, G, S....)
# taxid...........taxonomy ID of this clade
# name............name of this clade


# Kreport2mpa: Only with no minimizer data!
kreport2mpa.py -r kreports/SRR14092160_no_minimizer_data.k2report \
  -o kreports/SRR14092160_no_minimizer_data.mpa.txt
kreport2mpa.py -r kreports/SRR14092310_no_minimizer_data.k2report \
  -o kreports/SRR14092310_no_minimizer_data.mpa.txt
kreport2mpa.py -r kreports/SRR14143424_no_minimizer_data.k2report \
  -o kreports/SRR14143424_no_minimizer_data.mpa.txt

# On bracken also works!
kreport2mpa.py -r breports/SRR14092160_no_minimizer_data.breport \
  -o breports/SRR14092160_breport_no_minimizer_data.mpa.txt
kreport2mpa.py -r breports/SRR14092310_no_minimizer_data.breport \
  -o breports/SRR14092310_breport_no_minimizer_data.mpa.txt
kreport2mpa.py -r breports/SRR14143424_no_minimizer_data.breport \
  -o breports/SRR14143424_breport_no_minimizer_data.mpa.txt


# combine_mpa.py combines multiple outputs from kreport2mpa.py: can also use on bracken
combine_mpa.py -i kreports/*_no_minimizer_data.mpa.txt -o kreports/combined_mpa.txt
combine_mpa.py -i breports/*_no_minimizer_data.mpa.txt -o breports/bracken_combined_mpa.txt


# Calculate alpha-diversity
# with conda:
alpha_diversity.py -f bracken_outputs/SRR14143424.bracken -a BP
alpha_diversity.py -f bracken_outputs/SRR14143424.bracken -a Sh
alpha_diversity.py -f bracken_outputs/SRR14143424.bracken -a F
alpha_diversity.py -f bracken_outputs/SRR14143424.bracken -a Si
alpha_diversity.py -f bracken_outputs/SRR14143424.bracken -a ISi
alpha_diversity.py -f bracken_outputs/SRR14092160.bracken -a BP
alpha_diversity.py -f bracken_outputs/SRR14092160.bracken -a Sh
alpha_diversity.py -f bracken_outputs/SRR14092160.bracken -a F
alpha_diversity.py -f bracken_outputs/SRR14092160.bracken -a Si
alpha_diversity.py -f bracken_outputs/SRR14092160.bracken -a ISi
alpha_diversity.py -f bracken_outputs/SRR14092310.bracken -a BP
alpha_diversity.py -f bracken_outputs/SRR14092310.bracken -a Sh
alpha_diversity.py -f bracken_outputs/SRR14092310.bracken -a F
alpha_diversity.py -f bracken_outputs/SRR14092310.bracken -a Si
alpha_diversity.py -f bracken_outputs/SRR14092310.bracken -a ISi

# Calculate beta-diversity
beta_diversity.py -i bracken_outputs/SRR14092160.bracken \
bracken_outputs/SRR14092310.bracken \
bracken_outputs/SRR14143424.bracken --type bracken

# Generate Krona plots
mkdir b_krona_txt
mkdir krona_html
kreport2krona.py -r breports/SRR14143424.breport \
-o b_krona_txt/SRR14143424.b.krona.txt --no-intermediate-ranks
kreport2krona.py -r breports/SRR14143424_no_minimizer_data.breport \
-o b_krona_txt/SRR14143424_no_minimizer_data.b.krona.txt --no-intermediate-ranks
ktImportText b_krona_txt/SRR14143424.b.krona.txt \
-o krona_html/SRR14143424.krona.html
kreport2krona.py -r breports/SRR14092160.breport \
-o b_krona_txt/SRR14092160.b.krona.txt --no-intermediate-ranks
ktImportText b_krona_txt/SRR14092160.b.krona.txt \
-o krona_html/SRR14092160.krona.html
kreport2krona.py -r breports/SRR14092310.breport \
-o b_krona_txt/SRR14092310.b.krona.txt --no-intermediate-ranks
ktImportText b_krona_txt/SRR14092310.b.krona.txt \
-o krona_html/SRR14092310.krona.html