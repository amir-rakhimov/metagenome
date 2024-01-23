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
kraken2 --db k2protocol_db --threads 8 --report kreports/SRR14143424.k2report \
  --report-minimizer-data --minimum-hit-groups 3 m_samples/SRR14143424_1.fastq \
  m_samples/SRR14143424_2.fastq > kraken_outputs/SRR14143424.kraken2
'
Loading database information... done.
60431111 sequences (7148.66 Mbp) processed in 187.147s (19374.4 Kseq/m, 2291.88 Mbp/m).
17964649 sequences classified (29.73%)
42466462 sequences unclassified (70.27%)
'
kraken2 --db k2protocol_db --threads 8 --report kreports/SRR14092160.k2report \
  --report-minimizer-data --minimum-hit-groups 3 m_samples/SRR14092160_1.fastq \
  m_samples/SRR14092160_2.fastq > kraken_outputs/SRR14092160.kraken2
'
Loading database information... done.
73021500 sequences (8907.00 Mbp) processed in 236.398s (18533.5 Kseq/m, 2260.67 Mbp/m).
11996432 sequences classified (16.43%)
61025068 sequences unclassified (83.57%) 
'
kraken2 --db k2protocol_db --threads 8 --report kreports/SRR14092310.k2report \
  --report-minimizer-data --minimum-hit-groups 3 m_samples/SRR14092310_1.fastq \
  m_samples/SRR14092310_2.fastq > kraken_outputs/SRR14092310.kraken2
'
Loading database information... done.
53998936 sequences (6634.00 Mbp) processed in 156.566s (20693.7 Kseq/m, 2542.31 Mbp/m).
36170969 sequences classified (66.98%)
17827967 sequences unclassified (33.02%)
'
# 3. Run bracken for abundance estimation of microbiome samples
' Generic Bracken invocation
bracken -d kraken_database -i sample.k2report -r read_length \
-l taxonomic_level -t read_threshold -o sample.bracken -w sample.breport\
'
bracken -d k2protocol_db -i kreports/SRR14143424.k2report -r 100 -l S -t 10 \
  -o bracken_outputs/SRR14143424.bracken -w breports/SRR14143424.breport
'
BRACKEN SUMMARY (Kraken report: kreports/SRR14143424.k2report)
>>> Threshold: 10
>>> Number of species in sample: 1570
	>> Number of species with reads > threshold: 500
	>> Number of species with reads < threshold: 1070
>>> Total reads in sample: 60431111
	>> Total reads kept at species level (reads > threshold): 15283051
	>> Total reads discarded (species reads < threshold): 2787
	>> Reads distributed: 2677396
	>> Reads not distributed (eg. no species above threshold): 1415
	>> Unclassified reads: 42466462
'

bracken -d k2protocol_db -i kreports/SRR14092160.k2report -r 100 -l S -t 10 \
  -o bracken_outputs/SRR14092160.bracken -w breports/SRR14092160.breport
'
>>> Threshold: 10
>>> Number of species in sample: 1722
	>> Number of species with reads > threshold: 554
	>> Number of species with reads < threshold: 1168
>>> Total reads in sample: 73021500
	>> Total reads kept at species level (reads > threshold): 9277455
	>> Total reads discarded (species reads < threshold): 3165
	>> Reads distributed: 2714446
	>> Reads not distributed (eg. no species above threshold): 1366
	>> Unclassified reads: 61025068
'
bracken -d k2protocol_db -i kreports/SRR14092310.k2report -r 100 -l S -t 10 \
-o bracken_outputs/SRR14092310.bracken -w breports/SRR14092310.breport
'
>>> Threshold: 10
>>> Number of species in sample: 1013
	>> Number of species with reads > threshold: 201
	>> Number of species with reads < threshold: 812
>>> Total reads in sample: 53998936
	>> Total reads kept at species level (reads > threshold): 23421667
	>> Total reads discarded (species reads < threshold): 2065
	>> Reads distributed: 12747044
	>> Reads not distributed (eg. no species above threshold): 193
	>> Unclassified reads: 17827967 
'
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