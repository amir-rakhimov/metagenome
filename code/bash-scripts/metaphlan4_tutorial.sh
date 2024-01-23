# Set conda channel priority
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels biobakery

# Create a conda environment
conda create --name humann-tools-3.7 -c biobakery python=3.7 humann metaphlan
# We use the environment for HUMAnN because it uses metaphlan
conda activate humann-tools-3.7 
# Create the environment
conda create --name metaphlan-tools-4.0.2 -c bioconda python=3.7 metaphlan=4.0.2 biopython

conda activate metaphlan-tools-4.0.2

# 1. Download data
# Create a directory for the protocol
mkdir metaphlan_protocol/
cd metaphlan_protocol/
conda activate metaphlan-tools
mkdir data
mkdir output
mkdir data/fastq
cd fastq/
wget https://github.com/biobakery/MetaPhlAn/releases/download/4.0.2/SRS014476-Supragingival_plaque.fasta.gz
wget https://github.com/biobakery/MetaPhlAn/releases/download/4.0.2/SRS014494-Posterior_fornix.fasta.gz
wget https://github.com/biobakery/MetaPhlAn/releases/download/4.0.2/SRS014459-Stool.fasta.gz
wget https://github.com/biobakery/MetaPhlAn/releases/download/4.0.2/SRS014464-Anterior_nares.fasta.gz
wget https://github.com/biobakery/MetaPhlAn/releases/download/4.0.2/SRS014470-Tongue_dorsum.fasta.gz
wget https://github.com/biobakery/MetaPhlAn/releases/download/4.0.2/SRS014472-Buccal_mucosa.fasta.gz

metaphlan --install --index mpa_vJan21_CHOCOPhlAnSGB_202103
'Decompressing /home/amir/miniconda3/envs/metaphlan-tools-4.0.2/lib/python3.7/site-packages/
metaphlan/metaphlan_databases/mpa_vJan21_CHOCOPhlAnSGB_202103_SGB.fna.bz2 into 
/home/amir/miniconda3/envs/metaphlan-tools/lib/python3.7/site-packages/metaphlan/metaphlan_databases/
mpa_vJan21_CHOCOPhlAnSGB_202103_SGB.fna 
'
bowtie2-build ~/miniconda3/envs/metaphlan-tools-4.0.2/lib/python3.7/site-packages/metaphlan/metaphlan_databases/mpa_vJan21_CHOCOPhlAnSGB_202103.fna \
~/miniconda3/envs/metaphlan-tools-4.0.2/lib/python3.7/site-packages/metaphlan/metaphlan_databases/mpa_vJan21_CHOCOPhlAnSGB_202103

# 2. Running a single FASTQ sample
metaphlan data/fastq/SRS014476-Supragingival_plaque.fasta.gz \
	--input_type fasta \
	--bowtie2out output/bowtie2_output/SRS014476-Supragingival_plaque.fasta.gz.bowtie2out.txt \
	-o output/profiles/SRS014476-Supragingival_plaque_profile.txt \
	--nproc 4
metaphlan data/fastq/SRS014494-Posterior_fornix.fasta.gz --input_type fasta \
	-o output/profiles/SRS014494-Posterior_fornix_profile.txt \
	--bowtie2out output/bowtie2_output/SRS014494-Posterior_fornix.fasta.gz.bowtie2out.txt 

# Try the exact database and metaphlan version
metaphlan data/fastq/SRS014476-Supragingival_plaque.fasta.gz --input_type fasta \
 --bowtie2out output/bowtie2_output/SRS014476-Supragingival_plaque_202103.fasta.gz.bowtie2out.txt  \
 -o output/profiles/SRS014476-Supragingival_plaque_profile_202103.txt \
 --nproc 4 \
 -x mpa_vJan21_CHOCOPhlAnSGB_202103 \
  --bowtie2db ~/miniconda3/envs/metaphlan-tools-4.0.2/lib/python3.7/site-packages/metaphlan/metaphlan_databases/
# Re-profiling a sample using its bowtie2out file
# We change the --input_type since we are not starting from raw reads
metaphlan output/bowtie2_output/SRS014476-Supragingival_plaque.fasta.gz.bowtie2out.txt \
  --input_type bowtie2out \
  -o output/profiles/SRS014476-Supragingival_plaque_profile.txt
'WARNING: The metagenome profile contains clades that represent multiple species merged into a single representant.
An additional column listing the merged species is added to the MetaPhlAn output'

# Analysing multiple samples
# Use the basename command to remove directory from file name
directory="data/fastq/"
for file in "$directory"*.fasta.gz; do
    # Extract the base name without extension
    base_name=$(basename "$file" .fasta.gz);
       
    echo "$base_name"
	metaphlan ${directory}${base_name}.fasta.gz --input_type fasta --nproc 8 --bowtie2out output/bowtie2_output/${base_name}.fasta.gz.bowtie2out.txt -o output/profiles/${base_name}_profile.txt
done

# 3. Merging MetaPhlAn profiles
# Create a single tab-delimited table from a set of sample-specific abundance profiles 
# (isolating the sample names, feature taxonomies, and relative abundances)
merge_metaphlan_tables.py output/profiles/*_profile.txt > output/merged_abundance_table.txt

# 4. Visualisation
conda install -c biobakery hclust2

# Generate a species-only abundance table