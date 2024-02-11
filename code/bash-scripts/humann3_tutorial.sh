nthreads=8
source ~/miniconda3/etc/profile.d/conda.sh

# Set conda channel priority
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --add channels biobakery

# Create a conda environment
conda create --name humann-tools-3.7 -c biobakery python=3.7 humann=3.7 metaphlan=3.1

conda activate humann-tools-3.7 
# Testing installation: run HUMAnN unit tests
humann_test
# Download demo databases
humann_databases --download chocophlan DEMO humann_dbs
humann_databases --download uniref DEMO_diamond humann_dbs
# HUMAnN uses MetaPhlAn to detect organisms in the community. It should be installed
metaphlan --version

# 1. Metagenome functional profiling
mkdir data 
mkdir output
cd data
wget https://github.com/biobakery/humann/raw/master/examples/demo.fastq.gz
wget https://github.com/biobakery/humann/raw/master/examples/demo.m8
wget https://github.com/biobakery/humann/raw/master/examples/demo.sam
cd ..

# in humann_dbs
cd humann_dbs/chocophlan
gunzip -k *
cat * >refdb.ffn
bowtie2-build refdb.ffn refdb
cd ..
mkdir index
mv humann_dbs/chocophlan/refdb* index/

# Running HUMAnN
humann --input data/demo.fastq.gz \
  --output output/demo_fastq \
  --threads ${nthreads} 

# For custom metaphlan database
# --metaphlan-options "--index $index --bowtie2db your_path_to_db/MetaPhlan4"
# where index is an index from MetaPhlan database installed 
# (check the name of the file with “.tar” extention and copy it without extension).
# humann --input demo.fastq.gz --output output/demo_fastq \
#   --metaphlan-options "--index refdb --bowtie2db index/" # didn't work

#  column -t -s $'\t' file.tsv | less -S

# Need to concatenate paired end reads
cat sample_R1.fq sample_R2.fq > merge_sample.fq
# For multiple samples: go to metaphlan_protocol
for f in data/fastq/*.fasta.gz; do humann -i $f -o hmp_subset; done

humann_join_tables -i hmp_subset -o hmp_subset_genefamilies.tsv --file_name genefamilies
### Notice how we don't specify each sample but the name of directory 
### in front of _genefamilies.tsv
### file_name is also just genefamilies

# Normalizing RPKs to relative abundance
humann_renorm_table --input demo_fastq/demo_genefamilies.tsv \
    --output demo_fastq/demo_genefamilies-cpm.tsv --units cpm --update-snames
# Or for multiple samples
humann_renorm_table -i hmp_subset_genefamilies.tsv \
  -o hmp_subset_genefamilies-cpm.tsv \
  --units cpm \
  --update-snames

# Regrouping genes to other functional categories
humann_regroup_table --input demo_fastq/demo_genefamilies-cpm.tsv \
    --output demo_fastq/rxn-cpm.tsv --groups uniref90_rxn

# Attaching names to features
humann_rename_table --input demo_fastq/rxn-cpm.tsv \
    --output demo_fastq/rxn-cpm-named.tsv --names metacyc-rxn

