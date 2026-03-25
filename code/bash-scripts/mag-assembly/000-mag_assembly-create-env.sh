#!/usr/bin/env bash
#SBATCH -t 10:00:00
#SBATCH -J 20250608_16-30-mag_assembly-create-env
#SBATCH --output jobreports/20250608_16-30-mag_assembly-create-env-out.txt
#SBATCH --error jobreports/20250608_16-30-mag_assembly-create-env-out.txt
# This script performs the following:
# 1. Creates a conda environment for MAG assembly

# 2. Downloads CheckM2 database.

# Output directory: `~/projects/metagenome/data/checkm2_db/`

# 3. Dowloads GTDB-TK reference database

# Output directory: `~/miniconda3/envs/mag_assembly-tools/share/gtdbtk-1.7.0/db/`

# 4. Installs QUAST

source ~/miniconda3/etc/profile.d/conda.sh
# Install all tools with conda
## Setup channels for conda installation
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# ## Create a conda environment
# mamba create -yqn mag_assembly-tools -c conda-forge -c bioconda \
#      megahit prodigal checkm2 gtdbtk metaeuk \
#     pandas bbmap augustus sepp hmmer biopython r-base r-ggplot2 \
#     matplotlib prokka drep mmseqs2 coverm comparem anaconda::libtiff

# # Build MetaBAT locally
# # Install openssl first. Then, install metabat2
# cd ~
# wget https://bitbucket.org/berkeleylab/metabat/get/master.tar.gz
# mv master.tar.gz metabat2_master.tar.gz
# tar xzvf master.tar.gz
# cd berkeleylab-metabat-*
# mkdir build && cd build 
# cmake \
#   -DCMAKE_INSTALL_PREFIX=$HOME/berkeleylab-metabat-a51200652a89 \
#   -DBoost_NO_SYSTEM_PATHS=ON \
#   -DBoost_INCLUDE_DIR=$HOME/boost_1_88_0/include \
#   -DBoost_LIBRARY_DIR=$HOME/boost_1_88_0/lib \
#   -DBoost_PROGRAM_OPTIONS_LIBRARY_RELEASE=$HOME/boost_1_88_0/lib/libboost_program_options.so \
#   -DBoost_GRAPH_LIBRARY_RELEASE=$HOME/boost_1_88_0/lib/libboost_graph.so \
#   -DBoost_IOSTREAMS_LIBRARY_RELEASE=$HOME/boost_1_88_0/lib/libboost_iostreams.so \
#   -DBoost_DEBUG=ON \
#   ..

# make
# make test 
# make install

# conda activate mag_assembly-tools
# cd ~/projects/metagenome
# # Dowload GTDB-TK reference database
# export GTDBTK_DATA_PATH=~/miniconda3/envs/mag_assembly-tools/share/gtdbtk-2.4.0/db/release226
# download-db.sh

# # Install QUAST
# wget https://github.com/ablab/quast/releases/download/quast_5.2.0/quast-5.2.0.tar.gz -O ~/quast-5.2.0.tar.gz
# tar -xzf ~/quast-5.2.0.tar.gz -C ~/

# For CheckM2
conda activate mag_assembly-tools
wget https://ani.jgi.doe.gov/download_files/ANIcalculator_v1.tgz -O ~/ANIcalculator_v1.tgz
tar -xzf ~/ANIcalculator_v1.tgz -C ~/
export PATH="/home/rakhimov/ANIcalculator_v1/ANIcalculator:$PATH" 
conda deactivate 

# Get KofamScan
mkdir -p ~/kofamscan/db
cd ~/kofamscan/db
wget ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz 
wget ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz 
gunzip ko_list.gz 
tar xvzf profiles.tar.gz 
mkdir -p ~/kofamscan/bin
cd ~/kofamscan/bin
wget ftp://ftp.genome.jp/pub/tools/kofam_scan/kofam_scan-1.3.0.tar.gz
tar xvzf kofam_scan-1.3.0.tar.gz

# Get HMMer for KofamScan
cd ~/kofamscan 
mkdir hmmer src 
cd src 
wget http://eddylab.org/software/hmmer/hmmer.tar.gz 
tar xvzf hmmer.tar.gz 
cd hmmer-3.4
./configure --prefix=$HOME/kofamscan/hmmer 
make 
make install 

# Get ruby
cd ~/kofamscan
mkdir ruby
cd src
wget https://cache.ruby-lang.org/pub/ruby/3.4/ruby-3.4.4.tar.gz
cd ~/kofamscan/src 
tar xvzf ruby-3.4.4.tar.gz
cd ruby-3.4.4
./configure --prefix=$HOME/kofamscan/ruby 
make 
make install 
export PATH=$HOME/kofamscan/ruby/bin:$PATH 

# Create a config.yml
cd ~/kofamscan/bin/kofam_scan-1.3.0
cp config-template.yml config.yml
echo "profile: /home/rakhimov/kofamscan/db/profiles " >> config.yml
echo "ko_list: /home/rakhimov/kofamscan/db/ko_list " >> config.yml 
echo "hmmsearch: /home/rakhimov/kofamscan/hmmer/bin/hmmsearch  " >> config.yml 
echo "parallel: /usr/bin/parallel" >> config.yml 

conda install -n mag_assembly-tools  -c conda-forge -c bioconda -c defaults comparem

# Install coverm separately
mamba create -yqn coverm-env -c conda-forge -c bioconda \
    coverm samtools strobealign minimap2 bwa-mem2 skani fastani dashing



# Get the envrionment packages
conda env export mag_assembly-tools >  mag_assembly-tools-env.yaml

# Metagenomic thermometer
mamba create -yqn meta-thermo-env -c conda-forge -c bioconda \
    fastp seqkit prodigal python=3.6.5 numpy=1.19.0
