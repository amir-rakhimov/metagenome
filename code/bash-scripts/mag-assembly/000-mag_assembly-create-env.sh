#!/usr/bin/env bash
#SBATCH -t 10:00:00
#SBATCH -J 20250608_16-30-mag_assembly-create-env
#SBATCH --output jobreports/20250608_16-30-mag_assembly-create-env-out.txt
#SBATCH --error jobreports/20250608_16-30-mag_assembly-create-env-out.txt
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
#     matplotlib prokka drep mmseqs2 coverm comparem

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

conda activate mag_assembly-tools
# conda install -n mag_assembly-tools  -c conda-forge -c bioconda -c defaults drep
# conda env export environment_name > environment.yaml
# conda install -n mag_assembly-tools anaconda::libtiff
wget https://ani.jgi.doe.gov/download_files/ANIcalculator_v1.tgz -O ~/ANIcalculator_v1.tgz
tar -xzf ~/ANIcalculator_v1.tgz -C ~/
export PATH="/home/rakhimov/ANIcalculator_v1/ANIcalculator:$PATH" 

# Get KofamScan
git clone https://github.com/takaram/kofam_scan ~/kofam_scan
mkdir -p data/kofam_scan_data
# wget ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz -O data/kofam_scan_data/profiles.tar.gz
#  wget ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz -O data/kofam_scan_data/ko_list.gz

wget https://www.genome.jp/ftp/db/kofam/profiles.tar.gz -O data/kofam_scan_data/profiles.tar.gz
wget https://www.genome.jp/ftp/db/kofam/ko_list.gz -O data/kofam_scan_data/ko_list.gz

tar -xzf data/kofam_scan_data/profiles.tar.gz -C data/kofam_scan_data
gunzip ko_list.gz 

cp ~/kofam_scan/config-template.yml ~/kofam_scan/config.yml 
echo "profile: /lustre10/home/rakhimov/projects/metagenome/data/kofam_scan_data/profiles " >> ~/kofam_scan/config.yml 
echo "ko_list: /lustre10/home/rakhimov/projects/metagenome/data/kofam_scan_data/ko_list " >> ~/kofam_scan/config.yml 
echo "hmmsearch: /lustre10/home/rakhimov/miniconda3/envs/mag_assembly-tools/bin/hmmsearch " >> ~/kofam_scan/config.yml 
echo "parallel: /lustre10/home/rakhimov/miniconda3/envs/mag_assembly-tools/bin/parallel " >> ~/kofam_scan/config.yml 

conda install -n mag_assembly-tools  -c conda-forge -c bioconda -c defaults comparem
