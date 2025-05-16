#!/usr/bin/env bash
#SBATCH -t 10:00:00
#SBATCH -J 20250514_15-52-mag_assembly-create-env
#SBATCH --output jobreports/20250514-15-52-mag_assembly-create-env-out.txt
#SBATCH --error jobreports/20250514-15-52-mag_assembly-create-env-out.txt
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
#     matplotlib prokka
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
conda install -n mag_assembly-tools  -c conda-forge -c bioconda -c defaults prokka
