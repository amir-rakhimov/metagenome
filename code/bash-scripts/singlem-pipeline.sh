#!/usr/bin/env bash
### Get the current date in yyyy-mm-dd format with date -I
### Remove the hyphen with sed (use global search with /g to remove all hyphens)
### Store as a variable date_var
singlem_db_date=20240926
# time format is 16_28_52 
time_var=$(date +%T |sed 's/:/_/g' )
date_var=$(date -I|sed 's/-//g')
date_time=${date_var}_${time_var}
# date_time=20240926_17_32_40
nthreads=2
# mem_req=8G
project_home_dir=~/projects/metagenome
fastq_dir=data/fastq/yasuda-fastq
singlem_db_dir=data/singlem_db/singlem_data_${singlem_db_date}
singlem_output_dir=output/singlem_pipeline/singlem_output
singlem_profiles_dir=output/singlem_pipeline/singlem_output/singlem_profiles
singlem_otu_tables_dir=output/singlem_pipeline/singlem_output/singlem_otu_tables
#
# mkdir data/singlem_db
# mkdir data/singlem_db/singlem_data_${singlem_db_date}
# mkdir output/singlem_pipeline
# mkdir output/singlem_pipeline/singlem_output
# mkdir output/singlem_pipeline/singlem_output/singlem_profiles
# mkdir output/singlem_pipeline/singlem_output/singlem_otu_tables


# Installation of SingleM
source ~/miniconda3/etc/profile.d/conda.sh
# conda config --set channel_priority flexible
# mamba create -c bioconda -c conda-forge --name singlem singlem'>='0.18.3
# Test installation
conda activate singlem
# singlem -h
# Download data with 'singlem data' command (1.4 GB)
# 09/26/2024 03:30:22 PM INFO: SingleM v0.18.3
# 09/26/2024 03:30:22 PM INFO: Downloading data version 4.3.0 with ZenodoBackpack from 10.5281/zenodo.5739611 ..
# 09/26/2024 03:30:25 PM INFO: Downloading https://zenodo.org/api/records/11323477/files/S4.3.0.GTDB_r220.metapackage_20240523.smpkg.zb.tar.gz/content to data/singlem_db/singlem_data_20240926/S4.3.0.GTDB_r220.metapackage_20240523.smpkg.zb.tar.gz.
# 100%|█████████████████████████████████████████████████████████████████████████████████| 1.40G/1.40G [08:02<00:00, 2.89MiB/s]
# 09/26/2024 03:38:31 PM INFO: Extracting files from archive...
# 09/26/2024 03:39:13 PM INFO: Verifying version and checksums...
# 09/26/2024 03:39:24 PM INFO: Verification success.
# 09/26/2024 03:39:24 PM INFO: Finished downloading data
# 09/26/2024 03:39:24 PM INFO: The environment variable SINGLEM_METAPACKAGE_PATH can now be set to data/singlem_db/singlem_data_20240926
# 09/26/2024 03:39:24 PM INFO: For instance, the following can be included in your .bashrc (requires logout and login after inclusion):
# 09/26/2024 03:39:24 PM INFO: export SINGLEM_METAPACKAGE_PATH='/home/rakhimov/projects/metagenome/data/singlem_db/singlem_data_20240926/S4.3.0.GTDB_r220.metapackage_20240523.smpkg.zb'
# singlem data --output-directory ${singlem_db_dir}
# You must set the environment variable to use the data! Add this line to .bashrc file
# export SINGLEM_METAPACKAGE_PATH='~/projects/metagenome/data/singlem_db/singlem_data_20240926/S3.0.5.metapackage20220806.smpkg.zb'
# for supercomputer :
# export SINGLEM_METAPACKAGE_PATH='/lustre7/home/rakhimov/projects/metagenome/data/singlem_db/singlem_data_20240926/S4.3.0.GTDB_r220.metapackage_20240523.smpkg.zb'

# To test the main subcommand of SingleM, pipe works, download a minimal dataset (2.5 MB)
# and generate a taxonomic profile like so:
# wget -P ${singlem_output_dir} 'https://github.com/wwood/singlem/raw/44e1f81404c12931742259088999290edbb271b3/test/data/methanobacteria/genomes/GCA_000309865.1_genomic.fna'
# singlem pipe -1 ${singlem_output_dir}/GCA_000309865.1_genomic.fna -p ${singlem_output_dir}/testdata_output.tsv
# singlem pipe -1 ${fastq_dir}/2D10_wms_L3_1.fq.gz -2 ${fastq_dir}/2D10_wms_L3_2.fq.gz \
#     --otu-table ${singlem_output_dir}/2D10_wms_otu_table.tsv \
#     --output-extras \
#     -p ${singlem_output_dir}/2D10_wms_profile.tsv \
#     --threads ${nthreads}

# singlem pipe -1  data/2D10_raw_subset_R1.fq.gz  \
#     -2 data/2D10_raw_subset_R2.fq.gz  \
#     --otu-table ${singlem_output_dir}/2D10_subset_otu_table.tsv \
#     --output-extras \
#     -p ${singlem_output_dir}/2D10_subset_profile.tsv \
#     --threads ${nthreads}

for FILE in ${fastq_dir}/*wms_L3_1.fq.gz
do 
 SAMPLE=$(echo ${FILE} | sed "s/_wms_L3_1\.fq\.gz//")
 base_name=$(basename "$SAMPLE" )
 singlem pipe -1 ${fastq_dir}/${base_name}_wms_L3_1.fq.gz -2 ${fastq_dir}/${base_name}_wms_L3_2.fq.gz \
  -p ${singlem_profiles_dir}/${date_time}_${base_name}_wms_profile.tsv \
  --otu-table ${singlem_otu_tables_dir}/${date_time}_${base_name}_otu_table.tsv \
  --output-extras \
  --threads ${nthreads};
done   2>&1 |tee ${date_time}_singlem_pipeline_extras_stdout.txt   


# -1, --forward, --reads, --sequences sequence_file [sequence_file ...]
#     nucleotide read sequence(s) (forward or unpaired) to be searched. 
#     Can be FASTA or FASTQ format, GZIP-compressed or not.

# -2, --reverse sequence_file [sequence_file ...]
#     reverse reads to be searched. Can be FASTA or FASTQ format, GZIP-compressed or not.

# -p, --taxonomic-profile FILE
#     output a 'condensed' taxonomic profile for each sample based on the OTU table. 
#     Taxonomic profiles output can be further converted to other formats using singlem summarise.

# --threads num_threads
#     number of CPUS to use [default: 1]

# Convert the output into relative abundances
for FILE in ${singlem_profiles_dir}/*wms_profile.tsv
do 
 SAMPLE=$(echo ${FILE} | sed "s/_wms_profile\.tsv//" |sed "s/${date_time}_//")
 base_name=$(basename "$SAMPLE" )
 singlem summarise --input-taxonomic-profile ${singlem_profiles_dir}/${date_time}_${base_name}_wms_profile.tsv \
    --output-taxonomic-profile-with-extras ${singlem_profiles_dir}/${date_time}_${base_name}_wms_profile_with_extras.tsv;
done   2>&1 |tee ${date_time}_singlem_sumarise_extras_stdout.txt   


# singlem summarise --input-taxonomic-profile ${singlem_output_dir}/2D10_wms_profile.tsv \
#     --output-species-by-site-relative-abundance-prefix ${singlem_output_dir}/2D10_wms_profile

# singlem summarise --input-taxonomic-profile ${singlem_output_dir}/2D10_wms_profile.tsv \
#     --output-taxonomic-profile-with-extras ${singlem_output_dir}/2D10_wms_profile.with_extras.tsv

# singlem summarise --input-taxonomic-profile  ${singlem_output_dir}/2D10_wms_profile.tsv \
#     --output-taxonomic-profile-krona  ${singlem_output_dir}/2D10_wms_profile.html

# Split the taxonomy information into separate columns
for FILE in ${singlem_profiles_dir}/*wms_profile_with_extras.tsv
do 
 SAMPLE=$(echo ${FILE} | sed "s/_wms_profile_with_extras\.tsv//" |sed "s/${date_time}_//")
 base_name=$(basename "$SAMPLE" )
 python3 code/python-scripts/parse_singlem_report.py \
  --input_file ${singlem_profiles_dir}/${date_time}_${base_name}_wms_profile_with_extras.tsv \
  --output_file ${singlem_profiles_dir}/${date_time}_${base_name}_wms_profile_clean.tsv;
done  

