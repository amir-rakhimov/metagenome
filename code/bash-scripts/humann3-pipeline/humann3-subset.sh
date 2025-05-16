#!/usr/bin/env bash
### Get the current date in yyyy-mm-dd format with date -I
### Remove the hyphen with sed (use global search with /g to remove all hyphens)
### Store as a variable date_var
date_var=$(date -I|sed 's/-//g')
nthreads=10
### directory with decontaminated fastq files
bowtie2_decontam_fastq_dir=data/bowtie2_decontam_fastq
# mkdir output/humann3_pipeline
humann_out_dir=output/humann3_pipeline
# mkdir data/humann3_db
# mkdir data/humann3_db/humann3_full_${date_var}
# humann3_db_dir=data/humann3_db/humann3_full_${date_var}
humann3_db_dir=data/humann3_db/humann3_full_20240208
# Set conda channel priority
# conda config --add channels defaults
# conda config --add channels bioconda
# conda config --add channels conda-forge
# conda config --add channels biobakery
source ~/miniconda3/etc/profile.d/conda.sh

# Create a conda environment
# conda create -yqn humann-tools-3.7 -c biobakery python=3.7 humann=3.7 metaphlan
conda activate humann-tools-3.7 
## Testing installation: run HUMAnN unit tests
# humann_test
## Download the ChocoPhlAn database (16.4 GB)
# humann_databases --download chocophlan full ${humann3_db_dir}
# NOTE: The humann config file will be updated to point to this location for the default 
# chocophlan database. If you move this database, please use the "humann_config" command 
# to update the default location of this database. Alternatively you can always provide 
# the location of the chocophlan database you would like to use with the
#  "--nucleotide-database " option to humann.

## Download a translated search database
## To download the full UniRef90 database (20.7GB, recommended)
# humann_databases --download uniref uniref90_diamond ${humann3_db_dir}
# NOTE: The humann config file will be updated to point to this location for the default 
# uniref database. If you move this database, please use the humann_config command to
#  update the default location of this database. Alternatively you can always provide the 
#  location of the uniref database you would like to use with the 
#  --protein-database <uniref> option to humann.

# Note: Please make sure that all the databases in the directory are of same version.
#  Different versions of database must be separated in different directory.

# Install metaphlan database
# metaphlan --install

# # Take a subset of data
for FILE in ${bowtie2_decontam_fastq_dir}/*decontam_R1.fastq.gz
do 
 SAMPLE=$(echo ${FILE} | sed "s/_decontam_R1\.fastq\.gz//")
 base_name=$(basename "$SAMPLE" )
 zcat ${bowtie2_decontam_fastq_dir}/${base_name}_decontam_R1.fastq.gz |head -n 8 | \
 gzip -1 --fast > ${bowtie2_decontam_fastq_dir}/${base_name}_subset_R1.fastq.gz
 zcat ${bowtie2_decontam_fastq_dir}/${base_name}_decontam_R2.fastq.gz |head -n 8 | \
 gzip -1 --fast > ${bowtie2_decontam_fastq_dir}/${base_name}_subset_R2.fastq.gz;
done
# Merge subset data
for FILE in ${bowtie2_decontam_fastq_dir}/*subset_R1.fastq.gz
do 
 SAMPLE=$(echo ${FILE} | sed "s/_subset_R1\.fastq\.gz//")
 base_name=$(basename "$SAMPLE" )
 zcat ${bowtie2_decontam_fastq_dir}/${base_name}_subset_R1.fastq.gz \
 ${bowtie2_decontam_fastq_dir}/${base_name}_subset_R2.fastq.gz | \
 gzip -1 --fast > ${bowtie2_decontam_fastq_dir}/${base_name}_subset_merged.fastq.gz;
done

# Humann
for FILE in ${bowtie2_decontam_fastq_dir}/*merged.fastq.gz
do 
 SAMPLE=$(echo ${FILE} | sed "s/_merged\.fastq\.gz//")
 base_name=$(basename "$SAMPLE" )
 humann  --threads ${nthreads} \
	-i ${bowtie2_decontam_fastq_dir}/${base_name}_merged.fastq.gz \
	-o ${humann_out_dir}/humann_out_${date_var} \
	--nucleotide-database ${humann3_db_dir}/chocophlan \
    --protein-database ${humann3_db_dir}/uniref;
done

