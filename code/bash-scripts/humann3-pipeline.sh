#!/usr/bin/env bash
### Get the current date in yyyy-mm-dd format with date -I
### Remove the hyphen with sed (use global search with /g to remove all hyphens)
### Store as a variable date_var
date_var=$(date -I|sed 's/-//g')
time_var=$(date +%T |sed 's/:/_/g' )
date_time=${date_var}_${time_var}
nthreads=8
project_home_dir=~/projects/metagenome
### directory with decontaminated fastq files
bowtie2_decontam_fastq_dir=data/bowtie2_decontam_fastq
# mkdir output/humann3_pipeline
humann_out_dir=output/humann3_pipeline
# mkdir data/humann3_db
humann3_db_dir=data/humann3_db/humann3_full_20240513
source ~/miniconda3/etc/profile.d/conda.sh
# # Set conda channel priority
# conda config --add channels defaults
# conda config --add channels bioconda
# conda config --add channels conda-forge
# conda config --add channels biobakery

# Create a conda environment
# conda create -yqn humann-tools-3.7 -c biobakery python=3.7 humann=3.7 metaphlan
conda activate humann-tools-3.7 
# Testing installation: run HUMAnN unit tests
# humann_test
## Download the ChocoPhlAn database (16.4 GB)
# humann_databases --download chocophlan full ${humann3_db_dir}
# NOTE: The humann config file will be updated to point to this location for the default 
# chocophlan database. If you move this database, please use the "humann_config" command 
# to update the default location of this database. Alternatively you can always provide 
# the location of the chocophlan database you would like to use with the
#  "--nucleotide-database " option to humann.
# humann_config --update database_folders protein $DIR 
# humann_config --update database_folders protein \
	# /lustre7/home/rakhimov/projects/metagenome/data/humann3_db/humann3_full_20240513/uniref
# humann_config --update database_folders nucleotide \
	# /lustre7/home/rakhimov/projects/metagenome/data/humann3_db/humann3_full_20240513/chocophlan

## Download a translated search database
## To download the full UniRef90 database (20.7GB, recommended)
# humann_databases --download uniref uniref90_diamond ${humann3_db_dir}
# NOTE: The humann config file will be updated to point to this location for the default 
# uniref database. If you move this database, please use the humann_config command to
#  update the default location of this database. Alternatively you can always provide the 
#  location of the uniref database you would like to use with the 
#  --protein-database <uniref> option to humann.

# Note: Please make sure that all the databases in the directory are of same version.
#  Different versions of database must be separated in different directory

# Install metaphlan database
# metaphlan --install

# bowtie2-build --large-index --threads 12 \
# 	~/miniconda3/envs/humann-tools-3.7/lib/python3.7/site-packages/metaphlan/metaphlan_databases/mpa_vOct22_CHOCOPhlAnSGB_202212_SGB.fna \
# 	~/miniconda3/envs/humann-tools-3.7/lib/python3.7/site-packages/metaphlan/metaphlan_databases/mpa_vOct22_CHOCOPhlAnSGB_202212_SGB

# Download the ChocoPhlAn database (16.4 GB)
# humann_databases --download chocophlan full ${humann3_db_dir}
# Take a subset of data
# for FILE in ${bowtie2_decontam_fastq_dir}/*decontam_R1.fastq.gz
# do 
#  SAMPLE=$(echo ${FILE} | sed "s/_decontam_R1\.fastq\.gz//")
#  base_name=$(basename "$SAMPLE" )
#  zcat ${bowtie2_decontam_fastq_dir}/${base_name}_decontam_R1.fastq.gz |head -n 40| \
#  gzip -9 --best > ${bowtie2_decontam_fastq_dir}/${base_name}_subset_R1.fastq.gz
#  zcat ${bowtie2_decontam_fastq_dir}/${base_name}_decontam_R2.fastq.gz |head -n 40| \
#  gzip -9 --best > ${bowtie2_decontam_fastq_dir}/${base_name}_subset_R2.fastq.gz;
# done
# # Merge subset data
# for FILE in ${bowtie2_decontam_fastq_dir}/*subset_R1.fastq.gz
# do 
#  SAMPLE=$(echo ${FILE} | sed "s/_subset_R1\.fastq\.gz//")
#  base_name=$(basename "$SAMPLE" )
#  zcat ${bowtie2_decontam_fastq_dir}/${base_name}_subset_R1.fastq.gz \
#  ${bowtie2_decontam_fastq_dir}/${base_name}_subset_R2.fastq.gz | \
#  gzip -9 --best > ${bowtie2_decontam_fastq_dir}/${base_name}_subset_merged.fastq.gz;
# done

# Full data
# Need to merge paired end reads by concatenation
# for FILE in ${bowtie2_decontam_fastq_dir}/*decontam_R1.fastq.gz
# do 
#  SAMPLE=$(echo ${FILE} | sed "s/_decontam_R1\.fastq\.gz//")
#  base_name=$(basename "$SAMPLE" )
#  zcat ${bowtie2_decontam_fastq_dir}/${base_name}_decontam_R1.fastq.gz \
#   ${bowtie2_decontam_fastq_dir}/${base_name}_decontam_R2.fastq.gz | \
#   gzip -9 --best > ${bowtie2_decontam_fastq_dir}/${base_name}_merged.fastq.gz;
# done

# Humann
### We skip prescreen with Metaphlan because our community is not well-characterised. 
### We avoid loss of data.
### Humann will create a custom ChocoPhlan database, run bowtie2-build, bowtie2,
### then perform nucleotide alignment, and align reads to uniref90 with diamond
# for FILE in ${bowtie2_decontam_fastq_dir}/*merged.fastq.gz
# do 
#  SAMPLE=$(echo ${FILE} | sed "s/_merged\.fastq\.gz//")
#  base_name=$(basename "$SAMPLE" )
#  humann  --threads ${nthreads} \
# 	-i ${bowtie2_decontam_fastq_dir}/${base_name}_merged.fastq.gz \
# 	-o ${humann_out_dir}/humann_out_${date_time}\
# 	--nucleotide-database ${humann3_db_dir}/chocophlan \
#     --protein-database ${humann3_db_dir}/uniref \
# 	--bypass-prescreen
# #  rm ${bowtie2_decontam_fastq_dir}/${base_name}_merged.fastq.gz
# #  cp ${humann_out_dir}/humann_out_${date_time}/${base_name}_merged_humann_temp/${base_name}_merged_metaphlan_bugs_list.tsv  ${humann_out_dir}/humann_out_${date_time}/
#  cp ${humann_out_dir}/humann_out_${date_time}/${base_name}_merged_humann_temp/${base_name}_merged.log  ${humann_out_dir}/humann_out_${date_time}/
#  rm -rf ${humann_out_dir}/humann_out_${date_time}/${base_name}_merged_humann_temp;
# done

# humann -i data/bowtie2_decontam_fastq/2D10_subset_merged.fastq.gz \
# 	-o output/humann3_pipeline/humann_out_${date_time}\
# 	--protein-database data/humann3_db/humann3_full_20240513/uniref \
# 	--bypass-prescreen

# 	--nucleotide-database data/humann3_db/humann3_full_20240513/chocophlan \
# We can skip creation of custom Chocophlan database if you already have it using
#  --nucleotide-database $OUTPUT_DIR/$SAMPLE_1_humann_temp/ --bypass-nucleotide-index
# (I run out of memory during bowtie2 build-index and I don't want to create the database 
# again
for FILE in ${bowtie2_decontam_fastq_dir}/*merged.fastq.gz
do 
 SAMPLE=$(echo ${FILE} | sed "s/_merged\.fastq\.gz//")
 base_name=$(basename "$SAMPLE" )
 humann  --threads ${nthreads} \
	-i ${bowtie2_decontam_fastq_dir}/${base_name}_merged.fastq.gz \
	-o ${humann_out_dir}/humann_out_${date_time}\
	--nucleotide-database ${humann3_db_dir}/chocophlan \
	--protein-database ${humann3_db_dir}/uniref \
	--bypass-prescreen 
 cp ${humann_out_dir}/humann_out_${date_time}/${base_name}_merged_humann_temp/${base_name}_merged.log  ${humann_out_dir}/humann_out_${date_time}
 cd ${humann_out_dir}/humann_out_${date_time}/
 rm ${base_name}_merged_humann_temp/*bt2l
 rm ${base_name}_merged_humann_temp/${base_name}_merged_custom_chocophlan_database.ffn
 cd ${project_home_dir}/${humann_out_dir}/humann_out_${date_time}
 tar -zcvf ${base_name}_merged_humann_temp.tar ${base_name}_merged_humann_temp/.
 gzip -9 --best ${base_name}_merged_humann_temp.tar
 cd ${project_home_dir};
done

# SAMPLE is data/bowtie2_decontam_fastq/2D14_trim
# base_name is 2D14_trim

# for FILE in ${bowtie2_decontam_fastq_dir}/*merged.fastq.gz
# do 
#  SAMPLE=$(echo ${FILE} | sed "s/_merged\.fastq\.gz//")
#  base_name=$(basename "$SAMPLE" )
#  humann  --threads ${nthreads} \
# 	-i ${bowtie2_decontam_fastq_dir}/${base_name}_merged.fastq.gz \
# 	-o ${humann_out_dir}/humann_out_${date_time}\
#     --protein-database ${humann3_db_dir}/uniref \
# 	--bypass-prescreen \
# 	--nucleotide-database ${humann_out_dir}/humann_out_20240516_17_20_40/2D10_trim_merged_humann_temp/ \
# 	--bypass-nucleotide-index
#  cp ${humann_out_dir}/humann_out_${date_time}/${base_name}_merged_humann_temp/${base_name}_merged.log  ${humann_out_dir}/humann_out_${date_time};
# done
# humann  --threads ${nthreads} \
# 	-i ${bowtie2_decontam_fastq_dir}/2D10_trim_subset_merged.fastq.gz \
# 	-o ${humann_out_dir}/humann_out_trial

### Each fastq input file will have its own output directory.
### The intermediate files for that input file will be in the sub-directory
### But the three main analysis types (gene families, pathway abundances, 
### and pathway coverages) will be in the general output directory,
### i.e. outside sample sub-directory

# SRS014459-Stool_humann_temp # temp files are in a sub-directory
# SRS014472-Buccal_mucosa_humann_temp # temp files are in a sub-directory

# SRS014459-Stool_genefamilies.tsv # main output is in general output directory
# SRS014472-Buccal_mucosa_genefamilies.tsv # main output is in general output directory
# SRS014459-Stool_pathabundance.tsv # main output is in general output directory
# SRS014472-Buccal_mucosa_pathabundance.tsv # main output is in general output directory
# SRS014459-Stool_pathcoverage.tsv # main output is in general output directory
# SRS014472-Buccal_mucosa_pathcoverage.tsv # main output is in general output directory
#mkdir ${humann_out_dir}/humann_out_${date_time}_merged
#humann_merged_out_dir=${humann_out_dir}/humann_out_${date_time}_merged
### Let's join gene family abundance outputs into a single table
#humann_join_tables -i ${humann_out_dir}/humann_out_${date_time} \
#  -o ${humann_merged_out_dir}/${date_time}_humann_out_genefamilies.tsv \
#  --file_name genefamilies
### Notice how we don't specify each sample but the name of directory 
### in front of _genefamilies.tsv
### file_name is also just genefamilies: only join tables with this string included in the file name
### Result: ${humann_out_dir}/humann_out_${date_time}_merged/${date_time}_humann_out_genefamilies.tsv
#humann_join_tables -i ${humann_out_dir}/humann_out_${date_time} \
#  -o ${humann_merged_out_dir}/${date_time}_humann_out_pathabundance.tsv \
#  --file_name pathabundance

# Normalizing RPKs to relative abundance
#humann_renorm_table --input ${humann_merged_out_dir}/${date_time}_humann_out_genefamilies.tsv \
#    --output ${humann_merged_out_dir}/${date_time}_humann_out_genefamilies-cpm.tsv \
#	--units cpm \
#	--update-snames
# for pathabundance renorm is enough
#humann_renorm_table --input ${humann_merged_out_dir}/${date_time}_humann_out_pathabundance.tsv \
#    --output ${humann_merged_out_dir}/${date_time}_humann_out_pathabundance-cpm.tsv \
#	--units cpm \
#	--update-snames
# Regrouping genes to other functional categories
#humann_regroup_table --input ${humann_merged_out_dir}/${date_time}_humann_out_genefamilies-cpm.tsv \
#    --output ${humann_merged_out_dir}/${date_time}_humann_out_genefamilies-cpm-rxn.tsv \
#	--groups uniref90_rxn

# Attaching names to features
#humann_rename_table --input ${humann_merged_out_dir}/${date_time}_humann_out_genefamilies-cpm-rxn.tsv \
#    --output ${humann_merged_out_dir}/${date_time}_humann_out_genefamilies-cpm-rxn-named.tsv \
#	--names metacyc-rxn


# Remove unmapped and ungrouped data from gene families
grep -v -e UNMAPPED -e UNGROUPED ${humann_merged_out_dir}/${date_time}_humann_out_genefamilies-cpm-rxn-named.tsv >${humann_merged_out_dir}/${date_time}_humann_out_genefamilies-cpm-rxn-named-filtered.tsv
# Extract stratifications
grep "|" ${humann_merged_out_dir}/${date_time}_humann_out_genefamilies-cpm-rxn-named-filtered.tsv >${humann_merged_out_dir}/${date_time}_humann_out_genefamilies-cpm-rxn-named-filtered-strat.tsv
# Remove stratifications
grep -v "|" ${humann_merged_out_dir}/${date_time}_humann_out_genefamilies-cpm-rxn-named-filtered.tsv >${humann_merged_out_dir}/${date_time}_humann_out_genefamilies-cpm-rxn-named-filtered-total.tsv

# Remove unmapped and ungrouped data from pathabundance
grep -v -e UNMAPPED -e UNINTEGRATED ${humann_merged_out_dir}/${date_time}_humann_out_pathabundance-cpm.tsv >${humann_merged_out_dir}/${date_time}_humann_out_pathabundance-cpm-filtered.tsv
# Extract stratifications
grep "|" ${humann_merged_out_dir}/${date_time}_humann_out_pathabundance-cpm-filtered.tsv >${humann_merged_out_dir}/${date_time}_humann_out_pathabundance-cpm-filtered-strat.tsv
# Remove stratifications
grep -v "|" ${humann_merged_out_dir}/${date_time}_humann_out_pathabundance-cpm-filtered.tsv >${humann_merged_out_dir}/${date_time}_humann_out_pathabundance-cpm-filtered-total.tsv

# Do it in reverse: extract stratifications and totals, then filter (gene families)
grep "|" ${humann_merged_out_dir}/${date_time}_humann_out_genefamilies-cpm-rxn-named.tsv >${humann_merged_out_dir}/${date_time}_humann_out_genefamilies-cpm-rxn-named-strat.tsv
grep -v "|" ${humann_merged_out_dir}/${date_time}_humann_out_genefamilies-cpm-rxn-named.tsv >${humann_merged_out_dir}/${date_time}_humann_out_genefamilies-cpm-rxn-named-total.tsv
grep -v -e UNMAPPED -e UNGROUPED ${humann_merged_out_dir}/${date_time}_humann_out_genefamilies-cpm-rxn-named-total.tsv >${humann_merged_out_dir}/${date_time}_humann_out_genefamilies-cpm-rxn-named-total-filtered.tsv

# Do it in reverse: extract stratifications and totals, then filter (pathabundance)
grep "|" ${humann_merged_out_dir}/${date_time}_humann_out_pathabundance-cpm.tsv > ${humann_merged_out_dir}/${date_time}_humann_out_pathabundance-cpm-strat.tsv
grep -v "|" ${humann_merged_out_dir}/${date_time}_humann_out_pathabundance-cpm.tsv  > ${humann_merged_out_dir}/${date_time}_humann_out_pathabundance-cpm-total.tsv
grep -v -e UNMAPPED -e UNINTEGRATED ${humann_merged_out_dir}/${date_time}_humann_out_pathabundance-cpm-total.tsv > ${humann_merged_out_dir}/${date_time}_humann_out_pathabundance-cpm-total-filtered.tsv

# Renormalizing after filtering
humann_renorm_table --input ${humann_merged_out_dir}/${date_time}_humann_out_genefamilies-cpm-rxn-named-filtered-total.tsv \
   --output ${humann_merged_out_dir}/${date_time}_humann_out_genefamilies-cpm-rxn-named-filtered-total-renorm.tsv \
	--units cpm 
humann_renorm_table --input ${humann_merged_out_dir}/${date_time}_humann_out_pathabundance-cpm-filtered-total.tsv \
   --output ${humann_merged_out_dir}/${date_time}_humann_out_pathabundance-cpm-filtered-total-renorm.tsv \
	--units cpm 


humann_renorm_table --input ${humann_merged_out_dir}/${date_time}_humann_out_pathabundance-cpm-total-filtered.tsv \
   --output ${humann_merged_out_dir}/${date_time}_humann_out_pathabundance-cpm-total-filtered-renorm.tsv \
	--units cpm 

