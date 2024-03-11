#!/usr/bin/env bash
### Get the current date in yyyy-mm-dd format with date -I
### Remove the hyphen with sed (use global search with /g to remove all hyphens)
### Store as a variable date_var
date_var=$(date -I|sed 's/-//g')
date_var=20240217
nthreads=10
### directory with decontaminated fastq files
bowtie2_decontam_fastq_dir=data/bowtie2_decontam_fastq
# mkdir output/humann3_pipeline
humann_out_dir=output/humann3_pipeline
# mkdir data/humann3_db
humann3_db_dir=data/humann3_db/humann3_full_20240208
humann_merged_out_dir=${humann_out_dir}/humann_merged_${date_var}
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

## Download a translated search database
## To download the full UniRef90 database (20.7GB, recommended)
# humann_databases --download uniref uniref90_diamond ${humann3_db_dir}
# NOTE: The humann config file will be updated to point to this location for the default 
# uniref database. If you move this database, please use the humann_config command to
#  update the default location of this database. Alternatively you can always provide the 
#  location of the uniref database you would like to use with the 
#  --protein-database <uniref> option to humann.

# wget https://huttenhower.sph.harvard.edu/humann_data/uniprot/uniref_annotated/uniref90_annotated_v201901b_full.tar.gz \
#  -O uniref90_annotated_v201901b_full.tar.gz 
# tar -xzf uniref90_annotated_v201901b_full.tar.gz --directory uniref/ 
# Note: Please make sure that all the databases in the directory are of same version.
#  Different versions of database must be separated in different directory

# Install metaphlan database
# metaphlan --install
# scp mpa_vOct22_CHOCOPhlAnSGB_202212* \
#  rakhimov@gw.ddbj.nig.ac.jp:~/miniconda3/envs/humann-tools-3.7/lib/python3.7/site-packages/metaphlan/metaphlan_databases
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
#  gzip -1 --fast > ${bowtie2_decontam_fastq_dir}/${base_name}_subset_R1.fastq.gz
#  zcat ${bowtie2_decontam_fastq_dir}/${base_name}_decontam_R2.fastq.gz |head -n 40| \
#  gzip -1 --fast > ${bowtie2_decontam_fastq_dir}/${base_name}_subset_R2.fastq.gz;
# done
# # Merge subset data
# for FILE in ${bowtie2_decontam_fastq_dir}/*subset_R1.fastq.gz
# do 
#  SAMPLE=$(echo ${FILE} | sed "s/_subset_R1\.fastq\.gz//")
#  base_name=$(basename "$SAMPLE" )
#  zcat ${bowtie2_decontam_fastq_dir}/${base_name}_subset_R1.fastq.gz \
#  ${bowtie2_decontam_fastq_dir}/${base_name}_subset_R2.fastq.gz | \
#  gzip -1 --fast > ${bowtie2_decontam_fastq_dir}/${base_name}_subset_merged.fastq.gz;
# done

# Full data
# Need to merge paired end reads by concatenation
for FILE in ${bowtie2_decontam_fastq_dir}/*decontam_R1.fastq.gz
do 
 SAMPLE=$(echo ${FILE} | sed "s/_decontam_R1\.fastq\.gz//")
 base_name=$(basename "$SAMPLE" )
 zcat ${bowtie2_decontam_fastq_dir}/${base_name}_decontam_R1.fastq.gz \
 ${bowtie2_decontam_fastq_dir}/${base_name}_decontam_R2.fastq.gz | \
 gzip -1 --fast > ${bowtie2_decontam_fastq_dir}/${base_name}_merged.fastq.gz;
done

# Humann
for FILE in ${bowtie2_decontam_fastq_dir}/*merged.fastq.gz
do 
 SAMPLE=$(echo ${FILE} | sed "s/_merged\.fastq\.gz//")
 base_name=$(basename "$SAMPLE" )
 humann  --threads ${nthreads} \
	-i ${bowtie2_decontam_fastq_dir}/${base_name}_merged.fastq.gz \
	-o ${humann_out_dir}/humann_out_${date_var}\
	--nucleotide-database ${humann3_db_dir}/chocophlan \
    --protein-database ${humann3_db_dir}/uniref
 rm ${bowtie2_decontam_fastq_dir}/${base_name}_merged.fastq.gz
 cp ${humann_out_dir}/humann_out_${date_var}/${base_name}_merged_humann_temp/${base_name}_merged_metaphlan_bugs_list.tsv  ${humann_out_dir}/humann_out_${date_var}/
 cp ${humann_out_dir}/humann_out_${date_var}/${base_name}_merged_humann_temp/${base_name}_merged.log  ${humann_out_dir}/humann_out_${date_var}/
 rm -rf ${humann_out_dir}/humann_out_${date_var}/${base_name}_merged_humann_temp;
done

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
mkdir ${humann_out_dir}/humann_merged_${date_var}
### Let's join gene family abundance outputs into a single table
humann_join_tables -i ${humann_out_dir}/humann_out_${date_var} \
  -o ${humann_merged_out_dir}/${date_var}_humann_out_genefamilies.tsv \
  --file_name genefamilies

# humann_join_tables -i ${humann_out_dir}/humann_all_${date_var} \
#   -o ${humann_merged_out_dir}/${date_var}_humann_out_genefamilies.tsv \
#   --file_name genefamilies
### Notice how we don't specify each sample but the name of directory 
### in front of _genefamilies.tsv
### file_name is also just genefamilies: only join tables with this string included in the file name
### Result: ${humann_out_dir}/humann_merged_${date_var}/${date_var}_humann_out_genefamilies.tsv
# Join pathabundance
humann_join_tables -i ${humann_out_dir}/humann_out_${date_var} \
  -o ${humann_merged_out_dir}/${date_var}_humann_out_pathabundance.tsv \
  --file_name pathabundance

# Normalizing RPKs to relative abundance
humann_renorm_table --input ${humann_merged_out_dir}/${date_var}_humann_out_genefamilies.tsv \
    --output ${humann_merged_out_dir}/${date_var}_humann_out_genefamilies-cpm.tsv \
	--units cpm \
	--update-snames
# for pathabundance renorm is enough
humann_renorm_table --input ${humann_merged_out_dir}/${date_var}_humann_out_pathabundance.tsv \
    --output ${humann_merged_out_dir}/${date_var}_humann_out_pathabundance-cpm.tsv \
	--units cpm \
	--update-snames
# Regrouping genes to other functional categories
humann_regroup_table --input ${humann_merged_out_dir}/${date_var}_humann_out_genefamilies-cpm.tsv \
    --output ${humann_merged_out_dir}/${date_var}_humann_out_genefamilies-cpm-rxn.tsv \
	--groups uniref90_rxn

# Attaching names to features
humann_rename_table --input ${humann_merged_out_dir}/${date_var}_humann_out_genefamilies-cpm-rxn.tsv \
    --output ${humann_merged_out_dir}/${date_var}_humann_out_genefamilies-cpm-rxn-named.tsv \
	--names metacyc-rxn


maaslin2_file=output/rtables/maaslin2-humann3-signif-gene_family-list-NMR.tsv
IFS=$'\r\n' GLOBIGNORE='*' command eval  'maaslin2_array=($(cat ${maaslin2_file}))'

for maaslin2_feaure in "${maaslin2_array[@]}"
do
    echo "${maaslin2_feaure}"
    humann_barplot --input ${humann_merged_out_dir}/${date_var}_humann_out_genefamilies-cpm-rxn-named_copy.tsv \
	    --focal-metadata age_group --last-metadata age_group \
      --output images/humann3_images/gene_family/gene_family_${maaslin2_feaure}.png --focal-feature ${maaslin2_feaure} \
	    --sort sum metadata
done



maaslin2_file=output/rtables/maaslin2-humann3-signif-pathabundance-NMR.tsv
IFS=$'\r\n' GLOBIGNORE='*' command eval  'maaslin2_array=($(cat ${maaslin2_file}))'
for maaslin2_feaure in "${maaslin2_array[@]}"
do
    echo "${maaslin2_feaure}"
    humann_barplot --input ${humann_merged_out_dir}/${date_var}_humann_out_pathabundance-cpm_copy.tsv \
	    --focal-metadata age_group --last-metadata age_group \
        --output images/humann3_images/pathway/pathway_${maaslin2_feaure}.png --focal-feature ${maaslin2_feaure} \
	    --sort sum metadata;
done




### TIFF
maaslin2_file=output/rtables/maaslin2-humann3-signif-gene_family-list-NMR.tsv
IFS=$'\r\n' GLOBIGNORE='*' command eval  'maaslin2_array=($(cat ${maaslin2_file}))'

for maaslin2_feaure in "${maaslin2_array[@]}"
do
    echo "${maaslin2_feaure}"
    humann_barplot --input ${humann_merged_out_dir}/${date_var}_humann_out_genefamilies-cpm-rxn-named_copy.tsv \
	    --focal-metadata age_group --last-metadata age_group \
      --output images/humann3_images/gene_family/gene_family_${maaslin2_feaure}.tiff --focal-feature ${maaslin2_feaure} \
	    --sort sum metadata;
done



maaslin2_file=output/rtables/maaslin2-humann3-signif-pathabundance-NMR.tsv
IFS=$'\r\n' GLOBIGNORE='*' command eval  'maaslin2_array=($(cat ${maaslin2_file}))'
for maaslin2_feaure in "${maaslin2_array[@]}"
do
    echo "${maaslin2_feaure}"
    humann_barplot --input ${humann_merged_out_dir}/${date_var}_humann_out_pathabundance-cpm_copy.tsv \
	    --focal-metadata age_group --last-metadata age_group \
        --output images/humann3_images/pathway/pathway_${maaslin2_feaure}.tiff --focal-feature ${maaslin2_feaure} \
	    --sort sum metadata;
done
