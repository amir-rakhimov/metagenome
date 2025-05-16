#!/usr/bin/env bash
#SBATCH -t 240:00:00
#SBATCH -N 4
#SBATCH -n 10
#SBATCH --mem-per-cpu 8g
#SBATCH -J 20250401_16-07-dbcan
#SBATCH --output jobreports/20250401_16-33-dbcan-out.txt
#SBATCH --error jobreports/20250401_16-33-dbcan-out.txt
# I am requesting 4 nodes containing 10 CPUs, with 8 GB memory per CPU. Total: 80 GB

source ~/miniconda3/etc/profile.d/conda.sh
shopt -s nullglob
# When nullglob is enabled, if a glob pattern does not match any files,
# it expands to nothing (an empty string) instead of returning the pattern itself.
# So, if no matches are found, the script will skip the file

date_var=$(date -I|sed 's/-//g')
time_var=$(date +%T |sed 's/:/_/g' )
date_time=${date_var}_${time_var}
start_date_time=$(date +"%F %H:%M:%S")
# For binning with MetaBAT2
prodigal_date_time=20250303_17_54_26
nthreads=10
mem_req=8G
project_home_dir=~/projects/metagenome
output_dir=~/projects/metagenome/output/mag_assembly
prodigal_output_dir=output/mag_assembly/prodigal_output
dbcan_output_dir=output/mag_assembly/dbcan_output

cd ${project_home_dir}
mkdir data/dbcan_db
mkdir output/mag_assembly/dbcan_output

# 0. Show the current time for logging
echo "${start_date_time}"

# Use dbCAN to predict CAZymes
# conda create --name dbcan-tools -c conda-forge -c bioconda python=3.8 dbcan 
conda activate dbcan-tools
# Download the database (7.9 GB)
# dbcan_build --cpus "${nthreads}" --db-dir data/dbcan_db --clean
# run_dbcan -h

for FILE_DIR in "${prodigal_output_dir}"/"${prodigal_date_time}"*_prodigal
do 
 SAMPLE=$(echo "${FILE_DIR}" | sed "s/_prodigal//" |sed "s/${prodigal_date_time}_//")
 base_name=$(basename "$SAMPLE" )
 echo "Running dbCAN on ${base_name}"
 intermediate_date_time=$(date +"%F %H:%M:%S")
 echo "${intermediate_date_time}"
 run_dbcan "${prodigal_output_dir}"/"${prodigal_date_time}"_"${base_name}"_prodigal/"${prodigal_date_time}"_"${base_name}"_prodigal_proteins.faa \
   protein \
   --out_dir "${dbcan_output_dir}"/"${date_time}"_"${base_name}"_dbcan_proteins \
   --db_dir data/dbcan_db \
   --dia_cpu "${nthreads}" \
   --hmm_cpu "${nthreads}" \
   --stp_cpu "${nthreads}" \
   2>&1 |tee "${dbcan_output_dir}"/"${date_time}"_"${base_name}"_dbcan_report.txt;
done
end_date_time=$(date +"%F %H:%M:%S")
echo "${end_date_time}"



# Use dbCAN to predict CAZymes
conda create --name dbcan-tools -c conda-forge -c bioconda python=3.8 dbcan 
conda activate dbcan-tools
mkdir data/dbcan_db
mkdir output/mag_assembly/dbcan_output
# Download the database (7.9 GB)
dbcan_build --cpus "${nthreads}" --db-dir data/dbcan_db --clean
run_dbcan -h
date_var=$(date -I|sed 's/-//g')
time_var=$(date +%T |sed 's/:/_/g' )
date_time=${date_var}_${time_var}
dbcan_output_dir=output/mag_assembly/dbcan_output
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
for FILE_DIR in "${prodigal_output_dir}"/"${prodigal_date_time}"*_prodigal
do 
 SAMPLE=$(echo "${FILE_DIR}" | sed "s/_prodigal//" |sed "s/${prodigal_date_time}_//")
 base_name=$(basename "$SAMPLE" )
 echo "Running dbCAN on ${base_name}"
 run_dbcan "${prodigal_output_dir}"/"${prodigal_date_time}"_"${base_name}"_prodigal/"${prodigal_date_time}"_"${base_name}"_prodigal_proteins.faa \
   --out_dir "${dbcan_output_dir}"/"${date_time}"_"${base_name}"_dbcan_proteins \
   --db_dir data/dbcan_db \
   --dia_cpu "${nthreads}" \
   --hmm_cpu "${nthreads}" \
   --stp_cpu "${nthreads}" \
   2>&1 |tee "${dbcan_output_dir}"/"${date_time}"_"${base_name}"_dbcan_report.txt;
done


run_dbcan output/mag_assembly/prodigal_output/20250303_17_54_26_2D10_prodigal/20250303_17_54_26_2D10_prodigal_proteins.faa \
 protein \
 --out_dir output/mag_assembly/dbcan_output/"${date_time}"_2D10_dbcan_proteins \
 --db_dir data/dbcan_db \
 --dia_cpu "${nthreads}" \
 --hmm_cpu "${nthreads}" \
 --stp_cpu "${nthreads}"

end_date_time=$(date +"%F %H:%M:%S")
echo "${end_date_time}"