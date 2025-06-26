#!/usr/bin/env bash
#SBATCH -t 360:00:00
#SBATCH -N 4
#SBATCH -n 40
#SBATCH --mem-per-cpu 8g
#SBATCH -J 20250625_18-44-genome-annotation
#SBATCH --output jobreports/20250625_18-44-genome-annotation-out.txt
#SBATCH --error jobreports/20250625_18-44-genome-annotation-out.txt
#I am requesting 4 nodes containing 40 CPUs, with 8 GB memory per CPU. Total: 320 GB
source ~/miniconda3/etc/profile.d/conda.sh
shopt -s nullglob
# When nullglob is enabled, if a glob pattern does not match any files,
# it expands to nothing (an empty string) instead of returning the pattern itself.
# So, if no matches are found, the script will skip the file

# Prodigal finds protein-coding features (CDS)
# Prokka annotates MAGs using Prodigal at the start. Finds more features
date_var=$(date -I|sed 's/-//g')
time_var=$(date +%T |sed 's/:/_/g' )
date_time=${date_var}_${time_var}
start_date_time=$(date +"%F %H:%M:%S")
megahit_date_time=20250303_17_54_26
prodigal_date_time=20250303_17_54_26
# prokka_date_time=20250522_10_37_23
prokka_date_time="${date_time}"
mmseqs_easy_cluster_date_time="${date_time}"
dbcan_date_time="${date_time}"
kofam_scan_date_time="${date_time}"
nthreads=40
nthreads_sort=35
mem_req=8G
mem_req_sort=4G
project_home_dir=~/projects/metagenome
output_dir=~/projects/metagenome/output/mag_assembly
megahit_output_dir=output/mag_assembly/megahit_output
prodigal_output_dir=output/mag_assembly/prodigal_output
prokka_output_dir=output/mag_assembly/prokka_output
mmseqs_easy_cluster_output_dir=output/mag_assembly/mmseqs_output/easy_cluster
dbcan_output_dir=output/mag_assembly/dbcan_output
kofam_scan_output_dir=output/mag_assembly/kofam_scan_output


cd ${project_home_dir}
mkdir -p data/dbcan_db
mkdir -p output/mag_assembly/megahit_output/alignedreads
mkdir -p output/mag_assembly/metabat2_reports
mkdir -p output/mag_assembly/samtools_reports
mkdir -p output/mag_assembly/checkm2_output
mkdir -p output/mag_assembly/gtdbtk_output
mkdir -p output/mag_assembly/bam_contig_depths
mkdir -p output/mag_assembly/metaquast_output
mkdir -p output/mag_assembly/prodigal_output
mkdir -p output/mag_assembly/prokka_output
mkdir -p output/mag_assembly/contig_coverages
mkdir -p output/mag_assembly/mmseqs_output/easy_cluster
mkdir -p output/mag_assembly/dbcan_output
mkdir -p output/mag_assembly/kofam_scan_output


# 1. Activate conda environment
conda activate mag_assembly-tools

# 4. Prodigal
# Anonymous mode (for metagenomes)
# should I use -p meta or -p anon?
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
# for FILE_DIR in ${megahit_output_dir}/"${megahit_date_time}"*megahit_asm
# do 
#  SAMPLE=$(echo "${FILE_DIR}" | sed "s/\.megahit_asm//" |sed "s/${megahit_date_time}_//")
#  base_name=$(basename "$SAMPLE" )
#  mkdir "${prodigal_output_dir}"/"${prodigal_date_time}"_"${base_name}"_prodigal
#  prodigal \
#     -i "${megahit_output_dir}"/"${megahit_date_time}"_"${base_name}".megahit_asm/"${megahit_date_time}"_"${base_name}"_final.contigs.fa \
#     -o "${prodigal_output_dir}"/"${prodigal_date_time}"_"${base_name}"_prodigal/"${prodigal_date_time}"_"${base_name}"_prodigal_coords.gbk \
#     -a "${prodigal_output_dir}"/"${prodigal_date_time}"_"${base_name}"_prodigal/"${prodigal_date_time}"_"${base_name}"_prodigal_proteins.faa \
#     -p meta 2>&1 |tee "${prodigal_output_dir}"/"${prodigal_date_time}"_"${base_name}"_prodigal_report.txt;
# done


# Prokka
# Input is contigs. Metagenome mode
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
echo "Running Prokka"
for FILE_DIR in ${megahit_output_dir}/"${megahit_date_time}"*megahit_asm
do 
 SAMPLE=$(echo "${FILE_DIR}" | sed "s/\.megahit_asm//" |sed "s/${megahit_date_time}_//")
 base_name=$(basename "$SAMPLE" )
 prokka \
    --outdir "${prokka_output_dir}"/"${prokka_date_time}"_"${base_name}"_prokka \
    --prefix "${prokka_date_time}"_"${base_name}"_pred_prokka \
    "${megahit_output_dir}"/"${megahit_date_time}"_"${base_name}".megahit_asm/"${megahit_date_time}"_"${base_name}"_final.contigs.fa \
    --metagenome \
    --cpus "${nthreads_sort}" 2>&1 |tee "${prokka_output_dir}"/"${prokka_date_time}"_"${base_name}"_prokka_report.txt;
done

# Combine all the proteins into one file for clustering
cat  "${prokka_output_dir}"/"${prokka_date_time}"_*_prokka/*pred_prokka.faa > "${prokka_output_dir}"/"${prokka_date_time}"_all_proteins.faa
# There are 3871810 predicted proteins in total.
grep ">" "${prokka_output_dir}"/"${prokka_date_time}"_all_proteins.faa | \
 awk ' { count++ }
    END {print "There are " count " predicted proteins in total." }
    '  

# Filter proteins by length
conda deactivate
conda activate qc-tools
seqkit seq -m 33  "${prokka_output_dir}"/"${prokka_date_time}"_all_proteins.faa > \
 "${prokka_output_dir}"/"${prokka_date_time}"_all_proteins_filtered.faa

# There are 3852959 predicted proteins >=100 bp (i.e. 33 aa).
grep ">" "${prokka_output_dir}"/"${prokka_date_time}"_all_proteins_filtered.faa | \
 awk ' { count++ }
    END {print "There are " count " predicted proteins >=100 bp (i.e. 33 aa)." }
    ' 

conda deactivate
conda activate mag_assembly-tools

# Remove redundant genes with MMseqs2
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
echo "Running MMseqs2"
# mmseqs module input_db output_db args [options]
# cluster the genes with the criteria of identity ≥95% and overlap ≥90%. 
# Need to create output directory first
mkdir -p "${mmseqs_easy_cluster_output_dir}"/"${mmseqs_easy_cluster_date_time}"
mmseqs easy-cluster "${prokka_output_dir}"/"${prokka_date_time}"_all_proteins_filtered.faa \
 "${mmseqs_easy_cluster_output_dir}"/"${mmseqs_easy_cluster_date_time}"/"${mmseqs_easy_cluster_date_time}"_prokka_nr_prot \
 "${mmseqs_easy_cluster_output_dir}"/"${mmseqs_easy_cluster_date_time}"/"${mmseqs_easy_cluster_date_time}"_prokka_nr_prot_temp \
 --min-seq-id 0.95 \
 -c 0.9 \
 --cov-mode 0 \
 --threads "${nthreads_sort}" 2>&1 |tee "${mmseqs_easy_cluster_output_dir}"/"${mmseqs_easy_cluster_date_time}"_mmseqs_easy_cluster_report.txt
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"

# There are 1951683 nonredundant predicted proteins in total.
grep ">" "${mmseqs_easy_cluster_output_dir}"/"${mmseqs_easy_cluster_date_time}"/"${mmseqs_easy_cluster_date_time}"_prokka_nr_prot_rep_seq.fasta | \
 awk ' { count++ }
    END {print "There are " count " nonredundant predicted proteins in total." }
    ' 

# Use dbCAN to predict CAZymes
# conda create --name dbcan-tools -c conda-forge -c bioconda python=3.8 dbcan 
conda deactivate
conda activate dbcan-tools
# Download the database (7.9 GB)
# dbcan_build --cpus "${nthreads}" --db-dir data/dbcan_db --clean
# run_dbcan -h

echo "Running dbCAN on nonredundant gene catalog"
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
run_dbcan "${mmseqs_easy_cluster_output_dir}"/"${mmseqs_easy_cluster_date_time}"/"${mmseqs_easy_cluster_date_time}"_prokka_nr_prot_rep_seq.fasta \
   protein \
   --out_dir "${dbcan_output_dir}"/"${dbcan_date_time}"_nr_dbcan_proteins \
   --db_dir data/dbcan_db \
   --dia_cpu "${nthreads}" \
   --hmm_cpu "${nthreads}" \
   --stp_cpu "${nthreads}" \
   2>&1 |tee "${dbcan_output_dir}"/"${dbcan_date_time}"_nr_dbcan_report.txt

# for FILE_DIR in "${prodigal_output_dir}"/"${prodigal_date_time}"*_prodigal
# do 
#  SAMPLE=$(echo "${FILE_DIR}" | sed "s/_prodigal//" |sed "s/${prodigal_date_time}_//")
#  base_name=$(basename "$SAMPLE" )
#  echo "Running dbCAN on ${base_name}"
#  intermediate_date_time=$(date +"%F %H:%M:%S")
#  echo "${intermediate_date_time}"
#  run_dbcan "${prodigal_output_dir}"/"${prodigal_date_time}"_"${base_name}"_prodigal/"${prodigal_date_time}"_"${base_name}"_prodigal_proteins.faa \
#    protein \
#    --out_dir "${dbcan_output_dir}"/"${date_time}"_"${base_name}"_dbcan_proteins \
#    --db_dir data/dbcan_db \
#    --dia_cpu "${nthreads}" \
#    --hmm_cpu "${nthreads}" \
#    --stp_cpu "${nthreads}" \
#    2>&1 |tee "${dbcan_output_dir}"/"${date_time}"_"${base_name}"_dbcan_report.txt;
# done

conda activate mag_assembly-tools
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
echo "Running KofamScan"
~/kofam_scan/exec_annotation \
  --cpu "${nthreads}"
  -o "${kofam_scan_output_dir}"/"${kofam_scan_date_time}"_nr_prot_kofam_scan.txt \
  "${mmseqs_easy_cluster_output_dir}"/"${mmseqs_easy_cluster_date_time}"/"${mmseqs_easy_cluster_date_time}"_prokka_nr_prot_rep_seq.fasta
