#!/usr/bin/env bash
#SBATCH -t 360:00:00
#SBATCH -N 4
#SBATCH -n 40
#SBATCH --mem-per-cpu 8g
#SBATCH -J 20250717_18-47-genome-annotation
#SBATCH --output jobreports/20250717_18-47-genome-annotation-out.txt
#SBATCH --error jobreports/20250717_18-47-genome-annotation-out.txt
#I am requesting 4 nodes containing 40 CPUs, with 8 GB memory per CPU. Total: 320 GB
source ~/miniconda3/etc/profile.d/conda.sh
shopt -s nullglob
# When nullglob is enabled, if a glob pattern does not match any files,
# it expands to nothing (an empty string) instead of returning the pattern itself.
# So, if no matches are found, the script will skip the file

# This script annotates contigs from MEGAHIT using several tools. Prodigal finds protein-coding features (CDS)
# (not used later). Prokka annotates MAGs using Prodigal at the start, finds more features. 
# Predicted proteins from all samples are combined into a single fasta file, and short sequences are 
# removed to avoid false-positives. Redundant sequences are removed
# using MMseqs2 to create a non-redundant gene catalog (multi-fasta file). This catalog is further annotated using
# database of carbohydrate-active enzymes (CAZymes) with dbCAN3. Additionally, this catalog is annotated using
# KEGG orthologs with KofamScan.
# Main output: Proteins predicted by PROKKA (fasta for each file), non-redundant protein sequences (multi-fasta), 
# CAZyme annotation (tab-separated file), KEGG Ortholog annotation.
date_var=$(date -I|sed 's/-//g')
time_var=$(date +%T |sed 's/:/_/g' )
date_time=${date_var}_${time_var}
start_date_time=$(date +"%F %H:%M:%S")
megahit_date_time=20250303_17_54_26
prodigal_date_time=20250303_17_54_26
# prokka_date_time=20250522_10_37_23
prokka_date_time=20250626_22_11_43
mmseqs_easy_cluster_date_time=20250626_22_11_43
dbcan_date_time=20250626_22_11_43
kofam_scan_date_time=20250718_09_21_41
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
mkdir -p output/mag_assembly/prodigal_output
mkdir -p output/mag_assembly/prokka_output
mkdir -p output/mag_assembly/mmseqs_output/easy_cluster
mkdir -p output/mag_assembly/dbcan_output
mkdir -p output/mag_assembly/kofam_scan_output


# 1. Activate conda environment
conda activate mag_assembly-tools

# 2. Prodigal
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


# # 3. Prokka
# # Input is contigs. Metagenome mode
# intermediate_date_time=$(date +"%F %H:%M:%S")
# echo "${intermediate_date_time}"
# echo "Running Prokka"
# for FILE_DIR in ${megahit_output_dir}/"${megahit_date_time}"*megahit_asm
# do 
#  SAMPLE=$(echo "${FILE_DIR}" | sed "s/\.megahit_asm//" |sed "s/${megahit_date_time}_//")
#  base_name=$(basename "$SAMPLE" )
#  prokka \
#     --outdir "${prokka_output_dir}"/"${prokka_date_time}"_"${base_name}"_prokka \
#     --prefix "${prokka_date_time}"_"${base_name}"_pred_prokka \
#     "${megahit_output_dir}"/"${megahit_date_time}"_"${base_name}".megahit_asm/"${megahit_date_time}"_"${base_name}"_final.contigs.fa \
#     --metagenome \
#     --cpus "${nthreads_sort}" 2>&1 |tee "${prokka_output_dir}"/"${prokka_date_time}"_"${base_name}"_prokka_report.txt;
# done

# # 3.1 Combine all the CDS into one file for clustering
# cat  "${prokka_output_dir}"/"${prokka_date_time}"_*_prokka/*pred_prokka.faa > "${prokka_output_dir}"/"${prokka_date_time}"_all_proteins.faa
# # There are 3771950 predicted CDS in total.
# grep ">" "${prokka_output_dir}"/"${prokka_date_time}"_all_proteins.faa | \
#  awk ' { count++ }
#     END {print "There are " count " predicted CDS in total." }
#     '  

# # 3.2 Filter proteins by length: >=100 bp (i.e. 33 aa).
# conda deactivate
# conda activate qc-tools
# seqkit seq -m 33  "${prokka_output_dir}"/"${prokka_date_time}"_all_proteins.faa > \
#  "${prokka_output_dir}"/"${prokka_date_time}"_all_proteins_filtered.faa

# # There are 3764752 predicted proteins >=100 bp (i.e. 33 aa).
# grep ">" "${prokka_output_dir}"/"${prokka_date_time}"_all_proteins_filtered.faa | \
#  awk ' { count++ }
#     END {print "There are " count " predicted proteins >=100 bp (i.e. 33 aa)." }
#     ' 

# conda deactivate
# conda activate mag_assembly-tools

# # 4. Remove redundant proteins with MMseqs2
# intermediate_date_time=$(date +"%F %H:%M:%S")
# echo "${intermediate_date_time}"
# echo "Running MMseqs2"
# # mmseqs module input_db output_db args [options]
# # cluster the proteins with the criteria of identity ≥95% and overlap ≥90%. 
# # Need to create output directory first
# mkdir -p "${mmseqs_easy_cluster_output_dir}"/"${mmseqs_easy_cluster_date_time}"
# mmseqs easy-cluster "${prokka_output_dir}"/"${prokka_date_time}"_all_proteins_filtered.faa \
#  "${mmseqs_easy_cluster_output_dir}"/"${mmseqs_easy_cluster_date_time}"/"${mmseqs_easy_cluster_date_time}"_prokka_nr_prot \
#  "${mmseqs_easy_cluster_output_dir}"/"${mmseqs_easy_cluster_date_time}"/"${mmseqs_easy_cluster_date_time}"_prokka_nr_prot_temp \
#  --min-seq-id 0.95 \
#  -c 0.9 \
#  --cov-mode 0 \
#  --threads "${nthreads_sort}" 2>&1 |tee "${mmseqs_easy_cluster_output_dir}"/"${mmseqs_easy_cluster_date_time}"_mmseqs_easy_cluster_report.txt
# intermediate_date_time=$(date +"%F %H:%M:%S")
# echo "${intermediate_date_time}"

# # 4.1 There are 1922857 nonredundant predicted proteins in total.
# grep ">" "${mmseqs_easy_cluster_output_dir}"/"${mmseqs_easy_cluster_date_time}"/"${mmseqs_easy_cluster_date_time}"_prokka_nr_prot_rep_seq.fasta | \
#  awk ' { count++ }
#     END {print "There are " count " nonredundant predicted proteins in total." }
#     ' 

# # 5. Use dbCAN to predict CAZymes
# # conda create --name dbcan-tools -c conda-forge -c bioconda python=3.8 dbcan 
# conda deactivate
# conda activate dbcan-tools
# # Download the database (7.9 GB)
# # dbcan_build --cpus "${nthreads}" --db-dir data/dbcan_db --clean
# # run_dbcan -h

# echo "Running dbCAN on nonredundant gene catalog"
# intermediate_date_time=$(date +"%F %H:%M:%S")
# echo "${intermediate_date_time}"
# run_dbcan "${mmseqs_easy_cluster_output_dir}"/"${mmseqs_easy_cluster_date_time}"/"${mmseqs_easy_cluster_date_time}"_prokka_nr_prot_rep_seq.fasta \
#    protein \
#    --out_dir "${dbcan_output_dir}"/"${dbcan_date_time}"_nr_dbcan_proteins \
#    --db_dir data/dbcan_db \
#    --dia_cpu "${nthreads}" \
#    --hmm_cpu "${nthreads}" \
#    --stp_cpu "${nthreads}" \
#    2>&1 |tee "${dbcan_output_dir}"/"${dbcan_date_time}"_nr_dbcan_report.txt
# Rename the output file for comprehension
# cp "${dbcan_output_dir}"/"${dbcan_date_time}"_nr_dbcan_proteins/overview.txt \
#   "${dbcan_output_dir}"/"${dbcan_date_time}"_nr_dbcan_proteins/"${dbcan_date_time}"_nr_dbcan_overview.txt 

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
conda deactivate 

# 6. Run KofamScan to predict KEGG orthologs
# conda activate mag_assembly-tools
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
echo "Running KofamScan"
mkdir -p "${kofam_scan_output_dir}"/"${kofam_scan_date_time}"_tmp_dir
~/kofamscan/bin/kofam_scan-1.3.0/exec_annotation \
  -o "${kofam_scan_output_dir}"/"${kofam_scan_date_time}"_nr_prot_kofam_scan.txt \
  -p ~/kofamscan/db/profiles \
  -k ~/kofamscan/db/ko_list \
  --tmp-dir "${kofam_scan_output_dir}"/"${kofam_scan_date_time}"_tmp_dir \
  --cpu "${nthreads}" \
  -f detail-tsv \
  "${mmseqs_easy_cluster_output_dir}"/"${mmseqs_easy_cluster_date_time}"/"${mmseqs_easy_cluster_date_time}"_prokka_nr_prot_rep_seq.fasta
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"

# The first field has asterisks, so we don't use it.
# Remove entries where threshold is empty or score is less than threshold. Then, sort by score and evalue
awk 'BEGIN {FS=OFS="\t"} 
  FNR>2{
    gene = $2;
    ko = $3;
    threshold = $4;
    score = $5;
    evalue = $6;
    ko_def = $7;
    if ((threshold !="") && (score > threshold)) {print gene, ko, threshold, score, evalue, ko_def}
    
}' "${kofam_scan_output_dir}"/"${kofam_scan_date_time}"_nr_prot_kofam_scan.txt | \
  sort -t$'\t' -k1,1 -k4,4nr -k5,5g | awk '!seen[$1]++' | \
  awk 'BEGIN {FS=OFS="\t"; print "gene_name","ko","threshold","score","evalue","ko_definition"} {print}' > \
   "${kofam_scan_output_dir}"/"${kofam_scan_date_time}"_nr_prot_kofam_scan_top_hits.txt


# #######
# # awk -v ORS="" 'FNR>1 {print  "-p", "\x27"$3, "\x27 "}' output/rtables/gene-annotation-eukaryotes-rrna.tsv |cat -A
# # seqkit grep -n -r -p $(awk -v ORS="" 'FNR>1 {print  "-p", "\x27"$3, "\x27 "}' output/rtables/gene-annotation-eukaryotes-rrna.tsv)  output/rtables/gene-annotation-eukaryotes-cds-ids.txt 

# # Get protein sequences
# awk 'FNR>1 {print $3, ""}' output/rtables/gene-annotation-eukaryotes-cds.tsv  > \
#   output/rtables/gene-annotation-eukaryotes-cds-ids.txt
# seqkit grep -n -r -f output/rtables/gene-annotation-eukaryotes-cds-ids.txt \
#   output/mag_assembly/prokka_output/20250626_22_11_43_all_proteins.faa  > \
#   output/fasta/eukaryotic-cds-pep.faa

# # Get gene sequences
# awk '$1=="2D10" {print $0}' output/rtables/gene-annotation-eukaryotes-rrna.tsv |awk '{print $3,""}' > \
#   output/rtables/2D10-rrna-ids.txt

# awk '$1=="2D14" {print $0}' output/rtables/gene-annotation-eukaryotes-rrna.tsv |awk '{print $3,""}' > \
#   output/rtables/2D14-rrna-ids.txt

# awk '$1=="G14" {print $0}' output/rtables/gene-annotation-eukaryotes-rrna.tsv |awk '{print $3,""}' > \
#   output/rtables/G14-rrna-ids.txt

# awk '$1=="G18" {print $0}' output/rtables/gene-annotation-eukaryotes-rrna.tsv |awk '{print $3,""}' > \
#   output/rtables/G18-rrna-ids.txt

# awk '$1=="H15" {print $0}' output/rtables/gene-annotation-eukaryotes-rrna.tsv |awk '{print $3,""}' > \
#   output/rtables/H15-rrna-ids.txt

# awk '$1=="H3" {print $0}' output/rtables/gene-annotation-eukaryotes-rrna.tsv |awk '{print $3,""}' > \
#   output/rtables/H3-rrna-ids.txt

# awk '$1=="H4" {print $0}' output/rtables/gene-annotation-eukaryotes-rrna.tsv |awk '{print $3,""}' > \
#   output/rtables/H4-rrna-ids.txt

# awk '$1=="H21" {print $0}' output/rtables/gene-annotation-eukaryotes-rrna.tsv |awk '{print $3,""}' > \
#   output/rtables/H21-rrna-ids.txt

# awk '$1=="O15" {print $0}' output/rtables/gene-annotation-eukaryotes-rrna.tsv |awk '{print $3,""}' > \
#   output/rtables/O15-rrna-ids.txt

# awk '$1=="Y51b" {print $0}' output/rtables/gene-annotation-eukaryotes-rrna.tsv |awk '{print $3,""}' > \
#   output/rtables/Y51b-rrna-ids.txt

# awk '$1=="Y66b" {print $0}' output/rtables/gene-annotation-eukaryotes-rrna.tsv |awk '{print $3,""}' > \
#   output/rtables/Y66b-rrna-ids.txt


# touch output/fasta/eukaryotic-rrna-gene.fna
# for fname in output/rtables/*-rrna-ids.txt
# do 
#   SAMPLE=$(basename $fname |sed "s/-rrna-ids\.txt//")
#   seqkit grep -n -r -f $fname \
#   output/mag_assembly/prokka_output/20250626_22_11_43_${SAMPLE}_prokka/20250626_22_11_43_${SAMPLE}_pred_prokka.ffn >> \
#   output/fasta/eukaryotic-rrna-gene.fna;
# done

# # Get contigs
# for FILE_DIR in output/mag_assembly/prokka_output/20250626_22_11_43_*_prokka
# do
#   SAMPLE=$(basename $FILE_DIR | sed "s/20250626_22_11_43_//"| sed "s/_prokka//")
#   grep "Eukaryota" output/mag_assembly/blastn_output/20250701_05_24_10/20250701_05_24_10_blastn_all_samples.tsv | \
#     awk -v sample_name=${SAMPLE} '$1 == sample_name {print $2}' > output/rtables/${SAMPLE}_euk_contig_ids.txt
#   seqkit grep -n  -f output/rtables/${SAMPLE}_euk_contig_ids.txt \
#     output/mag_assembly/prokka_output/20250626_22_11_43_${SAMPLE}_prokka/20250626_22_11_43_${SAMPLE}_pred_prokka.fna > \
#     output/fasta/${SAMPLE}_euk_contigs.fna;
# done
# cat output/fasta/*_euk_contigs.fna > output/fasta/all_euk_contigs_merged.fna
