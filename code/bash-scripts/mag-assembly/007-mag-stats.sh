#!/usr/bin/env bash
#SBATCH -t 360:00:00
#SBATCH -N 4
#SBATCH -n 45
#SBATCH --mem-per-cpu 8g
#SBATCH -J 20250730_18-04-mag-stats
#SBATCH --output jobreports/20250730_18-04-mag-stats-out.txt
#SBATCH --error jobreports/20250730_18-04-mag-stats-out.txt
#I am requesting 4 nodes containing 45 CPUs, with 8 GB memory per CPU. Total: 360 GB
source ~/miniconda3/etc/profile.d/conda.sh
shopt -s nullglob

# This script calculates statistics of assembled bins from MetaBAT2. CoverM is used to 
# calculate relative abundances of dereplicated bins and contigs (from MEGAHIT).
# Seqkit is used to calculate statistics such as contig length, GC%, AT% for all contigs 
# in MEGAHIT output and in dereplicated bins, and saves as tab-separated files for each sample.
date_var=$(date -I|sed 's/-//g')
time_var=$(date +%T |sed 's/:/_/g' )
date_time=${date_var}_${time_var}
# metabat2_date_time=20250417_23_21_29
metabat2_date_time=20250612_13_37_47
drep_date_time=20250627_19_56_42
megahit_date_time=20250303_17_54_26
seqkit_stats_date_time=20250712_18_07_37
coverm_date_time=20250803_05_07_53
coverm_drep_date_time=20250803_05_07_53
coverm_all_bins_date_time=20250803_05_07_53
coverm_contigs_date_time=20250803_05_07_53
nthreads=40
nthreads_sort=35
mem_req=8G
mem_req_sort=4G
project_home_dir=~/projects/metagenome
megahit_output_dir=output/mag_assembly/megahit_output
bowtie2_decontam_fastq_dir=data/bowtie2_decontam_fastq
output_dir=~/projects/metagenome/output/mag_assembly
metabat2_output_dir=output/mag_assembly/metabat2_output
drep_output_dir=output/mag_assembly/drep_output
seqkit_output_dir=output/mag_assembly/seqkit_output
coverm_output_dir=output/mag_assembly/coverm_output

cd ${project_home_dir}
mkdir -p output/mag_assembly/seqkit_output
mkdir -p output/mag_assembly/coverm_output

# 0. Activate environment
conda activate coverm-env

intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"

# 1. CoverM calculates abundance of dereplicated bins
for FILE in "${bowtie2_decontam_fastq_dir}"/*decontam_R1.fastq.gz
do 
 SAMPLE=$(basename "${FILE}" | sed "s/_trim_decontam_R1\.fastq\.gz//")
 #  base_name=$(basename "$SAMPLE" )
 ## At species level (95% ANI)
 intermediate_date_time=$(date +"%F %H:%M:%S")
 echo "${intermediate_date_time}"
 echo "Running CoverM on ${SAMPLE} at species level (95% ANI)"
 mkdir -p "${coverm_output_dir}"/"${coverm_drep_date_time}"_"${SAMPLE}"_drep95_tmp
 TMPDIR="${coverm_output_dir}"/"${coverm_drep_date_time}"_"${SAMPLE}"_drep95_tmp
 coverm genome \
  --coupled "${bowtie2_decontam_fastq_dir}"/"${SAMPLE}"_trim_decontam_R1.fastq.gz \
    "${bowtie2_decontam_fastq_dir}"/"${SAMPLE}"_trim_decontam_R2.fastq.gz \
  --genome-fasta-files \
  "${drep_output_dir}"/"${drep_date_time}"_sa_95perc/dereplicated_genomes/"${metabat2_date_time}"_"${SAMPLE}"_bin.*.fa \
  -t "${nthreads}" \
  --min-covered-fraction 0 \
  -m relative_abundance mean trimmed_mean covered_fraction covered_bases variance length count \
    reads_per_base rpkm tpm \
  -o "${coverm_output_dir}"/"${coverm_drep_date_time}"_"${SAMPLE}"_sa_95perc_coverm_output.tsv 
 intermediate_date_time=$(date +"%F %H:%M:%S")
 echo "${intermediate_date_time}"
 ## At strain level (99% ANI)
 mkdir -p "${coverm_output_dir}"/"${coverm_drep_date_time}"_"${SAMPLE}"_drep99_tmp
 TMPDIR="${coverm_output_dir}"/"${coverm_drep_date_time}"_"${SAMPLE}"_drep99_tmp
 echo "Running CoverM on ${SAMPLE} at strain level (99% ANI)"
 coverm genome \
  --coupled "${bowtie2_decontam_fastq_dir}"/"${SAMPLE}"_trim_decontam_R1.fastq.gz \
    "${bowtie2_decontam_fastq_dir}"/"${SAMPLE}"_trim_decontam_R2.fastq.gz \
  --genome-fasta-files \
  "${drep_output_dir}"/"${drep_date_time}"_sa_99perc/dereplicated_genomes/"${metabat2_date_time}"_"${SAMPLE}"_bin.*.fa \
  -t "${nthreads}" \
  --min-covered-fraction 0 \
  -m relative_abundance mean trimmed_mean covered_fraction covered_bases variance length count \
    reads_per_base rpkm tpm \
  -o "${coverm_output_dir}"/"${coverm_drep_date_time}"_"${SAMPLE}"_sa_99perc_coverm_output.tsv 
 intermediate_date_time=$(date +"%F %H:%M:%S")
 echo "${intermediate_date_time}"
 # Running CoverM on all bins
 mkdir -p "${coverm_output_dir}"/"${coverm_all_bins_date_time}"_"${SAMPLE}"_all_bins_tmp
 TMPDIR="${coverm_output_dir}"/"${coverm_all_bins_date_time}"_"${SAMPLE}"_all_bins_tmp
 echo "Running CoverM on all bins from ${SAMPLE}"
 coverm genome \
  --coupled "${bowtie2_decontam_fastq_dir}"/"${SAMPLE}"_trim_decontam_R1.fastq.gz \
    "${bowtie2_decontam_fastq_dir}"/"${SAMPLE}"_trim_decontam_R2.fastq.gz \
  --genome-fasta-files \
  "${metabat2_output_dir}"/"${metabat2_date_time}"_"${SAMPLE}"_bins/"${metabat2_date_time}"_"${SAMPLE}"_bin.*.fa \
  -t "${nthreads}" \
  --min-covered-fraction 0 \
  -m relative_abundance mean trimmed_mean covered_fraction covered_bases variance length count \
    reads_per_base rpkm tpm \
  -o "${coverm_output_dir}"/"${coverm_all_bins_date_time}"_"${SAMPLE}"_all_bins_coverm_output.tsv;
done 2>&1 |tee "${coverm_output_dir}"/"${coverm_date_time}"_all_coverm_report.txt
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"


for FILE in "${bowtie2_decontam_fastq_dir}"/*decontam_R1.fastq.gz
do 
 SAMPLE=$(basename "${FILE}" | sed "s/_trim_decontam_R1\.fastq\.gz//")
 intermediate_date_time=$(date +"%F %H:%M:%S")
 echo "${intermediate_date_time}"
 # On contigs
 mkdir -p "${coverm_output_dir}"/"${coverm_contigs_date_time}"_"${SAMPLE}"_contigs_tmp
 TMPDIR="${coverm_output_dir}"/"${coverm_contigs_date_time}"_"${SAMPLE}"_contigs_tmp
 echo "Running CoverM on ${SAMPLE} at contig level"
 coverm contig \
  --coupled "${bowtie2_decontam_fastq_dir}"/"${SAMPLE}"_trim_decontam_R1.fastq.gz \
    "${bowtie2_decontam_fastq_dir}"/"${SAMPLE}"_trim_decontam_R2.fastq.gz \
  --reference \
  "${megahit_output_dir}"/"${megahit_date_time}"_"${SAMPLE}".megahit_asm/"${megahit_date_time}"_"${SAMPLE}"_final.contigs.fa \
  -t "${nthreads}" \
  --min-covered-fraction 0 \
  -m mean trimmed_mean covered_fraction covered_bases variance length count reads_per_base rpkm tpm \
  -o "${coverm_output_dir}"/"${coverm_contigs_date_time}"_"${SAMPLE}"_megahit_coverm_output.tsv;
done 2>&1 |tee "${coverm_output_dir}"/"${coverm_date_time}"_all_coverm_contig_report.txt



# Fix the unmapped value to merge
for FILE in "${coverm_output_dir}"/"${coverm_drep_date_time}"_*_sa_95perc_coverm_output.tsv
do
 SAMPLE=$(basename "${FILE}" | sed "s/${coverm_drep_date_time}_//" |sed "s/_sa_95perc_coverm_output\.tsv//")
 awk -F'\t' -v sample_name="${SAMPLE}" -v OFS='\t' '{
  if ($1=="unmapped") {$1= "unmapped_" sample_name;}
      {print }
  }'  \
  "${coverm_output_dir}"/"${coverm_drep_date_time}"_"${SAMPLE}"_sa_95perc_coverm_output.tsv > \
  "${coverm_output_dir}"/"${coverm_drep_date_time}"_"${SAMPLE}"_sa_95perc_coverm_output_unmap_fix.tsv 
 # Fix for the other file (99%)
 awk -F'\t' -v sample_name="${SAMPLE}" -v OFS='\t' '{
  if ($1=="unmapped") {$1= "unmapped_" sample_name;} 
      {print }
  }'  \
  "${coverm_output_dir}"/"${coverm_drep_date_time}"_"${SAMPLE}"_sa_99perc_coverm_output.tsv > \
  "${coverm_output_dir}"/"${coverm_drep_date_time}"_"${SAMPLE}"_sa_99perc_coverm_output_unmap_fix.tsv;
done

# For all bins
for FILE in "${coverm_output_dir}"/"${coverm_all_bins_date_time}"_*all_bins_coverm_output.tsv
do
 SAMPLE=$(basename "${FILE}" | sed "s/${coverm_all_bins_date_time}_//" |sed "s/_all_bins_coverm_output\.tsv//")
 awk -F'\t' -v sample_name="${SAMPLE}" -v OFS='\t' '{
  if ($1=="unmapped") {$1= "unmapped_" sample_name;}
      {print }
  }'  \
  "${coverm_output_dir}"/"${coverm_all_bins_date_time}"_"${SAMPLE}"_all_bins_coverm_output.tsv > \
  "${coverm_output_dir}"/"${coverm_all_bins_date_time}"_"${SAMPLE}"_all_bins_coverm_output_unmap_fix.tsv;
done

# For all contigs: we need to add sample name to contig
for FILE in "${coverm_output_dir}"/"${coverm_contigs_date_time}"_*_megahit_coverm_output.tsv
do
 SAMPLE=$(basename "${FILE}" | sed "s/${coverm_contigs_date_time}_//" |sed "s/_megahit_coverm_output\.tsv//")
 awk -F'\t' -v sample_name="${SAMPLE}" -v OFS='\t' '{
  $1= sample_name "_" $1 ;
  print
  }'  \
  "${coverm_output_dir}"/"${coverm_contigs_date_time}"_"${SAMPLE}"_megahit_coverm_output.tsv > \
  "${coverm_output_dir}"/"${coverm_contigs_date_time}"_"${SAMPLE}"_megahit_coverm_output_fix.tsv;
done

conda activate qc-tools
# Merge all 95% ANI bins
awk 'FNR!=1'  "${coverm_output_dir}"/"${coverm_drep_date_time}"_*_sa_95perc_coverm_output_unmap_fix.tsv | \
  csvtk -t add-header -n genome,relative_abundance,mean,trimmed_mean,covered_fraction,covered_bases,variance,length,read_count,reads_per_base,rpkm,tpm > \
  "${coverm_output_dir}"/"${coverm_drep_date_time}"_all_samples_sa_95perc_coverm.tsv
# Merge all 99% ANI bins
awk 'FNR!=1'  "${coverm_output_dir}"/"${coverm_drep_date_time}"_*_sa_99perc_coverm_output_unmap_fix.tsv | \
  csvtk -t add-header -n genome,relative_abundance,mean,trimmed_mean,covered_fraction,covered_bases,variance,length,read_count,reads_per_base,rpkm,tpm > \
  "${coverm_output_dir}"/"${coverm_drep_date_time}"_all_samples_sa_99perc_coverm.tsv
# Merge all metabat bins
awk 'FNR!=1'  "${coverm_output_dir}"/"${coverm_all_bins_date_time}"_*_all_bins_coverm_output_unmap_fix.tsv | \
  csvtk -t add-header -n genome,relative_abundance,mean,trimmed_mean,covered_fraction,covered_bases,variance,length,read_count,reads_per_base,rpkm,tpm > \
  "${coverm_output_dir}"/"${coverm_all_bins_date_time}"_all_samples_bins_coverm.tsv
# Merge all megahit contigs
awk 'FNR!=1'  "${coverm_output_dir}"/"${coverm_contigs_date_time}"_*_megahit_coverm_output_fix.tsv | \
  csvtk -t add-header -n contig,mean,trimmed_mean,covered_fraction,covered_bases,variance,length,read_count,reads_per_base,rpkm,tpm > \
  "${coverm_output_dir}"/"${coverm_contigs_date_time}"_all_samples_megahit_coverm.tsv

# conda activate qc-tools
# # 2. Calculate statistics for each contig in each sample (MEGAHIT output): contig IDs, contig lengths, GC %, 
# # AT %. Output is a tab-separated file for each sample which is located in a
# # directory with all sample statistics.
# intermediate_date_time=$(date +"%F %H:%M:%S")
# echo "${intermediate_date_time}"
# echo "Calculating statistics on contigs"
# mkdir -p "${seqkit_output_dir}"/"${seqkit_stats_date_time}"_megahit_all_tsv
# seqkit_all_tsv_dir="${seqkit_output_dir}"/"${seqkit_stats_date_time}"_megahit_all_tsv
# for FILE_DIR in "${megahit_output_dir}"/"${megahit_date_time}"_*.megahit_asm 
# do 
#  intermediate_date_time=$(date +"%F %H:%M:%S")
#  echo "${intermediate_date_time}"
#  base_name=$(basename "$FILE_DIR" )
#  SAMPLE=$(echo "${base_name}" | sed "s/\.megahit_asm//" |sed "s/${megahit_date_time}_//")
#  contigs_fasta_file="${FILE_DIR}"/"${megahit_date_time}"_"${SAMPLE}"_final.contigs.fa
#  echo "Calculating statistics on ${SAMPLE} contigs"
#  seqkit fx2tab -l -g -n -i -H -B AT "${contigs_fasta_file}" | \
#     csvtk -t -C '&'  rename -f "#id",length -n contig_id,contig_length | \
#      awk -v sample="${SAMPLE}" 'BEGIN {FS=OFS="\t"}
#       FNR==1 {print "sample",$0};
#       FNR>1 {print sample,$0}'  > \
#       "${seqkit_all_tsv_dir}"/"${seqkit_stats_date_time}"_megahit_"${SAMPLE}"_seqkit.tsv;
# done
# intermediate_date_time=$(date +"%F %H:%M:%S")
# echo "${intermediate_date_time}"





# # Calculate statistics for each contig in each dereplicated bin: contig IDs, contig lengths, GC %, 
# # AT %. Output is a tab-separated file for each bin which is located in a
# # directory with all sample statistics.
# intermediate_date_time=$(date +"%F %H:%M:%S")
# echo "${intermediate_date_time}"
# echo "Calculating statistics on dereplicated bins"
# mkdir -p "${seqkit_output_dir}"/"${seqkit_stats_date_time}"_sa_95perc_all_tsv
# seqkit_all_tsv_dir="${seqkit_output_dir}"/"${seqkit_stats_date_time}"_sa_95perc_all_tsv

# for BIN_FASTA_FILE in "${drep_output_dir}"/"${drep_date_time}"_sa_95perc/dereplicated_genomes/"${metabat2_date_time}"_*_bin.*.fa
# do
#  base_name=$(basename "${BIN_FASTA_FILE}" )
#  base_name=$(echo "${base_name}" | sed "s/${metabat2_date_time}_//"| sed "s/\.fa//")
#  SAMPLE=$(echo "${base_name}" | sed "s/_bin\.[0-9]*//")
#  BIN_ID=$(echo "${base_name}" | sed "s/${SAMPLE}_bin\.//")
#  seqkit fx2tab -l -g -n -i -H -B AT "${BIN_FASTA_FILE}" | \
#     csvtk -t -C '&'  rename -f "#id",length -n contig_id,contig_length | \
#      awk -v bin_id="${BIN_ID}" -v sample="${SAMPLE}" 'BEGIN {FS=OFS="\t"}
#       FNR==1 {print "sample","bin_id",$0};
#       FNR>1 {print sample,bin_id,$0}'  > \
#       "${seqkit_all_tsv_dir}"/"${seqkit_stats_date_time}"_sa_95perc_"${SAMPLE}"_bin."${BIN_ID}"_seqkit.tsv;
# done
# echo "Done"

# # Also for metabat2 bins
# mkdir -p "${seqkit_output_dir}"/"${seqkit_stats_date_time}"_metabat2_bins_all_tsv
# seqkit_all_tsv_dir="${seqkit_output_dir}"/"${seqkit_stats_date_time}"_metabat2_bins_all_tsv

# for BIN_FASTA_FILE in "${metabat2_output_dir}"/"${metabat2_date_time}"*bins/"${metabat2_date_time}"_*bin.*.fa
# do 
#  base_name=$(basename "${BIN_FASTA_FILE}" )
#  base_name=$(echo "${base_name}" | sed "s/${metabat2_date_time}_//"| sed "s/\.fa//")
#  SAMPLE=$(echo "${base_name}" | sed "s/_bin\.[0-9]*//")
#  BIN_ID=$(echo "${base_name}" | sed "s/${SAMPLE}_bin\.//")
#  seqkit fx2tab -l -g -n -i -H -B AT "${BIN_FASTA_FILE}" | \
#     csvtk -t -C '&'  rename -f "#id",length -n contig_id,contig_length | \
#      awk -v bin_id="${BIN_ID}" -v sample="${SAMPLE}" 'BEGIN {FS=OFS="\t"}
#       FNR==1 {print "sample","bin_id",$0};
#       FNR>1 {print sample,bin_id,$0}'  > \
#       "${seqkit_all_tsv_dir}"/"${seqkit_stats_date_time}"_metabat2_all_"${SAMPLE}"_bin."${BIN_ID}"_seqkit.tsv;
# done
# echo "Done"



# [NOT RUN] This calculates stats for bins within each sample and creates a directory for each sample
# for FILE_DIR in "${metabat2_output_dir}"/"${metabat2_date_time}"*bins/"${metabat2_date_time}"_*bin.*.fa
# do 
#  SAMPLE=$(echo "${FILE_DIR}" | sed "s/${metabat2_date_time}_//" | sed "s/_bins//")
#  base_name=$(basename "$SAMPLE" )
#  mkdir -p "${seqkit_output_dir}"/"${metabat2_date_time}"_"${base_name}"_seqkit_stats
#  seqkit_stats_out_dir="${seqkit_output_dir}"/"${metabat2_date_time}"_"${base_name}"_seqkit_stats
#  echo "Calculating statistics on ${base_name}" 
#   for BIN_FASTA_FILE in "${metabat2_output_dir}"/"${metabat2_date_time}"_"${base_name}"_bins/*.fa
#   do
#     BIN_ID=$(basename "${BIN_FASTA_FILE}" )
#     BIN_ID=$(echo "${BIN_ID}" | sed "s/${metabat2_date_time}_//"| sed "s/\.fa//")
#     seqkit fx2tab -l -g -n -i -H -B AT "${BIN_FASTA_FILE}" | \
#      csvtk -t -C '&'  rename -f "#id",length -n contig_id,contig_length \
#      > "${seqkit_stats_out_dir}"/"${BIN_ID}"_seqkit.tsv;
#   done
#  echo "Done";
# done

# -l, --length    print sequence length
# -g, --gc    print GC content (%)
# -n, --name    only print names (no sequences and qualities)
# -i, --only-id    print ID instead of full head
# -H, --header-line    print header line (#id     length  GC)
# -B, --base-content strings    print base content (%). (case ignored, multiple values supported) e.g. -B AT -B N

# csvtk:
#   1. By default, csvtk assumes input files have header row, if not, switch flag "-H" on.
#   2. By default, csvtk handles CSV files, use flag "-t" for tab-delimited files.

# rename:
#   -f, --fields string   select only these fields. e.g -f 1,2 or -f columnA,columnB
#   -F, --fuzzy-fields    using fuzzy fields, e.g., -F -f "*name" or -F -f "id123*"
#   -n, --names string    comma separated new names



# Merge tables
# Join per-bin stats into a per-sample file: all bins of each sample into one table. 
# Add a column with sample name and bin number.
# Output is tab-separated files with stats for all bins in each sample (sample name,
# bin number, contig ID, contig length, GC content, AT content). Each sample has
# its own tab-separated file.

# If the row number is 1 (header), print the header: Sample, Bin, everything else.
# If the row number is >1, we split the file name into parts that are separated by
# slash (path to the file). Then, the variable 'filename' will contain only the last 
# part of the 'parts' object, which is the name of the file without parent directories.
# We remove the '_seqkit.tsv' with sub, and split the filename by underscore.
# The sample name is the first item ('sample_id'), the bin number is the second item ('bin_id').
# In the output file, the first column is sample_id, the second column is bin_id, and then
# we print the remaining contents (contig_id, length, GC, AT, etc).
# seqkit_all_tsv_dir="${seqkit_output_dir}"/"${seqkit_stats_date_time}"_sa_95perc_all_tsv
# # Merge files without the first line then add header with csvtk
# # https://stackoverflow.com/questions/24491801/concatenate-multiple-files-skipping-first-line
# awk FNR!=1 "${seqkit_all_tsv_dir}"/"${seqkit_stats_date_time}"_sa_95perc_*_seqkit.tsv | \
#  csvtk -t add-header  -n sample,bin_id,contig_id,contig_length,GC,AT > \
#   "${seqkit_output_dir}"/"${seqkit_stats_date_time}"_sa_95perc_combined_stats.tsv

# # Do the same for MetaBAT2 stats (all bins)
# seqkit_all_tsv_dir="${seqkit_output_dir}"/"${seqkit_stats_date_time}"_metabat2_bins_all_tsv
# awk FNR!=1 "${seqkit_all_tsv_dir}"/"${seqkit_stats_date_time}"_metabat2_all_*_seqkit.tsv | \
#  csvtk -t add-header  -n sample,bin_id,contig_id,contig_length,GC,AT > \
#   "${seqkit_output_dir}"/"${seqkit_stats_date_time}"_metabat2_combined_stats.tsv

# # Do the same for MEGAHIT stats (be careful of the different columns: no bin_id)
# seqkit_all_tsv_dir="${seqkit_output_dir}"/"${seqkit_stats_date_time}"_megahit_all_tsv
# awk FNR!=1 "${seqkit_all_tsv_dir}"/"${seqkit_stats_date_time}"_megahit_*_seqkit.tsv | \
#  csvtk -t add-header  -n sample,contig_id,contig_length,GC,AT > \
#   "${seqkit_output_dir}"/"${seqkit_stats_date_time}"_megahit_combined_stats.tsv

# Sanity check
# wc_count=0
# for FILE in "${seqkit_output_dir}"/"${seqkit_stats_date_time}"_sa_95perc_all_tsv/"${seqkit_stats_date_time}"_sa_95perc_*_seqkit.tsv
# do 
#  count=$(tail -n +2 $FILE  |wc -l )
#  wc_count=$((wc_count + count))
# done
# echo $wc_count
# tail -n +2 "${seqkit_output_dir}"/"${seqkit_stats_date_time}"_sa_95perc_combined_stats.tsv |wc -l

# wc_count=0
# for FILE in "${seqkit_output_dir}"/"${seqkit_stats_date_time}"_metabat2_bins_all_tsv/"${seqkit_stats_date_time}"_metabat2_all_*_seqkit.tsv
# do 
#  count=$(tail -n +2 $FILE  |wc -l )
#  wc_count=$((wc_count + count))
# done
# echo $wc_count
# tail -n +2 "${seqkit_output_dir}"/"${seqkit_stats_date_time}"_metabat2_combined_stats.tsv |wc -l

# wc_count=0
# for FILE in "${seqkit_output_dir}"/"${seqkit_stats_date_time}"_megahit_all_tsv/"${seqkit_stats_date_time}"_megahit_*_seqkit.tsv
# do 
#  count=$(tail -n +2 $FILE  |wc -l )
#  wc_count=$((wc_count + count))
# done
# echo $wc_count
# tail -n +2 "${seqkit_output_dir}"/"${seqkit_stats_date_time}"_megahit_combined_stats.tsv |wc -l




# #  seqkit_stats_out_dir="${seqkit_output_dir}"/"${seqkit_stats_date_time}"_"${base_name}"_seqkit_stats
# #  stats_files=("${seqkit_stats_out_dir}"/*.tsv) # Store matching files in an array
# #  if [[ ${#stats_files[@]} -eq 0 ]]; then
# #     echo "Error: No *.tsv files found in the directory ${seqkit_stats_out_dir}"
# #     continue
# #  else
# #     echo "Combining tables ${base_name}"
# #     awk 'BEGIN {FS=OFS="\t"} 
# #         FNR==1 && NR==1 {print "Sample", "Bin", $0} 
# #         FNR>1 {split(FILENAME, parts, "/"); filename=parts[length(parts)]; \
# #             sub(/_seqkit.tsv$/, "", filename); 
# #             split(filename, parts, "_"); sample_id=parts[1]; bin_id=parts[2];
# #             sub(/bin./, "", bin_id); \
# #             print sample_id, bin_id, $0}' \
# #             "${seqkit_stats_out_dir}"/*.tsv > \
# #         "${seqkit_output_dir}"/"${seqkit_stats_date_time}"_"${base_name}"_combined_stats.tsv
# #  fi
# #  echo "Done";
# # done

# # Join all sample stats into one tab-separated file. 
# # Output is a tab-separated file with stats for all bins in each sample (sample name,
# # bin number, contig ID, contig length, GC content, AT content). All samples are 
# # gathered in one file. To avoid duplications in contig IDs, I changed the ID column 
# # to a concatenation of Sample, Bin, and contig ID.

# # If the row number is 1 (header), we declare that 
# # the first three columns of the header are separated by tab (Sample, Bin, contig_id), 
# # then print these columns. Then, we print all the remaining columns from the bin stats file 
# # (contig_length, GC%, AT%, etc), which are also separated by tab 
# # (for (i = 4; i <= NF; i++) printf "\t%s", $i;).
# # We use a for loop because we don't know how many columns are remaining.
# # Then, we print a new line.

# # If the row number is >1, we declare that 
# # the first three columns are separated by tab ("%s\t%s\t%s"), then we 
# # print their contents: column 1 of the bin stats file (Sample), column 2 (bin),
# # concatenation of columns 1, 2, and 3 (Sample, bin, contig_id), and the 
# # remaining columns that will be separated by tab. The new column 3 is the new contig_id.
# # Then, print the new line and move to the next row.
# # awk 'BEGIN {FS=OFS="\t"} 
# #     FNR==1 && NR==1 {printf  "%s\t%s\t%s", "Sample", "Bin", "contig_id"; 
# #       for (i = 4; i <= NF; i++) printf "\t%s", $i;
# #       print "";} 
# #     FNR>1 {printf "%s\t%s\t%s", $1, $2, $1"_"$2"_"$3; for (i = 4; i <= NF; i++) printf "\t%s", $i;
# #       print "";}' \
# #         "${seqkit_output_dir}"/*_combined_stats.tsv > \
# #     "${seqkit_output_dir}"/"${seqkit_stats_date_time}"_all_samples_stats.tsv

# # Sanity check
# # wc_count=0
# # for FILE in "${seqkit_output_dir}"/*_combined_stats.tsv
# # do 
# #  count=$(grep -v  "Sample" $FILE  |wc -l )
# #  wc_count=$((wc_count + count))
# # done
# # echo $wc_count
# # wc "${seqkit_output_dir}"/"${seqkit_stats_date_time}"_all_samples_stats.tsv


# wc_count=0
# for FILE in "${seqkit_all_tsv_dir}"/"${seqkit_stats_date_time}"_sa_95perc_*_seqkit.tsv
# do 
#  count=$(tail -n +2 $FILE  |wc -l )
#  wc_count=$((wc_count + count))
# done
# echo $wc_count
# wc -l "${seqkit_output_dir}"/"${seqkit_stats_date_time}"_all_samples_stats.tsv

# # Calculate average GC% and AT% per bin and select bins with high average AT%
# # head -n 390 20250303_17_54_26_all_samples_stats.tsv >stats_head.tsv

# # Sample is column 1, bin is column 2, GC% is column 5, AT% is column 6
# # The gsub(/"/, "") function used to remove all double quotes (") from each 
# # line before processing the data because TSV and CSV files often enclose 
# # their columns with double quotes. Removing double quotes ensures correct
# # processing.
# # /"/: The regular expression that matches double quotes (").

# # I create sample_bin_GC and sample_bin_AT, which are associative arrays 
# # (also called hash tables or dictionaries in other programming languages) in awk.
# # In awk, arrays can have string keys, not just numeric indices.
# # The key in both sample_bin_GC and sample_bin_AT is key = $1 "_" $2 
# # (which combines Sample and bin into a single string).
# # The values stored in these arrays (count) are the cumulative sums
# #  of GC% ($5) and AT% ($6), respectively.
# # In the end, we print the key (e.g. SampleA_bin1), the mean GC%, and the mean AT%
# awk 'BEGIN {FS=OFS="\t"} 
#     FNR==1 && NR==1 {printf  "%s\t%s\t%s", "Sample_Bin", "MeanGC", "MeanAT"; 
#       print "";} 

#       FNR >1 {
#         gsub(/"/, "")
#         key = $1 "_" $2
#         sample_bin_GC[key] += $5
#         sample_bin_AT[key] += $6
#         ++count[key]
#       }
#       END {
#         for (key in sample_bin_GC)
#           print key, sample_bin_GC[key]/count[key], sample_bin_AT[key]/count[key]
#       }' "${seqkit_output_dir}"/"${seqkit_stats_date_time}"_all_samples_stats.tsv > \
#         "${seqkit_output_dir}"/"${seqkit_stats_date_time}"_average_gc_at_bins.tsv

# # Filter to keep rows with low GC (below 30) (column 2)
# awk -v threshold=30 'NR==1 || $2 < threshold'  "${seqkit_output_dir}"/"${seqkit_stats_date_time}"_average_gc_at_bins.tsv \
#         > "${seqkit_output_dir}"/"${seqkit_stats_date_time}"_low_gc_bins.tsv


# # Extract certain contig by its id
# # seqkit grep -n output/mag_assembly/metabat2_output/20250211_22_55_27_2D14/20250211_22_55_27_2D14_bin.13.fa \
# #  -p "k141_109001" |grep -v "^>"| tr -d '\n' > output/mag_assembly/2D14_13_k141_109001.fasta

# # seqkit grep -n output/mag_assembly/metabat2_output/20250211_22_55_27_2D14/20250211_22_55_27_2D14_bin.10.fa \
# #  -p "k141_237915" | seqkit subseq -r 1:10000 > output/mag_assembly/metabat2_output/2D14_bin.10.k141_237915.fasta

# # seqkit grep -n output/mag_assembly/megahit_output/20250211_22_55_27_2D14.megahit_asm/20250211_22_55_27_2D14_final.contigs.fa \
# #  -p "k141_0 flag=1 multi=1.0000 len=364"  > output/mag_assembly/metabat2_output/2D14_k141_0.fa

# How different are read counts in our samples?
seqkit stats data/bowtie2_decontam_fastq/*_trim_decontam_R*.fastq.gz \
  --tabular \
  --threads 4 > \
  output/mag_assembly/seqkit_output/20250728_14_27_54_trim_decontam_read_count.tsv 
