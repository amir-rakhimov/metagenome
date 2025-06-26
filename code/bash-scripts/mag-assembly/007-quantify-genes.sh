#!/usr/bin/env bash
#SBATCH -t 360:00:00
#SBATCH -N 4
#SBATCH -n 30
#SBATCH --mem-per-cpu 8g
#SBATCH -J 20250619_18-52-quantify-genes
#SBATCH --output jobreports/20250619_18-52-quantify-genes-out.txt
#SBATCH --error jobreports/20250619_18-52-quantify-genes-out.txt
#I am requesting 4 nodes containing 30 CPUs, with 8 GB memory per CPU. Total: 240 GB
source ~/miniconda3/etc/profile.d/conda.sh
shopt -s nullglob
# When nullglob is enabled, if a glob pattern does not match any files,
# it expands to nothing (an empty string) instead of returning the pattern itself.
# So, if no matches are found, the script will skip the file
date_var=$(date -I|sed 's/-//g')
time_var=$(date +%T |sed 's/:/_/g' )
date_time="${date_var}_${time_var}"
start_date_time=$(date +"%F %H:%M:%S")
megahit_date_time=20250303_17_54_26
bbwrap_date_time=20250326_08_18_26
picard_date_time=20250603_18_19_44
htseq_count_date_time=20250619_23_34_58
htseq_count_date_time="${date_time}"
# prokka_date_time=20250515_11_29_31
prokka_date_time=20250522_10_37_23
#
java_path=~/jdk-24.0.1/bin/java
picard_jar_path=~/picard.jar
project_home_dir=~/projects/metagenome
output_dir=~/projects/metagenome/output/mag_assembly
megahit_aligned_reads_dir=output/mag_assembly/megahit_output/alignedreads
picard_output_dir=output/mag_assembly/picard_output
prokka_output_dir=output/mag_assembly/prokka_output
prokka_gtf_dir=output/mag_assembly/prokka_gtf_files
gene_count_output_dir=output/mag_assembly/gene_count_output
gene_lengths_dir=output/mag_assembly/gene_lengths
tpm_files_dir=output/mag_assembly/tpm_files
# index_files_dir=output/mag_assembly/index_files
#
mkdir -p output/mag_assembly/picard_output
# mkdir -p output/mag_assembly/index_files
mkdir -p output/mag_assembly/prokka_gtf_files
mkdir -p output/mag_assembly/gene_count_output
mkdir -p output/mag_assembly/gene_lengths
mkdir -p output/mag_assembly/tpm_files

echo "${start_date_time}"
conda activate qc-tools

# Removing duplicates in the sorted and indexed BAM file
# https://www.biostars.org/p/15818/
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
# for FILE in "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"*_either_read_mapped_sorted.bam.gz
# do 
#  SAMPLE=$(echo "${FILE}" | sed "s/_either_read_mapped_sorted\.bam\.gz//" |sed "s/${bbwrap_date_time}_//")
#  base_name=$(basename "$SAMPLE" )
#  gunzip "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_either_read_mapped_sorted.bam.gz
#  # Add Read group because Picard can't run without it
#  echo "Adding read group to ${base_name}"
#  samtools addreplacerg \
#     -r 'ID:group1' \
#     -r 'LB:lib1' \
#     -r 'PL:ILLUMINA' \
#     -r 'PU:unit1' \
#     -r 'SM:sample1' \
#     -o "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_either_read_mapped_sorted_RG.bam \
#     "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_either_read_mapped_sorted.bam
#  gzip -9 --best "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_either_read_mapped_sorted.bam
#  echo "Running Picard MarkDuplicates on ${base_name}"
#  "${java_path}" -Xms2g -Xmx32g -jar "${picard_jar_path}" MarkDuplicates \
#     INPUT="${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_either_read_mapped_sorted_RG.bam \
#     OUTPUT="${picard_output_dir}"/"${picard_date_time}"_"${base_name}"_map.markdup.bam \
#     METRICS_FILE="${picard_output_dir}"/"${picard_date_time}"_"${base_name}"_map.markdup.metrics \
#     AS=TRUE \
#     VALIDATION_STRINGENCY=LENIENT \
#     MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
#     REMOVE_DUPLICATES=TRUE \
#     2>&1 |tee "${picard_output_dir}"/"${picard_date_time}"_"${base_name}"_picard_markduplicates_report.txt
#  gzip -9 --best "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_either_read_mapped_sorted_RG.bam;
# done
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"


# Counting mapped reads per gene
# Convert Prokka output into GTF then count the number of reads mapped to each gene with htseq
# https://metagenomics-workshop.readthedocs.io/en/latest/annotation/quantification.html
for FILE_DIR in "${prokka_output_dir}"/"${prokka_date_time}"*_prokka
do 
 SAMPLE=$(echo "${FILE_DIR}" | sed "s/_prokka//" |sed "s/${prokka_date_time}_//")
 base_name=$(basename "$SAMPLE" )
 # Select specific data from Prokka GFF output and convert into GTF. 
 # cut -d uses ";" as a delimiter, while -f1 selects the first field
 # (everything up to ID).
 # In the final GTF file, all features are CDS because we want to count all genes.
 # Real feature types are in the 9th column ($2).
 gunzip "${picard_output_dir}"/"${picard_date_time}"_"${base_name}"_map.markdup.bam.gz
 echo "Converting Prokka output of ${base_name} into GTF"
 grep -v "#" "${FILE_DIR}"/"${prokka_date_time}"_"${base_name}"_pred_prokka.gff | \
    grep "ID=" | cut -f1 -d ';' | sed 's/ID=//g' | cut -f1,3,4,5,7,9 | \
    awk -v OFS='\t' '{print $1,"PROKKA","CDS",$3,$4,".",$5,".","gene_id", "ID="$6";Note="$2}' > \
    "${prokka_gtf_dir}"/"${prokka_date_time}"_"${base_name}"_pred_prokka.map.gtf
 # Counting the number of reads mapped to genes with htseq. Here we have to tell
 # htseq that the file is sorted by alignment coordinate -r pos.
 # The output file has two columns, the first contains the gene names and the second
 # the number of reads mapped to each gene. The last 5 lines gives you some summary 
 # information from htseq: __no_feature, __ambiguous, __too_low_aQual, __not_aligned,
 # __alignment_not_unique
 echo "Quantifying the number of reads mapped to genes with htseq on ${base_name}"
 htseq-count -r pos -t CDS -f bam "${picard_output_dir}"/"${picard_date_time}"_"${base_name}"_map.markdup.bam \
    "${prokka_gtf_dir}"/"${prokka_date_time}"_"${base_name}"_pred_prokka.map.gtf > \
    "${gene_count_output_dir}"/"${htseq_count_date_time}"_"${base_name}"_prokka.count \
    2>&1 |tee "${gene_count_output_dir}"/"${htseq_count_date_time}"_"${base_name}"_htseq_count_report.txt
 grep "ID" "${gene_count_output_dir}"/"${htseq_count_date_time}"_"${base_name}"_prokka.count | \
   sed 's/ID=//g' > \
   "${gene_count_output_dir}"/"${htseq_count_date_time}"_"${base_name}"_prokka_filtered.count
 gzip -9 --best "${picard_output_dir}"/"${picard_date_time}"_"${base_name}"_map.markdup.bam;
done
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"

# Normalizing to Transcripts Per Million (TPM)
# The gene lengths we can get from the GTF file that you used with htseq:
for FILE in "${prokka_gtf_dir}"/"${prokka_date_time}"*_pred_prokka.map.gtf
do 
 SAMPLE=$(echo "${FILE}" | sed "s/_pred_prokka\.map\.gtf//" |sed "s/${prokka_date_time}_//")
 base_name=$(basename "$SAMPLE" )
 # Here we extract only the start (#4), stop(#5) and gene name fields(#10) from the file, 
 # then remove the 'ID=' and ';Note= ...' string, print the gene name first followed by the 
 # length of the gene, change the separator to tab and store the results in the .genelengths file.
 cut -f4,5,10 "${prokka_gtf_dir}"/"${prokka_date_time}"_"${base_name}"_pred_prokka.map.gtf | \
   sed 's/ID=//g' |sed 's/\;.*//g' | gawk '{print $3,$2-$1+1}' | tr ' ' '\t' > \
   "${gene_lengths_dir}"/"${prokka_date_time}"_"${base_name}".genelengths
 # Now we can calculate TPM values using the tpm_table.py script:
 python3 code/python-scripts/tpm_table.py \
   -n "${base_name}" \
   -c "${gene_count_output_dir}"/"${htseq_count_date_time}"_"${base_name}"_prokka_filtered.count \
   -i <(echo -e ""${base_name}"\t150") \
   -l "${gene_lengths_dir}"/"${prokka_date_time}"_"${base_name}".genelengths > \
      "${tpm_files_dir}"/"${prokka_date_time}"_"${base_name}".tpm;
# We now have coverage values for all genes predicted and annotated by the PROKKA pipeline;
done

