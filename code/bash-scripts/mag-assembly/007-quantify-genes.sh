#!/usr/bin/env bash
#SBATCH -t 360:00:00
#SBATCH -N 4
#SBATCH -n 30
#SBATCH --mem-per-cpu 8g
#SBATCH -J 20250516_16-48-quantify-genes
#SBATCH --output jobreports/20250516_16-48-quantify-genes-out.txt
#SBATCH --error jobreports/20250516_16-48-quantify-genes-out.txt
#I am requesting 4 nodes containing 30 CPUs, with 8 GB memory per CPU. Total: 240 GB
source ~/miniconda3/etc/profile.d/conda.sh
shopt -s nullglob
# When nullglob is enabled, if a glob pattern does not match any files,
# it expands to nothing (an empty string) instead of returning the pattern itself.
# So, if no matches are found, the script will skip the file
date_var=$(date -I|sed 's/-//g')
time_var=$(date +%T |sed 's/:/_/g' )
date_time=${date_var}_${time_var}
start_date_time=$(date +"%F %H:%M:%S")
megahit_date_time=20250303_17_54_26
bbwrap_date_time=20250326_08_18_26
picard_date_time="${date_time}"
htseq_count_date_time="${date_time}"
prokka_date_time=20250515_11_29_31
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
# index_files_dir=output/mag_assembly/index_files
#
mkdir -p output/mag_assembly/picard_output
# mkdir -p output/mag_assembly/index_files
mkdir -p output/mag_assembly/prokka_gtf_files
mkdir -p output/mag_assembly/gene_count_output

echo "${start_date_time}"
conda activate qc-tools

# Removing duplicates in the sorted and indexed BAM file
# https://www.biostars.org/p/15818/
for FILE in "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"*_either_read_mapped_sorted.bam.gz
do 
 SAMPLE=$(echo "${FILE}" | sed "s/_either_read_mapped_sorted\.bam\.gz//" |sed "s/${bbwrap_date_time}_//")
 base_name=$(basename "$SAMPLE" )
 gunzip "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_either_read_mapped_sorted.bam.gz
 # Add Read group because Picard can't run without it
 echo "Adding read group to ${base_name}"
 samtools addreplacerg \
    -r 'ID:group1' \
    -r 'LB:lib1' \
    -r 'PL:ILLUMINA' \
    -r 'PU:unit1' \
    -r 'SM:sample1' \
    -o "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_either_read_mapped_sorted_RG.bam \
    "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_either_read_mapped_sorted.bam
 gzip "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_either_read_mapped_sorted.bam
 echo "Running Picard MarkDuplicates on ${base_name}"
 "${java_path}" -Xms2g -Xmx32g -jar "${picard_jar_path}" MarkDuplicates \
    INPUT="${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_either_read_mapped_sorted_RG.bam \
    OUTPUT="${picard_output_dir}"/"${picard_date_time}"_"${base_name}"_map.markdup.bam \
    METRICS_FILE="${picard_output_dir}"/"${picard_date_time}"_"${base_name}"_map.markdup.metrics \
    AS=TRUE \
    VALIDATION_STRINGENCY=LENIENT \
    MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
    REMOVE_DUPLICATES=TRUE \
    2>&1 |tee "${picard_output_dir}"/"${picard_date_time}"_"${base_name}"_picard_markduplicates_report.txt;
done
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"


# Counting mapped reads per gene
# Convert Prokka output into GTF then count the number of reads mapped to each gene with htseq
# https://metagenomics-workshop.readthedocs.io/en/latest/annotation/quantification.html
for FILE_DIR in "${prokka_output_dir}"/"${prokka_date_time}"*_prokka
do 
 SAMPLE=$(echo "${FILE_DIR}" | sed "s/_prokka//" |sed "s/${prokka_date_time}_//")
 base_name=$(basename "$SAMPLE" )
 # Prokka into GTF. All features are CDS because we want to count all genes. Real feature
 # types are in the 9th column ($2).
 echo "Converting Prokka output of ${base_name} into GTF"
 grep -v "#" "${FILE_DIR}"/"${prokka_date_time}"_"${base_name}"_pred_prokka.gff | \
    grep "ID=" | cut -f1 -d ';' | sed 's/ID=//g' | cut -f1,3,4,5,7,9 | \
    awk -v OFS='\t' '{print $1,"PROKKA","CDS",$3,$4,".",$5,".","gene_id", $6";"$2}' > \
    "${prokka_gtf_dir}"/"${prokka_date_time}"_"${base_name}"_pred_prokka.map.gtf
 # Counting the number of reads mapped to genes with htseq
 htseq-count -r pos -t CDS -f bam "${picard_output_dir}"/"${picard_date_time}"_"${base_name}"_map.markdup.bam \
    "${prokka_gtf_dir}"/"${prokka_date_time}"_"${base_name}"_pred_prokka.map.gtf > \
    "${gene_count_output_dir}"/"${htseq_count_date_time}"_"${base_name}"_prokka.count \
    2>&1 |tee "${gene_count_output_dir}"/"${htseq_count_date_time}"_"${base_name}"_htseq_count_report.txt
 gzip "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_map.markdup.bam;
done
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"

# Normalizing to Transcripts Per Million (TPM)
# The gene lengths we can get from the GTF file that you used with htseq:
# cut -f4,5,9,10 ${prokka_gtf_dir}/${prokka_date_time}_"${base_name}"_pred_prokka.map.gtf | \
#  sed 's/gene_id //g' | gawk '{print $4,$2-$1+1}' | tr ' ' '\t' > $SAMPLE.genelengths
# cut -f4,5,10 ${prokka_gtf_dir}/${prokka_date_time}_"${base_name}"_pred_prokka.map.gtf | \
#  sed 's/\;.*//g' | gawk '{print $4,$2-$1+1}' | tr ' ' '\t' > $SAMPLE.genelengths

Now we can calculate TPM values using the tpm_table.py script:
#TODO: read length
# code/python-scripts/tpm_table.py -n $SAMPLE -c $SAMPLE.count -i <(echo -e "$SAMPLE\t150") -l $SAMPLE.genelengths > $SAMPLE.tpm
