#!/usr/bin/env bash
#SBATCH -t 360:00:00
#SBATCH -N 4
#SBATCH -n 40
#SBATCH --mem-per-cpu 8g
#SBATCH -J 20250515_11-06-mag_assembly
#SBATCH --output jobreports/20250515_11-06-mag_assembly-megahit-metabat2-pipeline-out.txt
#SBATCH --error jobreports/20250515_11-06-mag_assembly-megahit-metabat2-pipeline-out.txt
#I am requesting 4 nodes containing 40 CPUs, with 8 GB memory per CPU. Total: 320 GB. 280 GB for GTDB-TK
source ~/miniconda3/etc/profile.d/conda.sh
shopt -s nullglob
# When nullglob is enabled, if a glob pattern does not match any files,
# it expands to nothing (an empty string) instead of returning the pattern itself.
# So, if no matches are found, the script will skip the file

# MEGAHIT assembles reads into contigs
# Metaquast is a metagenome assembly evaluation tool
# MetaBAT2 is a binning software: group contigs into bins (genomes)
# Metaquast is also used to evaluate the quality of bins
date_var=$(date -I|sed 's/-//g')
time_var=$(date +%T |sed 's/:/_/g' )
date_time=${date_var}_${time_var}
start_date_time=$(date +"%F %H:%M:%S")
megahit_date_time=20250303_17_54_26
metaquast_megahit_date_time=20250303_17_54_26
bbwrap_date_time=20250326_08_18_26
nthreads=40
nthreads_sort=35
mem_req=8G
mem_req_sort=4G
project_home_dir=~/projects/metagenome
output_dir=~/projects/metagenome/output/mag_assembly
bowtie2_decontam_fastq_dir=data/bowtie2_decontam_fastq
megahit_output_dir=output/mag_assembly/megahit_output
megahit_aligned_reads_dir=output/mag_assembly/megahit_output/alignedreads
contig_coverages_dir=output/mag_assembly/contig_coverages
metaquast_output_dir=output/mag_assembly/metaquast_output
metaquast_script_dir=~/quast-5.2.0
samtools_reports_dir=output/mag_assembly/samtools_reports
bbwrap_refs_dir=output/mag_assembly/bbwrap_refs

cd ${project_home_dir}
mkdir -p output/mag_assembly/megahit_output/alignedreads
mkdir -p output/mag_assembly/samtools_reports
mkdir -p output/mag_assembly/bbwrap_refs
mkdir -p output/mag_assembly/metaquast_output
mkdir -p output/mag_assembly/contig_coverages

# 0. Show the current time for logging
echo "${start_date_time}"

# 1. Activate conda environment
conda activate mag_assembly-tools

# megahit --test  # run on a toy dataset
# cd ../..
# 2. Assemble reads into contigs with MegaHIT
# echo "Starting MegaHIT"
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
# for FILE in "${bowtie2_decontam_fastq_dir}"/*decontam_R1.fastq.gz
# do 
#  SAMPLE=$(echo "${FILE}" | sed "s/_trim_decontam_R1\.fastq\.gz//")
#  base_name=$(basename "$SAMPLE" )
#  echo "Running MegaHIT on ${base_name}"
#  megahit \
#     -1 "${bowtie2_decontam_fastq_dir}"/"${base_name}"_trim_decontam_R1.fastq.gz \
#     -2 "${bowtie2_decontam_fastq_dir}"/"${base_name}"_trim_decontam_R2.fastq.gz \
#     -o "${megahit_output_dir}"/"${date_time}"_"${base_name}".megahit_asm \
#     -t "${nthreads}" 2>&1 |tee "${megahit_output_dir}"/"${date_time}"_"${base_name}"_megahit_report.txt
#  mv "${megahit_output_dir}"/"${date_time}"_"${base_name}".megahit_asm/final.contigs.fa \
#   "${megahit_output_dir}"/"${date_time}"_"${base_name}".megahit_asm/"${date_time}"_"${base_name}"_final.contigs.fa;
# done

# This one actually has one contig
# megahit \
#     -1 ~/miniconda3/envs/mag_assembly-tools/share/megahit/test_data/r3_1.fa \
#     -2 ~/miniconda3/envs/mag_assembly-tools/share/megahit/test_data/r3_2.fa \
#     -o "${megahit_output_dir}"/"${date_time}"_r3.megahit_asm \
#     -t "${nthreads}" 2>&1 |tee "${megahit_output_dir}"/"${date_time}"_r3_megahit_report.txt
# generates a directory 2D10_subset.megahit_asm 
# which contains checkpoints.txt, done, final.contigs.fa, intermediate_contig/, log, options.json

# --k-min : minimum kmer size (<= 255), must be odd number [default 21] 
# --k-max : maximum kmer size (<= 255), must be odd number [default 141]
# --k-step : increment of kmer size of each iteration (<= 28), must be even number [default 12]

# Tip 1. for ultra complex metagenomics data such as soil, a larger kmin, say 27,
#  is recommended to reduce the complexity of the de Bruijn graph. Quality trimming is also recommended
# Tip 2. for high-depth generic data, large --k-min (25 to 31) is recommended
# Tip 3. smaller --k-step, say 10, is more friendly to low-coverage datasets

# 3. Contig QC with QUAST
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
# for FILE_DIR in "${megahit_output_dir}"/"${megahit_date_time}"*megahit_asm
# do 
#  SAMPLE=$(echo "${FILE_DIR}" | sed "s/\.megahit_asm//" |sed "s/${megahit_date_time}_//")
#  base_name=$(basename "$SAMPLE" )
#  echo "Running MetaQuast on ${base_name}" 
#  "${metaquast_script_dir}"/metaquast.py \
#     "${megahit_output_dir}"/"${megahit_date_time}"_"${base_name}".megahit_asm/"${megahit_date_time}"_"${base_name}"_final.contigs.fa  \
#     -o "${metaquast_output_dir}"/"${metaquast_megahit_date_time}"_metaquast_megahit_"${base_name}" \
#     --threads "${nthreads}" \
#     2>&1 |tee "${metaquast_output_dir}"/"${metaquast_megahit_date_time}"_"${base_name}"_metaquast_megahit_report.txt;
# done
# General usage:
#     python metaquast.py contigs_1 contigs_2 ... -r reference_1,reference_2,reference_3,...
# One can provide several files or directories with multiple reference files inside 
# with -r option. 
# Option -r may be specified multiple times or all references may be specified as
#  a comma-separated list (without spaces!) with a single -r option beforehand

# 4. Extract mapped and unmapped reads
# 4.0 Use BBmap and samtools
conda deactivate
conda activate qc-tools

# 4.1 Align reads with bbwrap.sh:
# Then, convert file _mapped_and_unmapped.sam to _mapped_and_unmapped.bam (samtools view -b)
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
# for FILE_DIR in "${megahit_output_dir}"/"${megahit_date_time}"*megahit_asm
# do 
#  SAMPLE=$(echo "${FILE_DIR}" | sed "s/\.megahit_asm//" |sed "s/${megahit_date_time}_//")
#  base_name=$(basename "$SAMPLE" )
#  echo "Running bbwrap.sh on ${base_name}"
#  mkdir -p "${bbwrap_refs_dir}"/"${date_time}"_"${base_name}"_refs
#  bbwrap.sh ref="${FILE_DIR}"/"${megahit_date_time}"_"${base_name}"_final.contigs.fa \
#     in="${bowtie2_decontam_fastq_dir}"/"${base_name}"_trim_decontam_R1.fastq.gz \
#     in2="${bowtie2_decontam_fastq_dir}"/"${base_name}"_trim_decontam_R2.fastq.gz \
#     out="${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_mapped_and_unmapped.sam \
#     kfilter=22 \
#     subfilter=15 \
#     path="${bbwrap_refs_dir}"/"${bbwrap_date_time}"_"${base_name}"_refs \
#     maxindel=80 2>&1 |tee "${megahit_output_dir}"/"${bbwrap_date_time}"_"${base_name}"_bbwrap_report.txt
#  echo "Running samtools view -b on ${base_name}"
#  samtools view -b "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_mapped_and_unmapped.sam.gz > \
#     "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_mapped_and_unmapped.bam \
#     2>&1 |tee "${samtools_reports_dir}"/"${bbwrap_date_time}"_"${base_name}"_samtools_report.txt;
# done

## TEST commands: work
# mv "${megahit_output_dir}"/"${date_time}"_r3.megahit_asm/final.contigs.fa \
#   "${megahit_output_dir}"/"${date_time}"_r3.megahit_asm/"${date_time}"_r3_final.contigs.fa
# bbwrap.sh ref="${megahit_output_dir}"/"${date_time}"_r3.megahit_asm/"${date_time}"_r3_final.contigs.fa \
#     in=~/miniconda3/envs/mag_assembly-tools/share/megahit/test_data/r3_1.fa \
#     in2=~/miniconda3/envs/mag_assembly-tools/share/megahit/test_data/r3_2.fa \
#     out="${megahit_aligned_reads_dir}"/"${date_time}"_r3_mapped_and_unmapped.sam.gz \
#     kfilter=22 \
#     subfilter=15 \
#     maxindel=80 2>&1 |tee "${megahit_output_dir}"/"${date_time}"_r3_bbwrap_report.txt
###


# 4.2 Output per contig coverage to cov.txt with pileup.sh:
# BBMap generates coverage information by internally using Pileup
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
for FILE in "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"*_mapped_and_unmapped.bam
do 
 SAMPLE=$(echo "${FILE}" | sed "s/_mapped_and_unmapped\.bam//" |sed "s/${bbwrap_date_time}_//")
 base_name=$(basename "$SAMPLE" )
 echo "Running pileup.sh on ${base_name}"
 pileup.sh in="${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_mapped_and_unmapped.bam \
    out="${contig_coverages_dir}"/"${bbwrap_date_time}"_"${base_name}"_cov.txt \
    2>&1 |tee "${contig_coverages_dir}"/"${bbwrap_date_time}"_"${base_name}"_pileup_report.txt;
done


# 4.3 Extract unmapped reads (SE to unmapped.se.fq and PE to unmapped.pe.fq) (samtools view -u -f4)
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
# for FILE in "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"*_mapped_and_unmapped.bam
# do 
#  SAMPLE=$(echo "${FILE}" | sed "s/_mapped_and_unmapped\.bam\//" |sed "s/${bbwrap_date_time}_//")
#  base_name=$(basename "$SAMPLE" )
#  echo "Running samtools bam2fq on ${base_name}"
#  samtools view -u -f4 "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_mapped_and_unmapped.bam | \
#     samtools bam2fq -s "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_unmapped.se.fq - > \
#     "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_unmapped.pe.fq \
#     2>&1 |tee "${samtools_reports_dir}"/"${bbwrap_date_time}"_"${base_name}"_samtools_report.txt
#  gzip -9 --best "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_unmapped.se.fq
#  gzip -9 --best "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_unmapped.pe.fq;
# done

# samtools view – views and converts SAM/BAM/CRAM files
# -u, --uncompressed : Output uncompressed data
# -f4 flag filters out reads that are unmapped.
# -b, --bam :Output in the BAM format.

# bam2fq converts a bam into FASTQ format
# -s FILE : Write singleton reads in FASTQ format to FILE instead of outputting them.

# The final output is then directed to the standard output (stdout) using the - symbol 
#   followed by the redirection operator >
# The - > combination redirects the output of the second command (samtools bam2fq)
#    to a file named output/unmapped.pe.fq. The hyphen (-) after the redirection operator
#    indicates that the output should be sent to stdout (standard output), 
#    which can then be captured or saved to a file as needed.

# tee -a appends instead of overwriting 
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
## 4.4 Keep mapped reads then sort
# To get only the mapped reads use the parameter F, which works like -v of grep and skips the 
# alignments for a specific flag (samtools view -b -F 4).
# Then, sort _either_read_mapped.bam file with mapped reads by coordinate (no options used)
# (samtools sort -m)
# for FILE in "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"*_mapped_and_unmapped.bam
# do 
#  SAMPLE=$(echo "${FILE}" | sed "s/_mapped_and_unmapped\.bam\.gz//" |sed "s/${bbwrap_date_time}_//")
#  base_name=$(basename "$SAMPLE" )
#  echo "Running samtools view -b -F 4 on ${base_name}"
#  samtools view -b -F 4 \
#    "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_mapped_and_unmapped.bam \
#    > "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_either_read_mapped.bam \
#    2>&1 |tee -a "${samtools_reports_dir}"/"${bbwrap_date_time}"_"${base_name}"_samtools_report.txt
#  echo "Running samtools sort on ${base_name}"
#  samtools sort -m "${mem_req_sort}" -@ "${nthreads_sort}" \
#    "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_either_read_mapped.bam \
#    -o "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_either_read_mapped_sorted.bam \
#    2>&1 |tee -a "${samtools_reports_dir}"/"${bbwrap_date_time}"_"${base_name}"_samtools_report.txt;
# done

# This line is only if you sort selected samples
# for FILE in "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_2D14_either_read_mapped.bam.gz \
#  "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_G14_either_read_mapped.bam.gz \
#  "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_H15_either_read_mapped.bam.gz 
# do 
#  SAMPLE=$(echo "${FILE}" | sed "s/_either_read_mapped\.bam\.gz//" |sed "s/${bbwrap_date_time}_//")
#  base_name=$(basename "$SAMPLE" )
#  gunzip "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_either_read_mapped.bam.gz
#  echo "Running samtools sort on ${base_name}"
#  samtools sort -m "${mem_req_sort}" -@ "${nthreads_sort}" \
#    "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_either_read_mapped.bam \
#    -o "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_either_read_mapped_sorted.bam \
#    2>&1 |tee -a "${samtools_reports_dir}"/"${bbwrap_date_time}"_"${base_name}"_samtools_report.txt
#  gzip -9 --best "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_either_read_mapped.bam;
# done
# -b output in BAM format
# -F 4: Exclude unmapped reads. Keep only mapped reads (evem if mate is unmapped)

## 4.5 Filter unmapped read pairs (unmapped to contigs) (samtools view -b -f 12 -F 256)
### SAMtools SAM-flag filter: get unmapped pairs (both reads R1 and R2 unmapped).
# Then, sort _bothReadsUnmapped.bam file by read name ( -n ) to have paired reads next to each other 
### ("${nthreads}" parallel threads, each using up to 5G memory) (samtools sort)
# Finally, split paired-end reads into separated fastq files .._R1 .._R2 (samtools fastq)
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
# for FILE in "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"*mapped_and_unmapped.bam 
# do 
#  SAMPLE=$(echo "${FILE}" | sed "s/_mapped_and_unmapped\.bam//" |sed "s/${bbwrap_date_time}_//")
#  base_name=$(basename "$SAMPLE" )
#  echo "Running samtools view -b -f 12 -F 256 on ${base_name}" 
#  samtools view -b -f 12 -F 256 \
#    "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_mapped_and_unmapped.bam \
#    > "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_bothReadsUnmapped.bam \
#    2>&1 |tee -a "${samtools_reports_dir}"/"${bbwrap_date_time}"_"${base_name}"_samtools_report.txt
#  echo "Running samtools sort -n and fastq on ${base_name}"
#  samtools sort -n -m "${mem_req_sort}" -@ "${nthreads_sort}" \
#    "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_bothReadsUnmapped.bam \
#    -o "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_bothReadsUnmapped_sorted.bam \
#    2>&1 |tee -a "${samtools_reports_dir}"/"${bbwrap_date_time}"_"${base_name}"_samtools_report.txt;
#  samtools fastq -@ "${nthreads_sort}" \
#    "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_bothReadsUnmapped_sorted.bam \
#    -1 "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_unmapped_R1.fastq.gz \
#   	-2 "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_unmapped_R2.fastq.gz \
#   	-0 /dev/null -s /dev/null -n \
#    2>&1 |tee -a "${samtools_reports_dir}"/"${bbwrap_date_time}"_"${base_name}"_samtools_report.txt;
# done
###-f  12    # Extract only ( -f ) alignments with both reads unmapped: <read unmapped><mate unmapped>
###-F 256    # Do not (  -F  ) extract alignments which are: <not primary alignment>

# samtools sort: Sort alignments by leftmost coordinates, 
#  by read name when -n or -N are used,
#  by tag contents with -t,
#  or a minimiser-based collation order with -M.
# -m : Approximately the maximum required memory per thread
# samtools fastq: Converts a BAM or CRAM into either FASTQ or FASTA format
#  depending on the command invoked
# [NOT RUN] 4.6 For Metabat2, we need to fix NM tags in the _either_read_mapped_sorted.bam file 
# First, index the contig FASTA file (faidx).
# Then fix NM tags in the BAM file.
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
# for FILE_DIR in "${megahit_output_dir}"/"${megahit_date_time}"*megahit_asm
# do 
#  SAMPLE=$(echo "${FILE_DIR}" | sed "s/\.megahit_asm//" |sed "s/${megahit_date_time}_//")
#  base_name=$(basename "$SAMPLE" )
#  echo "Indexing FASTA files of ${base_name}"
#  samtools faidx "${FILE_DIR}"/"${megahit_date_time}"_"${base_name}"_final.contigs.fa
#  samtools calmd -b "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_either_read_mapped_sorted.bam \
#    "${FILE_DIR}"/"${megahit_date_time}"_"${base_name}"_final.contigs.fa > \
#    "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_either_read_mapped_sorted_fixed.bam \
#    2>&1 |tee -a "${samtools_reports_dir}"/"${bbwrap_date_time}"_"${base_name}"_samtools_report.txt
#  mv "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_either_read_mapped_sorted_fixed.bam \
#    "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_either_read_mapped_sorted.bam;
# done
# 4.7 Compress bam and sam files to save space
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
# for FILE in "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"*_mapped_and_unmapped.bam
# do 
#  SAMPLE=$(echo "${FILE}" | sed "s/_mapped_and_unmapped\.bam//" |sed "s/${bbwrap_date_time}_//")
#  base_name=$(basename "$SAMPLE" )
#  echo "gzipping ${base_name}"
#  gzip -9 --best "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_either_read_mapped.bam
#  gzip -9 --best "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_bothReadsUnmapped.bam
#  gzip -9 --best "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_bothReadsUnmapped_sorted.bam
#  gzip -9 --best "${megahit_aligned_reads_dir}"/"${bbwrap_date_time}"_"${base_name}"_mapped_and_unmapped.bam;
# done
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
