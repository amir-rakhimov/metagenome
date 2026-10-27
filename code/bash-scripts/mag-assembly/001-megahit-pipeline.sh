#!/usr/bin/env bash
source $HOME/miniconda3/etc/profile.d/conda.sh
# 0. Activate conda environment
conda activate mag_assembly-tools
set -euo pipefail
shopt -s nullglob
source config/bash/config.sh
# When nullglob is enabled, if a glob pattern does not match any files,
# it expands to nothing (an empty string) instead of returning the pattern itself.
# So, if no matches are found, the script will skip the file

# This script assembles reads into contigs with MEGAHIT. Then, contig quality is evaluated with Metaquast.
# Reads are aligned with bbwrap from BBMap, resulting in a BAM file {sample}_mapped_and_unmapped.bam.
# After that, unmapped reads are extracted into separate files with samtools.
# Bam file where either of the paired reads were mapped is sorted. All BAM files are compressed with gzip.
# Main final output: directories with contigs for each sample (one FASTA files with all contigs), BAM files with mapped 
# and unmapped reads, sorted BAM files where either read is mapped to contigs.
echo "$(date +"%F %H:%M:%S")"

# Function that checks if a required file exists or not
require() { [[ -s "$1" ]] || { echo "ERROR: required file missing: $1" >&2; exit 1; }; }
# Function that marks a sample as done (in a hidden directory)
is_done()   { [[ -e "${mark_dir}/${2}.${1}.done" ]]; }   # is_done STEP SAMPLE
mark_done() { touch "${mark_dir}/${2}.${1}.done"; }      # mark_done STEP SAMPLE

# 1. Read the input file list
[[ -s "${SAMPLES_TSV}" ]] || { echo "ERROR: metadata $SAMPLES_TSV missing or empty" >&2; exit 1; }
# Remove CR LF
sed -i 's/\r$//' "${SAMPLES_TSV}"

declare -A R1 R2
SAMPLES=()
while IFS=$'\t' read -r sample_id \
    || [[ -n "${sample_id}" ]]
do
    [[ "${sample_id}" == "sample_id" || -z "${sample_id}" ]] && continue
    r1="${bowtie2_decontam_fastq_dir}/${sample_id}_trim_decontam_R1.fastq.gz"
    r2="${bowtie2_decontam_fastq_dir}/${sample_id}_trim_decontam_R2.fastq.gz"
    for f in "${r1}" "${r2}"; do
        [[ -s "${f}" ]] || { echo "ERROR: missing input for ${sample_id}: ${f}" >&2; exit 1; } # True if file exists and has a size greater than zero.
    done
    sample_idS+=("${sample_id}"); R1[$sample_id]="${r1}"; R2[$sample_id]="${r2}"
done <  "${SAMPLES_TSV}"

(( ${#SAMPLES[@]} > 0 )) || { echo "ERROR: no samples found" >&2; exit 1; } # Counts the number of elements in the final array of samples. If
                                                    # the array is empty, raises an error

# 2. Assemble reads into contigs with MegaHIT. megahit command generates a directory {sample}.megahit_asm 
# which contains checkpoints.txt, done, final.contigs.fa, intermediate_contig/, log, options.json. 
# We rename the final.contigs.fa to add sample name and date
mark_dir="${megahit_output_dir}/.done"
mkdir -p "${mark_dir}"
for SAMPLE_ID in "${SAMPLES[@]}" 
do
    if is_done megahit "${SAMPLE_ID}"; then
        echo "Skipping MegaHIT for ${SAMPLE_ID}"; continue
    fi
    asm="${megahit_output_dir}/${SAMPLE_ID}.megahit_asm"
    rm -rf "$asm"   # megahit refuses to run if the output dir exists
    echo "Running MegaHIT on ${SAMPLE_ID}"
    megahit \
        -1 "${R1[$SAMPLE_ID]}" \
        -2 "${R2[$SAMPLE_ID]}" \
        -o "${asm}" \
        -t "${nthreads_megahit}" 2>&1 |tee "${megahit_logs_dir}"/"${SAMPLE_ID}"_megahit.log
    mv "${asm}/final.contigs.fa" "${asm}/${SAMPLE_ID}_final.contigs.fa"
    mark_done megahit "${SAMPLE_ID}";
done

# --k-min : minimum kmer size (<= 255), must be odd number [default 21] 
# --k-max : maximum kmer size (<= 255), must be odd number [default 141]
# --k-step : increment of kmer size of each iteration (<= 28), must be even number [default 12]

# Tip 1. for ultra complex metagenomics data such as soil, a larger kmin, say 27,
#  is recommended to reduce the complexity of the de Bruijn graph. Quality trimming is also recommended
# Tip 2. for high-depth generic data, large --k-min (25 to 31) is recommended
# Tip 3. smaller --k-step, say 10, is more friendly to low-coverage datasets

# 3. Contig QC with QUAST
mark_dir="${metaquast_megahit_output_dir}/.done"
mkdir -p "${mark_dir}"
echo "$(date +"%F %H:%M:%S")"
for SAMPLE_ID in "${SAMPLES[@]}" 
do
    if is_done metaquast_megahit "${SAMPLE_ID}"; then
        echo "Skipping MetaQuast for ${SAMPLE_ID} contigs"; continue
    fi
    rm -rf "${metaquast_megahit_output_dir}"/metaquast_megahit_"${SAMPLE_ID}"
    MEGAHIT_CONTIGS_FASTA="${megahit_output_dir}/${SAMPLE_ID}.megahit_asm/${SAMPLE_ID}_final.contigs.fa"
    require "${MEGAHIT_CONTIGS_FASTA}"
    echo "Running MetaQuast on ${SAMPLE_ID}" 
    "${metaquast_script_dir}"/metaquast.py \
        "${MEGAHIT_CONTIGS_FASTA}"  \
        -o "${metaquast_megahit_output_dir}"/metaquast_megahit_"${SAMPLE_ID}" \
        --threads "${nthreads_megahit}" \
        2>&1 |tee "${megahit_logs_dir}"/"${SAMPLE_ID}"_metaquast_megahit.log
    mark_done metaquast_megahit "${SAMPLE_ID}";
done

# General usage:
#     python metaquast.py contigs_1 contigs_2 ... -r reference_1,reference_2,reference_3,...
# One can provide several files or directories with multiple reference files inside 
# with -r option. 
# Option -r may be specified multiple times or all references may be specified as
#  a comma-separated list (without spaces!) with a single -r option beforehand

# 4. Extract mapped and unmapped reads
# 4.0 Use BBmap and samtools
conda deactivate
set +u; conda activate qc-tools; set -u

# 4.1 Align reads with bbwrap.sh:
# Then, convert file _mapped_and_unmapped.sam to _mapped_and_unmapped.bam (samtools view -b)
echo "$(date +"%F %H:%M:%S")"
mark_dir="${bbwrap_refs_dir}/.done"
mkdir -p "${mark_dir}"
for SAMPLE_ID in "${SAMPLES[@]}" 
do
    if is_done bbwrap_refs "${SAMPLE_ID}"; then
        echo "Skipping bbwrap.sh for ${SAMPLE_ID}"; continue
    fi
    SAMPLE_ID_REFS="${bbwrap_refs_dir}/${SAMPLE_ID}_refs"
    rm -rf "${SAMPLE_ID_REFS}"
    MEGAHIT_CONTIGS_FASTA="${megahit_output_dir}/${SAMPLE_ID}.megahit_asm/${SAMPLE_ID}_final.contigs.fa"
    require "${MEGAHIT_CONTIGS_FASTA}"
    echo "Running bbwrap.sh on ${SAMPLE_ID}"
    mkdir -p "${SAMPLE_ID_REFS}"
    bbwrap.sh ref="${MEGAHIT_CONTIGS_FASTA}" \
        in="${R1[$SAMPLE_ID]}" \
        in2="${R2[$SAMPLE_ID]}" \
        out="${megahit_aligned_reads_dir}/${SAMPLE_ID}_mapped_and_unmapped.sam" \
        kfilter=22 \
        subfilter=15 \
        path="${SAMPLE_ID_REFS}" \
        maxindel=80 2>&1 |tee "${megahit_logs_dir}"/"${SAMPLE_ID}"_bbwrap.log
    echo "Running samtools view -b on ${SAMPLE_ID}"
    samtools view -b "${megahit_aligned_reads_dir}/${SAMPLE_ID}_mapped_and_unmapped.sam" > \
        "${megahit_aligned_reads_dir}/${SAMPLE_ID}_mapped_and_unmapped.bam" 
    mark_done bbwrap_refs "${SAMPLE_ID}";
done

# 4.2 Output per contig coverage to cov.txt with pileup.sh:
# BBMap generates coverage information by internally using Pileup
echo "$(date +"%F %H:%M:%S")"
mark_dir="${contig_coverages_dir}/.done"
mkdir -p "${mark_dir}"
for SAMPLE_ID in "${SAMPLES[@]}" 
do
    if is_done pileup "${SAMPLE_ID}"; then
        echo "Skipping pileup.sh for ${SAMPLE_ID}"; continue
    fi

    ALIGNED_READS_BAM="${megahit_aligned_reads_dir}/${SAMPLE_ID}_mapped_and_unmapped.bam"
    require "${ALIGNED_READS_BAM}"
    echo "Running pileup.sh on ${SAMPLE_ID}"
    pileup.sh in="${ALIGNED_READS_BAM}"\
        out="${contig_coverages_dir}"/"${SAMPLE_ID}"_cov.txt \
        2>&1 |tee "${megahit_logs_dir}"/"${SAMPLE_ID}"_pileup.log
    mark_done pileup "${SAMPLE_ID}";
done

# 4.3 Extract unmapped reads (SE to unmapped.se.fq and PE to unmapped.pe.fq) (samtools view -u -f4)
echo "$(date +"%F %H:%M:%S")"
mark_dir="${megahit_aligned_reads_dir}/.done"
mkdir -p "${mark_dir}"
for SAMPLE_ID in "${SAMPLES[@]}" 
do
    if is_done bam2fq "${SAMPLE_ID}"; then
        echo "Skipping samtools bam2fq for ${SAMPLE_ID}"; continue
    fi
    ALIGNED_READS_BAM="${megahit_aligned_reads_dir}/${SAMPLE_ID}_mapped_and_unmapped.bam"
    require "${ALIGNED_READS_BAM}"
    echo "Running samtools bam2fq on ${SAMPLE_ID}"
    samtools view -u -f4 "${megahit_aligned_reads_dir}/${SAMPLE_ID}_mapped_and_unmapped.bam" | \
        samtools bam2fq -s "${megahit_aligned_reads_dir}/${SAMPLE_ID}_unmapped.se.fq" - > \
        "${megahit_aligned_reads_dir}/${SAMPLE_ID}_unmapped.pe.fq" ;
    gzip -f -9 --best "${megahit_aligned_reads_dir}/${SAMPLE_ID}_unmapped.se.fq"
    gzip -f -9 --best "${megahit_aligned_reads_dir}/${SAMPLE_ID}_unmapped.pe.fq"
    
    mark_done bam2fq "${SAMPLE_ID}";
done

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
echo "$(date +"%F %H:%M:%S")"
## 4.4 Keep mapped reads then sort
# To get only the mapped reads use the parameter F, which works like -v of grep and skips the 
# alignments for a specific flag (samtools view -b -F 4).
# Then, sort _either_read_mapped.bam file with mapped reads by coordinate (no options used)
# (samtools sort -m)
mark_dir="${megahit_aligned_reads_dir}/.done"
mkdir -p "${mark_dir}"

for SAMPLE_ID in "${SAMPLES[@]}" 
do
    if is_done samtools_sort "${SAMPLE_ID}"; then
        echo "Skipping samtools sort for ${SAMPLE_ID}"; continue
    fi
    ALIGNED_READS_BAM="${megahit_aligned_reads_dir}/${SAMPLE_ID}_mapped_and_unmapped.bam"
    require "${ALIGNED_READS_BAM}"
    echo "Running samtools view -b -F 4 on ${SAMPLE_ID}"
    samtools view -b -F 4 \
        "${ALIGNED_READS_BAM}" \
        > "${megahit_aligned_reads_dir}/${SAMPLE_ID}_either_read_mapped.bam" \
        2>&1 |tee -a "${megahit_logs_dir}"/"${SAMPLE_ID}"_samtools.log
    echo "Running samtools sort on ${SAMPLE_ID}"
    samtools sort -m "${mem_req_megahit_sort}" -@ "${nthreads_megahit_sort}" \
        "${megahit_aligned_reads_dir}/${SAMPLE_ID}_either_read_mapped.bam" \
        -o "${megahit_aligned_reads_dir}"/"${SAMPLE_ID}"_either_read_mapped_sorted.bam 
        
    mark_done samtools_sort "${SAMPLE_ID}";
done

# -b output in BAM format
# -F 4: Exclude unmapped reads. Keep only mapped reads (evem if mate is unmapped)

## 4.5 Filter unmapped read pairs (unmapped to contigs) (samtools view -b -f 12 -F 256)
### SAMtools SAM-flag filter: get unmapped pairs (both reads R1 and R2 unmapped).
# Then, sort _bothReadsUnmapped.bam file by read name ( -n ) to have paired reads next to each other 
### ("${nthreads_megahit}" parallel threads, each using up to 5G memory) (samtools sort)
# Finally, split paired-end reads into separated fastq files .._R1 .._R2 (samtools fastq)
echo "$(date +"%F %H:%M:%S")"
mark_dir="${megahit_aligned_reads_dir}/.done"
mkdir -p "${mark_dir}"

for SAMPLE_ID in "${SAMPLES[@]}" 
do
    if is_done samtools_filter_unmapped "${SAMPLE_ID}"; then
        echo "Skipping samtools filtering for ${SAMPLE_ID}"; continue
    fi
    ALIGNED_READS_BAM="${megahit_aligned_reads_dir}/${SAMPLE_ID}_mapped_and_unmapped.bam"
    require "${ALIGNED_READS_BAM}"
    echo "Running samtools view -b -f 12 -F 256 on ${SAMPLE_ID}" 
    samtools view -b -f 12 -F 256 \
        "${ALIGNED_READS_BAM}" \
        > "${megahit_aligned_reads_dir}/${SAMPLE_ID}_bothReadsUnmapped.bam" 
    echo "Running samtools sort -n and fastq on ${SAMPLE_ID}"
    samtools sort -n -m "${mem_req_megahit_sort}" -@ "${nthreads_megahit_sort}" \
        "${megahit_aligned_reads_dir}/${SAMPLE_ID}_bothReadsUnmapped.bam" \
        -o "${megahit_aligned_reads_dir}/${SAMPLE_ID}_bothReadsUnmapped_sorted.bam" ;
    samtools fastq -@ "${nthreads_megahit_sort}" \
        "${megahit_aligned_reads_dir}/${SAMPLE_ID}_bothReadsUnmapped_sorted.bam" \
        -1 "${megahit_aligned_reads_dir}/${SAMPLE_ID}_unmapped_R1.fastq.gz" \
            -2 "${megahit_aligned_reads_dir}/${SAMPLE_ID}_unmapped_R2.fastq.gz" \
            -0 /dev/null -s /dev/null -n 
    
    mark_done samtools_filter_unmapped "${SAMPLE_ID}";
done
###-f  12    # Extract only ( -f ) alignments with both reads unmapped: <read unmapped><mate unmapped>
###-F 256    # Do not (  -F  ) extract alignments which are: <not primary alignment>

# samtools sort: Sort alignments by leftmost coordinates, 
#  by read name when -n or -N are used,
#  by tag contents with -t,
#  or a minimiser-based collation order with -M.
# -m : Approximately the maximum required memory per thread
# samtools fastq: Converts a BAM or CRAM into either FASTQ or FASTA format
#  depending on the command invoked

# 4.6 Delete sam files to save space
echo "$(date +"%F %H:%M:%S")"
rm -f "${megahit_aligned_reads_dir}"/*_mapped_and_unmapped.sam
rm -f "${megahit_aligned_reads_dir}"/*_bothReadsUnmapped.bam
rm -f "${megahit_aligned_reads_dir}"/*_either_read_mapped.bam
gzip -9 --best "${megahit_aligned_reads_dir}"/*.bam
echo "$(date +"%F %H:%M:%S")"
