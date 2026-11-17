#!/usr/bin/env bash
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate mag_assembly-tools
set -euo pipefail
shopt -s nullglob
source config/bash/config.sh
SAMPLES_TSV="${1}"
# When nullglob is enabled, if a glob pattern does not match any files,
# it expands to nothing (an empty string) instead of returning the pattern itself.
# So, if no matches are found, the script will skip the file

# This script quantifies genes.
echo "$(date +"%F %H:%M:%S")"
require() { [[ -s "$1" ]] || { echo "ERROR: required file missing: $1" >&2; exit 1; }; }
# Function that marks a sample as done (in a hidden directory)
is_done()   { [[ -e "${mark_dir}/${2}.${1}.done" ]]; }   # is_done STEP SAMPLE
mark_done() { touch "${mark_dir}/${2}.${1}.done"; }      # mark_done STEP SAMPLE

# 1. Read the input file list
[[ -s "${SAMPLES_TSV}" ]] || { echo "ERROR: metadata $SAMPLES_TSV missing or empty" >&2; exit 1; }
SAMPLES=()
while IFS=$'\t' read -r sample_id \
    || [[ -n "${sample_id}" ]]
do
    [[ "${sample_id}" == "sample_id" || -z "${sample_id}" ]] && continue
    SAMPLES+=("${sample_id}"); 
done <  "${SAMPLES_TSV}"

(( ${#SAMPLES[@]} > 0 )) || { echo "ERROR: no samples found" >&2; exit 1; } # Counts the number of elements in the final array of samples. If
                                                    # the array is empty, raises an error
echo "$(date +"%F %H:%M:%S")"

# 0. Activate conda environment
conda activate qc-tools

# 1. Removing duplicates in the sorted and indexed BAM file
# https://www.biostars.org/p/15818/
echo "$(date +"%F %H:%M:%S")"
for SAMPLE_ID in "${SAMPLES[@]}" 
do
    echo "Adding read group to ${SAMPLE_ID}"
    echo "$(date +"%F %H:%M:%S")"
    samtools addreplacerg \
        -r 'ID:group1' \
        -r 'LB:lib1' \
        -r 'PL:ILLUMINA' \
        -r 'PU:unit1' \
        -r 'SM:sample1' \
        -o "${megahit_aligned_reads_dir}"/"${SAMPLE_ID}"_either_read_mapped_sorted_RG.bam \
        "${megahit_aligned_reads_dir}"/"${SAMPLE_ID}"_either_read_mapped_sorted.bam
    echo "Running Picard MarkDuplicates on ${SAMPLE_ID}"
    echo "$(date +"%F %H:%M:%S")"
    "${java_path}" -Xms2g -Xmx32g -jar "${picard_jar_path}" MarkDuplicates \
        INPUT="${megahit_aligned_reads_dir}"/"${SAMPLE_ID}"_either_read_mapped_sorted_RG.bam \
        OUTPUT="${picard_output_dir}"/"${SAMPLE_ID}"_map.markdup.bam \
        METRICS_FILE="${picard_output_dir}"/"${SAMPLE_ID}"_map.markdup.metrics \
        AS=TRUE \
        VALIDATION_STRINGENCY=LENIENT \
        MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
        REMOVE_DUPLICATES=TRUE \
        2>&1 |tee "${picard_logs_dir}"/"${SAMPLE_ID}"_picard_markduplicates.log
done

echo "$(date +"%F %H:%M:%S")"

# 2. Extract all non-redundant protein locus tags (unique for each sample)
echo "Extracting all non-redundant protein locus tags (unique for each sample)"
grep "^>" \
    "${mmseqs_easy_cluster_output_dir}"/prokka_nr_prot_rep_seq.fasta | \
    sed "s/^>//"| awk -F' ' '{print $1}' > \
    "${mmseqs_easy_cluster_output_dir}"/prokka_nr_prot_ids.txt
echo "Number of non-redundant protein locus tags:" $(wc -l "${mmseqs_clustered_ids_file}")

mmseqs_clustered_ids_file="${mmseqs_easy_cluster_output_dir}"/prokka_nr_prot_ids.txt


# 3. Counting mapped reads per gene
# Convert Prokka output into GTF then count the number of reads mapped to each gene with htseq
# https://metagenomics-workshop.readthedocs.io/en/latest/annotation/quantification.html
for SAMPLE_ID in "${SAMPLES[@]}" 
do
    # Select specific data from Prokka GFF output and convert into GTF. 
    # cut -d uses ";" as a delimiter, while -f1 selects the first field
    # (everything up to ID).
    # In the final GTF file, all features are CDS because we want to count all genes.
    # Real feature types are in the 9th column ($2).
    echo "Converting Prokka output of ${SAMPLE_ID} into GTF"
    grep -v "#" "${FILE_DIR}"/"${SAMPLE_ID}"_pred_prokka.gff | \
        grep "ID=" | cut -f1 -d ';' | sed 's/ID=//g' | cut -f1,3,4,5,7,9 | \
        awk -v OFS='\t' '{print $1,"PROKKA","CDS",$3,$4,".",$5,".","gene_id", "ID="$6";Note="$2}' > \
        "${prokka_gtf_dir}"/"${SAMPLE_ID}"_pred_prokka.map.gtf
    echo "Filtering the GTF file by MMseqs clustered IDs"
    # I want to filter a GTF file using a list of IDs stored in the "ID=XXXX;Note=" section. 
        # getline means read the file and if a line is successfully read, return 1. 0 if awk reached 
        # end of the file. -1 if there is an error like permission denied or file doesnt exist.

        # match($10, /ID=([^;]+)/, prokka_gtf_ids) will search string for the longest, leftmost substring matched by the 
        # regular expression ID=([^;]+) and return the character position (index) at which that substring begins
        # (one, if it starts at the beginning of string). 
        # Extract all the IDs in the gtf field 10 and store in an array prokka_gtf_ids.
        # match(string, regexp [, array])

        # Regex explanation: The regexp argument may be either a regexp constant (/…/) or a string constant ("…").
        # First, match /ID=/ (literal match). Then, /([^;]+)/ captures one or more characters that are not semicolons 
        # So, (...) is the capturing group: everything inside is stored as a match in the array prokka_gtf_ids.
        # We use square brackets [...] to define character class: [abc] means "a single character that is a, b, or c"
        # Here, [^;] means "a single character that is not a semicolon".
        # + means one or more character.

        # If the prokka_gtf_ids entry is found in indices of mmseqs_ids, print.
        # By the way, we use prokka_gtf_ids[1] because index 1 is the first captured group (if there are other matches between
        # ID and semicolon in a string, they would be in prokka_gtf_ids[2] or prokka_gtf_ids[3], or others)

    awk -F'\t' -v mmseqs_file="${mmseqs_clustered_ids_file}" \
    'BEGIN {
        while ((getline < mmseqs_file ) > 0) {
        mmseqs_ids[$1] = 1
        } 
    }
    {        
        match($10, /ID=([^;]+)/, prokka_gtf_ids)

        if (prokka_gtf_ids[1] in mmseqs_ids)
            print $0
    }
    ' "${prokka_gtf_dir}"/"${SAMPLE_ID}"_pred_prokka.map.gtf > \
        "${prokka_gtf_dir}"/"${SAMPLE_ID}"_pred_prokka_mmseqs_filtered.map.gtf
    echo "Number of non-redundant protein locus tags in" "${SAMPLE_ID}" \
    $(wc -l "${prokka_gtf_dir}"/"${SAMPLE_ID}"_pred_prokka_mmseqs_filtered.map.gtf)
    # Counting the number of reads mapped to genes with htseq. Here we have to tell
    # htseq that the file is sorted by alignment coordinate -r pos.
    # The output file has two columns, the first contains the gene names and the second
    # the number of reads mapped to each gene. The last 5 lines gives you some summary 
    # information from htseq: __no_feature, __ambiguous, __too_low_aQual, __not_aligned,
    # __alignment_not_unique
    echo "Quantifying the number of reads mapped to genes with htseq on ${SAMPLE_ID}"
    echo "$(date +"%F %H:%M:%S")"
    htseq-count -r pos -t CDS -f bam "${picard_output_dir}"/"${SAMPLE_ID}"_map.markdup.bam \
        "${prokka_gtf_dir}"/"${SAMPLE_ID}"_pred_prokka.map.gtf > \
        "${htseq_gene_count_dir}"/"${SAMPLE_ID}"_prokka.count \
        2>&1 |tee "${htseq_logs_dir}"/"${SAMPLE_ID}"_htseq_count.log
    grep "ID" "${htseq_gene_count_dir}"/"${SAMPLE_ID}"_prokka.count | \
    sed 's/ID=//g' > "${htseq_gene_count_dir}"/"${SAMPLE_ID}"_prokka_filtered.count
    
    # Combine GTF data into a TSV file: CDS, rRNA, tRNA, partial genes, pseudogenes, etc.
    # Select contig ID ($1), gene start($4) and gene end ($5), and gene ID ($10)
    # But gene ID column also has gene type (Note=...)
    echo "$(date +"%F %H:%M:%S")"
    cut -f1,4,5,10 "${FILE}" |sed 's/ID=//g' | \
    awk -F'\t' 'BEGIN{OFS="\t"} {
            split($4, arr, /;Note=/)
            if (length(arr) == 2) {
                $4 = arr[1]
                $5 = arr[2]
            }
        print
    }
    ' "${prokka_gtf_dir}"/"${SAMPLE_ID}"_pred_prokka.map.gtf \
        | gawk -v sample_name="${SAMPLE_ID}"  '{print sample_name,$1,$4,$3-$2+1,$5}'| tr ' ' '\t' > \
    "${prokka_contig_gene_maps_dir}"/"${SAMPLE_ID}"_contig_gene_mapping.tsv;
done
echo "$(date +"%F %H:%M:%S")"

# All mappings (including redundancies): 3847253 entries, but only 3771950 CDS
cat "${prokka_contig_gene_maps_dir}"/*_contig_gene_mapping.tsv > \
   "${prokka_contig_gene_maps_dir}"/all_contig_gene_maps.tsv
# 3847253 entries with counts
cat "${htseq_gene_count_dir}"/*_prokka_filtered.count > \
   "${htseq_gene_count_dir}"/all_prokka_counts.tsv


# 4. Normalizing to Transcripts Per Million (TPM)
# The gene lengths we can get from the GTF file that you used with htseq:
for SAMPLE_ID in "${SAMPLES[@]}" 
do
    # Here we extract only the start (#4), stop(#5) and gene name fields(#10) from the file, 
    # then remove the 'ID=' and ';Note= ...' string, print the gene name first followed by the 
    # length of the gene, change the separator to tab and store the results in the .genelengths file.
    cut -f4,5,10 "${prokka_gtf_dir}"/"${SAMPLE_ID}"_pred_prokka.map.gtf | \
        sed 's/ID=//g' |sed 's/\;.*//g' | gawk '{print $3,$2-$1+1}' | tr ' ' '\t' > \
        "${gene_lengths_dir}"/"${base_name}".genelengths
        # Now we can calculate TPM values using the tpm_table.py script:
    python3 "${python_utils_scripts_dir}/tpm_table.py" \
        -n "${SAMPLE_ID}" \
        -c "${htseq_gene_count_dir}"/"${SAMPLE_ID}"_prokka_filtered.count \
        -i <(echo -e ""${SAMPLE_ID}"\t150") \
        -l "${gene_lengths_dir}"/"${SAMPLE_ID}".genelengths > \
            "${tpm_files_dir}"/"${SAMPLE_ID}".tpm;
    # We now have coverage values for all genes predicted and annotated by the PROKKA pipeline;

done

