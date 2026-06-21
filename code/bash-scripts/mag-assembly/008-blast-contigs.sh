#!/usr/bin/env bash
source $HOME/miniconda3/etc/profile.d/conda.sh
conda activate mag_assembly-tools
set -euo pipefail
shopt -s nullglob
source config/bash/config.sh
# When nullglob is enabled, if a glob pattern does not match any files,
# it expands to nothing (an empty string) instead of returning the pattern itself.
# So, if no matches are found, the script will skip the file

# This script uses BLASTN on SILVA SSU Ref NR99 database on contigs from MEGAHIT
# for a broad taxonomic classification.
echo "$(date +"%F %H:%M:%S")"

# 1. Get the SILVA SSU Ref NR99 database
echo "Downloading SILVA SSU Ref NR99 database"
wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz \
    -O "${silva_db_dir}"/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz
gunzip "${silva_db_dir}"/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz
echo "$(date +"%F %H:%M:%S")"

# 1.1 Build a reference database
echo "Building a BLAST reference database"
makeblastdb -in "${silva_db_dir}"/SILVA_138.2_SSURef_NR99_tax_silva.fasta \
    -parse_seqids \
    -out "${blastdb_path}" \
    -logfile "jobreports/SILVA_138_2_SSURef_NR99_tax_silva-makedblastdb-out.txt" \
    -title "SILVA 138_2 SSURef NR99 db" \
    -dbtype nucl

# 2. Run BLASTN for a broad taxonomic screening
echo "$(date +"%F %H:%M:%S")"
mkdir -p "${blastn_output_dir}"/
echo "Running BLASTN"
for FILE_DIR in "${megahit_output_dir}"/*.megahit_asm 
do 
    echo "$(date +"%F %H:%M:%S")"
    SAMPLE=$(echo "${FILE_DIR}" | sed "s/\.megahit_asm//" )
    base_name=$(basename "$SAMPLE" )
    echo "Running BLASTN on ${base_name} using SILVA 138_2 SSURef NR99 database"
    blastn -query "${megahit_output_dir}"/"${base_name}".megahit_asm/"${base_name}"_final.contigs.fa \
        -db "${blastdb_path}" \
        -out "${blastn_output_dir}"/blastn_"${base_name}"_SILVA_results.tsv \
        -outfmt "6 qacc sacc qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sframe" \
        -num_threads "${nthreads_blast}" \
        2>&1 |tee  "${blastn_output_dir}"/blastn_"${base_name}"_SILVA.log;
done
echo "$(date +"%F %H:%M:%S")"
# qacc means Query accession
# sacc means Subject accession
# qseqid means Query Seq-id
# sseqid means Subject Seq-id
# pident means Percentage of identical matches
# length means Alignment length
# mismatch means Number of mismatches
# gapopen means Number of gap openings
# qstart means Start of alignment in query
# qend means End of alignment in query
# sstart means Start of alignment in subject
# send means End of alignment in subject
# evalue means Expect value
# bitscore means Bit score
# sframe means Subject frame


# 2.1 Extract top BLAST hits with sequence identity > 90% and alignment length > 500 bp
# 1) Filter based on identity and length
# 2) Sort by query ID and descending bitscore
# 3) Pick the top hit per query
# Then, extract qseqid, sseqid, pident, length, bitscore
for FILE in "${blastn_output_dir}"/blastn_*_SILVA_results.tsv
do
    base_name=$(basename "$FILE" )
    base_name=$(echo "${base_name}" | sed "s/_SILVA_results\.tsv//" )
    echo "$(date +"%F %H:%M:%S")"
    # Count all contigs from BLAST
    echo "Counting total number of contigs that had BLAST hits in ${base_name}"
    cut -f 1 "${blastn_output_dir}"/blastn_"${base_name}"_SILVA_results.tsv | \
    uniq  | \
    awk -v sample_id="${base_name}" ' { count++ }
        END {print count " contigs had BLAST hits in sample " sample_id }
        ' 
    # Filter by pident and length (sequence identity > 90% and alignment length > 500 bp), 
    #  and sort by bitscore
    echo "Filtering BLAST hits in ${base_name}"
    awk -F'\t' -v OFS='\t' '$5>90 && $6>500' "${blastn_output_dir}"/blastn_"${base_name}"_SILVA_results.tsv | \
    sort -k1,1 -k14,14nr | awk '!seen[$1]++' > \
        "${blastn_output_dir}"/blastn_"${base_name}"_SILVA_results_filtered.tsv
    #  sort -k1,1 -k14,14nr:
    # sort by qacc (column 1),
    # then by bitscore (column 14), numerically (-n), descending (-r).

    awk -F'\t' -v sample_id="${base_name}" \
        ' { count++ }
        END {print count " contigs had BLAST hits with sequence identity > 90% and alignment length > 500 bp in sample " sample_id }
        ' "${blastn_output_dir}"/blastn_"${base_name}"_SILVA_results_filtered.tsv 
    # Then, extract qseqid, sseqid, pident, length, bitscore
    awk -F'\t' -v OFS='\t' '{print $3, $4, $5, $6, $14}' "${blastn_output_dir}"/blastn_"${base_name}"_SILVA_results_filtered.tsv > \
        "${blastn_output_dir}"/blastn_"${base_name}"_results_for_taxonomy.tsv;
done

# grep "^>" "${silva_db_dir}"/SILVA_138.2_SSURef_NR99_tax_silva.fasta |sed 's/>//' > "${silva_db_dir}"/SILVA_138.2_headers.txt
# Substitute the space between ID and taxonomy with tab
grep "^>" "${silva_db_dir}"/SILVA_138.2_SSURef_NR99_tax_silva.fasta |sed -E 's/^>(\S+)\s+(.*)/\1\t\2/'  > "${silva_db_dir}"/SILVA_138.2_headers.txt

# 2.2 Match BLAST results with SILVA headers
# Final table: 
# contig_id | silva_ref_id | perc_identity | alignment_length | bitscore | taxonomy 
# qseqid, sseqid, pident, length, bitscore
echo "$(date +"%F %H:%M:%S")"
for FILE in "${blastn_output_dir}"/blastn_*_results_for_taxonomy.tsv
do
    base_name=$(basename "$FILE" )
    base_name=$(echo "${base_name}" | sed "s/_results_for_taxonomy\.tsv//" )
    echo "$(date +"%F %H:%M:%S")"
    echo "Matching BLASTN results (sseqid) with SILVA database taxonomy on ${base_name}"
    awk -F'\t' -v OFS='\t' -v sample_name="${base_name}" ' 
        NR==FNR { silva[$1] = substr($0, index($0, $2)); next }
        $2 in silva { print sample_name "\t" $0 "\t" silva[$2] }
        ' "${silva_db_dir}"/SILVA_138.2_headers.txt \
            "${blastn_output_dir}"/blastn_"${base_name}"_results_for_taxonomy.tsv > \
            "${blastn_output_dir}"/blastn_"${base_name}"_blast_with_taxonomy.tsv;
done
echo "$(date +"%F %H:%M:%S")"


# Merge all BLAST taxonomies
cat "${blastn_output_dir}"/blastn_*_blast_with_taxonomy.tsv > \
   "${blastn_output_dir}"/blastn_all_samples.tsv