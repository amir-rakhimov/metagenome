#!/usr/bin/env bash
#SBATCH -t 360:00:00
#SBATCH -N 4
#SBATCH -n 20
#SBATCH --mem-per-cpu 8g
#SBATCH -J 20250809_19-19-search-cazymes-blast
#SBATCH --output jobreports/20250809_19-19-search-cazymes-blast-out.txt
#SBATCH --error jobreports/20250809_19-19-search-cazymes-blast-out.txt
source ~/miniconda3/etc/profile.d/conda.sh
date_var=$(date -I|sed 's/-//g')
time_var=$(date +%T |sed 's/:/_/g' )
date_time=${date_var}_${time_var}
start_date_time=$(date +"%F %H:%M:%S")
# tblastn_date_time="${date_time}"
# blastp_date_time=20250729_13_42_31
blastp_date_time="${date_time}"
blastp_output_dir=output/mag_assembly/blastp_output
# tblastn_output_dir=output/mag_assembly/tblastn_output
blastdb_path=data/blastdbs/cazy_db_pep_db/cazy_db_pep
# blastdb_path=data/blastdbs/B_coli_genome_nt_db/Balantioides_coli_genome_nt_db
nthreads=20

conda activate qc-tools

# Download information about taxonomy of reference CAZymes
# mkdir data/cazy_data
# wget https://www.cazy.org/IMG/cazy_data/cazy_data.zip -O data/cazy_data/cazy_data.zip
# unzip data/cazy_data/cazy_data.zip -d data/cazy_data/

# Download the peptide sequences of reference CAZymes
# wget https://bcb.unl.edu/dbCAN2/download/CAZyDB.07142024.fa -O data/cazy_data/CAZyDB.07142024.fa

# Extract peptide sequences of CAZymes in the tsv file for tblastn
# seqkit grep -f output/rtables/cazymes-for-blast.tsv \
#     output/mag_assembly/prokka_output/20250626_22_11_43_all_proteins.faa \
#     -o output/fasta/cazymes-for-blast.faa

# seqkit grep -f \
#     <(cut -f 1 output/mag_assembly/dbcan_output/20250626_22_11_43_nr_dbcan_proteins/overview.txt |tail -n +2) \
#     output/mag_assembly/prokka_output/20250626_22_11_43_all_proteins.faa \
#     -o output/fasta/all-cazymes.faa

# Change the IDs to seq_00001 because current IDs are too long. The original ID follows after space
# The sprintf function in the awk allows you to format strings in a customizable way.
# Here, %08d pads the integer with zeros to ensure an 8-digit representation.

# awk '/^>/ {
#     header=$0;
#     id = sprintf("seq%08d", ++i);
#     print ">" id " " substr(header, 2);  # keep full header after the space
#     next
# }
# {
#     print
# }' data/cazy_data/CAZyDB.07142024.fa >data/cazy_data/CAZyDB.07142024.seq.fa


# Use BLAST to find CAZymes in the reference CAZy database
# mkdir -p data/blastdbs/cazy_db_pep_db
# Build a database
# makeblastdb -in data/cazy_data/CAZyDB.07142024.seq.fa \
#  -parse_seqids \
#  -out "${blastdb_path}" \
#  -logfile "jobreports/"${blastp_date_time}"-cazy_db_pep-makedblastdb-out.txt" \
#  -title "CAZy database peptide db" \
#  -dbtype prot


mkdir -p "${blastp_output_dir}"/"${blastp_date_time}"
blastp -db "${blastdb_path}" \
    -query output/fasta/all-cazymes.faa \
    -out "${blastp_output_dir}"/"${blastp_date_time}"/"${blastp_date_time}"_blastp_cazy_db_cazymes_results.tsv \
    -outfmt "6 qacc sacc qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sframe" \
    -num_threads "${nthreads}" \
    2>&1 |tee  "${blastp_output_dir}"/"${blastp_date_time}"/"${blastp_date_time}"_blastp_cazy_db_cazymes.log


# Select top hits by bitscore
#  sort -k1,1 -k14,14nr:
# sort by qacc (column 1),
# then by bitscore (column 14), numerically (-n), descending (-r).

sort -k1,1 -k14,14nr "${blastp_output_dir}"/"${blastp_date_time}"/"${blastp_date_time}"_blastp_cazy_db_cazymes_results.tsv | \
    awk  '!seen[$1]++' > "${blastp_output_dir}"/"${blastp_date_time}"/"${blastp_date_time}"_blastp_cazy_db_cazymes_top_hits.tsv


# Search top hits
# First, find the sequences that match the second column of top hits (IDs in the blastdb).
# Then, get the IDs by matching ">" symbol.
# Then, split the string by space. The second column is the original ID. Then, split that ID by pipe (|) and take the 
# first column which is the Uniprot/Uniparc/etc ID of the protein. Store it in the txt file

seqkit grep -f <(cut -f 2 "${blastp_output_dir}"/"${blastp_date_time}"/"${blastp_date_time}"_blastp_cazy_db_cazymes_top_hits.tsv) \
    data/cazy_data/CAZyDB.07142024.seq.fa |grep ">" |awk -F ' ' '{print $2}' |awk -F '|' '{print $1}' >\
    "${blastp_output_dir}"/"${blastp_date_time}"/"${blastp_date_time}"_blastp_cazy_db_cazymes_top_ids.txt

# Search all hits
# seqkit grep -f  <(cut -f 2 "${blastp_output_dir}"/"${blastp_date_time}"/"${blastp_date_time}"_blastp_cazy_db_cazymes_results.tsv) \
#     data/cazy_data/CAZyDB.07142024.seq.fa |grep ">" |awk -F ' ' '{print $2}' |awk -F '|' '{print $1}' >\
#     "${blastp_output_dir}"/"${blastp_date_time}"/"${blastp_date_time}"_blastp_cazy_db_cazymes_all_ids.txt
# Find all hits in the db
# grep -f "${blastp_output_dir}"/"${blastp_date_time}"/"${blastp_date_time}"_blastp_cazy_db_cazymes_all_ids.txt \
#     data/cazy_data/cazy_data_20250707.txt


# Extract taxonomy information
# Take the ids from the blastp_cazy_db_cazymes_top_ids.txt and store them in an associative array
# cazy_data_20250707.txt is the target. We check its 4th column (contains sequence IDs in the fasta file).
# ids[$1]; next stores each ID in the associative array and moves to the next line
# $4 in ids: Checks the existence of ID in the 4th column of target file
awk -F'\t' -v OFS='\t' 'NR==FNR {ids[$1]; next} $4 in ids' \
    "${blastp_output_dir}"/"${blastp_date_time}"/"${blastp_date_time}"_blastp_cazy_db_cazymes_top_ids.txt \
    data/cazy_data/cazy_data_20250707.txt |sort |uniq > \
    "${blastp_output_dir}"/"${blastp_date_time}"/"${blastp_date_time}"_blastp_cazy_db_cazymes_top_taxa.tsv

# The BLASTP output uses our new IDs, so let's join new IDs with original ones
grep ">" data/cazy_data/CAZyDB.07142024.seq.fa | sed "s/>//"| sed "s/|/ /" | awk -F' ' -v  OFS='\t' '{$1=$1; print}' > \
    data/cazy_data/CAZyDB.07142024.seq.sep.tsv

# Join BLASTP results with the reference table
awk -F'\t' -v OFS='\t' 'NR==FNR {map[$1]=$0; next} 
    {
        if($2 in map){
            # output file1 row + rest of matching row from file2 (excluding key)
            split(map[$2], a, FS)
            # join everything: file1 fields + file2 fields starting from $2
            printf "%s", $0
            for (i=2; i<=length(a); i++) printf OFS a[i]
            print ""
        } else {
            # no match: just print file1 row
            print $0
        }
    }' data/cazy_data/CAZyDB.07142024.seq.sep.tsv \
    "${blastp_output_dir}"/"${blastp_date_time}"/"${blastp_date_time}"_blastp_cazy_db_cazymes_top_hits.tsv > \
    "${blastp_output_dir}"/"${blastp_date_time}"/"${blastp_date_time}"_blastp_cazy_db_cazymes_top_ref_table.tsv

# Join BLASTP output with taxonomy information
awk -F'\t' -v OFS='\t' 'NR==FNR { map[$4] = $0; next } 
    { 
        if($16 in map){
            # output file1 row + rest of matching row from file2 (excluding key)
            split(map[$16], a, FS)
            # join everything: file1 fields + file2 fields starting from $2
            printf "%s", $0
            for (i=2; i<=3; i++) printf OFS a[i]
            print ""
        } else {
            # no match: just print file1 row
            print $0
        }
    }' \
    "${blastp_output_dir}"/"${blastp_date_time}"/"${blastp_date_time}"_blastp_cazy_db_cazymes_top_taxa.tsv \
    "${blastp_output_dir}"/"${blastp_date_time}"/"${blastp_date_time}"_blastp_cazy_db_cazymes_top_ref_table.tsv > \
     "${blastp_output_dir}"/"${blastp_date_time}"/"${blastp_date_time}"_blastp_cazy_db_cazymes_top_hits_with_taxa.tsv
    # "${blastp_output_dir}"/"${blastp_date_time}"/"${blastp_date_time}"_blastp_cazy_db_cazymes_results.tsv \
