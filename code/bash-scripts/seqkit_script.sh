#!/usr/bin/env bash
shopt -s nullglob
# date_var=$(date -I|sed 's/-//g')
# time_var=$(date +%T |sed 's/:/_/g' )
# date_time=${date_var}_${time_var}
date_time=20250211_22_55_27
project_home_dir=~/projects/metagenome
seqkit_output_dir=output/mag_assembly/seqkit_output
metabat2_output_dir=output/mag_assembly/metabat2_output
conda activate qc-tools
# Calculate statistics for each contig of each bin
for FILE_DIR in "${metabat2_output_dir}"/"${date_time}"*
do 
 SAMPLE=$(echo "${FILE_DIR}" | sed "s/${date_time}_//")
 base_name=$(basename "$SAMPLE" )
 mkdir -p "${seqkit_output_dir}"/"${date_time}"_"${base_name}"_seqkit_stats
 seqkit_stats_out_dir="${seqkit_output_dir}"/"${date_time}"_"${base_name}"_seqkit_stats
 echo "Calculating statistics on ${base_name}" 
  for BIN_FASTA_FILE in "${metabat2_output_dir}"/"${date_time}"_"${base_name}"/*.fa
  do
    BIN_ID=$(basename "${BIN_FASTA_FILE}" )
    BIN_ID=$(echo "${BIN_ID}" | sed "s/${date_time}_//"| sed "s/\.fa//")
    seqkit fx2tab -l -g -n -i -H -B AT "${BIN_FASTA_FILE}" | \
     csvtk -t -C '&'  rename -f "#id",length -n contig_id,contig_length \
     > "${seqkit_stats_out_dir}"/"${BIN_ID}"_seqkit.tsv;
  done
 echo "Done";
done

-l, --length    print sequence length
-g, --gc    print GC content (%)
-n, --name    only print names (no sequences and qualities)
-i, --only-id    print ID instead of full head
-H, --header-line    print header line (#id     length  GC)
-B, --base-content strings    print base content (%). (case ignored, multiple values supported) e.g. -B AT -B N

csvtk:
  1. By default, csvtk assumes input files have header row, if not, switch flag "-H" on.
  2. By default, csvtk handles CSV files, use flag "-t" for tab-delimited files.

rename:
  -f, --fields string   select only these fields. e.g -f 1,2 or -f columnA,columnB
  -F, --fuzzy-fields    using fuzzy fields, e.g., -F -f "*name" or -F -f "id123*"
  -n, --names string    comma separated new names



# Merge tables
# awk 'BEGIN {FS=OFS="\t"} 
#     FNR==1 && NR==1 {print "Sample", "Bin", $0} 
#     FNR>2 {split(FILENAME, parts, "/"); filename=parts[length(parts)]; \
#       sub(/_seqkit.tsv$/, "", filename); 
#       split(filename, parts, "_"); sample_id=parts[1]; bin_id=parts[2];
#       sub(/bin./, "", bin_id); \
#       print sample_id, bin_id, $0}' \
#         "${seqkit_output_dir}"/"${date_time}"_2D14_seqkit_stats/*.tsv > \
#         "${seqkit_output_dir}"/2D14_combined_stats.tsv

# Join per-bin stats into a per-sample file: all bins of each sample into one table. 
# Add a column with sample name and bin ID.
for FILE_DIR in "${seqkit_output_dir}"/"${date_time}"*_seqkit_stats
do 
 SAMPLE=$(echo "${FILE_DIR}" | sed "s/${date_time}_//" |sed "s/_seqkit_stats//")
 base_name=$(basename "$SAMPLE" )
 seqkit_stats_out_dir="${seqkit_output_dir}"/"${date_time}"_"${base_name}"_seqkit_stats
 stats_files=("${seqkit_stats_out_dir}"/*.tsv) # Store matching files in an array
 if [[ ${#stats_files[@]} -eq 0 ]]; then
    echo "Error: No *.tsv files found in the directory ${seqkit_stats_out_dir}"
    continue
 else
    echo "Combining tables ${base_name}"
    awk 'BEGIN {FS=OFS="\t"} 
        FNR==1 && NR==1 {print "Sample", "Bin", $0} 
        FNR>1 {split(FILENAME, parts, "/"); filename=parts[length(parts)]; \
            sub(/_seqkit.tsv$/, "", filename); 
            split(filename, parts, "_"); sample_id=parts[1]; bin_id=parts[2];
            sub(/bin./, "", bin_id); \
            print sample_id, bin_id, $0}' \
            "${seqkit_stats_out_dir}"/*.tsv > \
        "${seqkit_output_dir}"/"${date_time}"_"${base_name}"_combined_stats.tsv
 fi
 echo "Done";
done

# Join everything. The ID column is a concatenation of Sample, Bin, and contig ID to avoid duplicates
awk 'BEGIN {FS=OFS="\t"} 
    FNR==1 && NR==1 {printf  "%s\t%s\t%s", "Sample", "Bin", "contig_id"; 
      for (i = 4; i <= NF; i++) printf "\t%s", $i;
      print "";} 
    FNR>1 {printf "%s\t%s\t%s", $1, $2, $1"_"$2"_"$3; for (i = 4; i <= NF; i++) printf "\t%s", $i;
      print "";}' \
        "${seqkit_output_dir}"/*_combined_stats.tsv > \
    "${seqkit_output_dir}"/"${date_time}"_all_samples_stats.tsv

# Sanity check
# wc_count=0
# for FILE in "${seqkit_output_dir}"/*_combined_stats.tsv
# do 
#  count=$(grep -v  "Sample" $FILE  |wc -l )
#  wc_count=$((wc_count + count))
# done
# echo $wc_count

# Extract certain contig by its id
seqkit grep -n output/mag_assembly/metabat2_output/20250211_22_55_27_2D14/20250211_22_55_27_2D14_bin.13.fa \
 -p "k141_109001" |grep -v "^>"| tr -d '\n' > output/mag_assembly/2D14_13_k141_109001.fasta

seqkit grep -n output/mag_assembly/metabat2_output/20250211_22_55_27_2D14/20250211_22_55_27_2D14_bin.10.fa \
 -p "k141_237915" | seqkit subseq -r 1:10000 > output/mag_assembly/metabat2_output/2D14_bin.10.k141_237915.fasta

seqkit grep -n output/mag_assembly/megahit_output/20250211_22_55_27_2D14.megahit_asm/20250211_22_55_27_2D14_final.contigs.fa \
 -p "k141_0 flag=1 multi=1.0000 len=364"  > output/mag_assembly/metabat2_output/2D14_k141_0.fa