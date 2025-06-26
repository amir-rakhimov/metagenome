#!/usr/bin/env bash
source ~/miniconda3/etc/profile.d/conda.sh
shopt -s nullglob
date_var=$(date -I|sed 's/-//g')
time_var=$(date +%T |sed 's/:/_/g' )
date_time=${date_var}_${time_var}
# metabat2_date_time=20250417_23_21_29
metabat2_date_time=20250612_13_37_47
drep_date_time="${date_time}"
seqkit_stats_date_time="${date_time}"
coverm_date_time="${date_time}"
nthreads=40
nthreads_sort=35
mem_req=8G
mem_req_sort=4G
project_home_dir=~/projects/metagenome
bowtie2_decontam_fastq_dir=data/bowtie2_decontam_fastq
output_dir=~/projects/metagenome/output/mag_assembly
metabat2_output_dir=output/mag_assembly/metabat2_output
drep_output_dir=output/mag_assembly/drep_output
seqkit_output_dir=output/mag_assembly/seqkit_output
coverm_output_dir=output/mag_assembly/coverm_output

cd ${project_home_dir}
mkdir -p output/mag_assembly/seqkit_output
mkdir -p output/mag_assembly/coverm_output
conda activate mag_assembly-tools
# CoverM calculates abundance of dereplicated bins
## At species level (95% ANI)
for FILE in "${bowtie2_decontam_fastq_dir}"/*decontam_R1.fastq.gz
do 
 SAMPLE=$(echo "${FILE}" | sed "s/_trim_decontam_R1\.fastq\.gz//")
 base_name=$(basename "$SAMPLE" )
 echo "Running CoverM on ${base_name}"
 coverm genome \
  --coupled "${bowtie2_decontam_fastq_dir}"/"${base_name}"_trim_decontam_R1.fastq.gz \
    "${bowtie2_decontam_fastq_dir}"/"${base_name}"_trim_decontam_R2.fastq.gz \
  --genome-fasta-files \
  "${drep_output_dir}"/"${drep_date_time}"_sa_95perc/dereplicated_genomes/"${drep_date_time}"_"${base_name}"_bin.*.fa \
  -t "${nthreads}" \
  -m mean relative_abundance covered_fraction \
  -o "${coverm_output_dir}"/"${coverm_date_time}"_"${base_name}"_coverm_output.tsv |
  2>&1 |tee "${coverm_output_dir}"/"${coverm_date_time}"_"${base_name}"_coverm_report.txt



 

done

coverm genome \
  --coupled sample_1.1.fq.gz sample_1.2.fq.gz \
  --genome-fasta-files \
  "${drep_output_dir}"/"${drep_date_time}"_sa_95perc/dereplicated_genomes/*fa \
  -t "${nthreads}" \
  -m mean relative_abundance covered_fraction \
  -o "${coverm_output_dir}"/"${coverm_date_time}"_output_coverm.tsv


## At strain level (99% ANI)
coverm genome \
   "${drep_output_dir}"/"${drep_date_time}"_sa_99perc/dereplicated_genomes/*fa



conda activate qc-tools
# Calculate statistics for each contig in each bin: contig IDs, contig lengths, GC %, 
# AT %. Output is a tab-separated file for each bin which is located in a
# directory that corresponds to the sample where the bin was generated.
for FILE_DIR in "${metabat2_output_dir}"/"${metabat2_date_time}"*bins
do 
 SAMPLE=$(echo "${FILE_DIR}" | sed "s/${metabat2_date_time}_//" | sed "s/_bins//")
 base_name=$(basename "$SAMPLE" )
 mkdir -p "${seqkit_output_dir}"/"${metabat2_date_time}"_"${base_name}"_seqkit_stats
 seqkit_stats_out_dir="${seqkit_output_dir}"/"${metabat2_date_time}"_"${base_name}"_seqkit_stats
 echo "Calculating statistics on ${base_name}" 
  for BIN_FASTA_FILE in "${metabat2_output_dir}"/"${metabat2_date_time}"_"${base_name}"_bins/*.fa
  do
    BIN_ID=$(basename "${BIN_FASTA_FILE}" )
    BIN_ID=$(echo "${BIN_ID}" | sed "s/${metabat2_date_time}_//"| sed "s/\.fa//")
    seqkit fx2tab -l -g -n -i -H -B AT "${BIN_FASTA_FILE}" | \
     csvtk -t -C '&'  rename -f "#id",length -n contig_id,contig_length \
     > "${seqkit_stats_out_dir}"/"${BIN_ID}"_seqkit.tsv;
  done
 echo "Done";
done

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

for FILE_DIR in "${seqkit_output_dir}"/"${metabat2_date_time}"*_seqkit_stats
do 
 SAMPLE=$(echo "${FILE_DIR}" | sed "s/${metabat2_date_time}_//" |sed "s/_seqkit_stats//")
 base_name=$(basename "$SAMPLE" )
 seqkit_stats_out_dir="${seqkit_output_dir}"/"${metabat2_date_time}"_"${base_name}"_seqkit_stats
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
        "${seqkit_output_dir}"/"${seqkit_stats_date_time}"_"${base_name}"_combined_stats.tsv
 fi
 echo "Done";
done

# Join all sample stats into one tab-separated file. 
# Output is a tab-separated file with stats for all bins in each sample (sample name,
# bin number, contig ID, contig length, GC content, AT content). All samples are 
# gathered in one file. To avoid duplications in contig IDs, I changed the ID column 
# to a concatenation of Sample, Bin, and contig ID.

# If the row number is 1 (header), we declare that 
# the first three columns of the header are separated by tab (Sample, Bin, contig_id), 
# then print these columns. Then, we print all the remaining columns from the bin stats file 
# (contig_length, GC%, AT%, etc), which are also separated by tab 
# (for (i = 4; i <= NF; i++) printf "\t%s", $i;).
# We use a for loop because we don't know how many columns are remaining.
# Then, we print a new line.

# If the row number is >1, we declare that 
# the first three columns are separated by tab ("%s\t%s\t%s"), then we 
# print their contents: column 1 of the bin stats file (Sample), column 2 (bin),
# concatenation of columns 1, 2, and 3 (Sample, bin, contig_id), and the 
# remaining columns that will be separated by tab. The new column 3 is the new contig_id.
# Then, print the new line and move to the next row.
awk 'BEGIN {FS=OFS="\t"} 
    FNR==1 && NR==1 {printf  "%s\t%s\t%s", "Sample", "Bin", "contig_id"; 
      for (i = 4; i <= NF; i++) printf "\t%s", $i;
      print "";} 
    FNR>1 {printf "%s\t%s\t%s", $1, $2, $1"_"$2"_"$3; for (i = 4; i <= NF; i++) printf "\t%s", $i;
      print "";}' \
        "${seqkit_output_dir}"/*_combined_stats.tsv > \
    "${seqkit_output_dir}"/"${seqkit_stats_date_time}"_all_samples_stats.tsv

# Sanity check
wc_count=0
for FILE in "${seqkit_output_dir}"/*_combined_stats.tsv
do 
 count=$(grep -v  "Sample" $FILE  |wc -l )
 wc_count=$((wc_count + count))
done
echo $wc_count
wc "${seqkit_output_dir}"/"${seqkit_stats_date_time}"_all_samples_stats.tsv

# Calculate average GC% and AT% per bin and select bins with high average AT%
# head -n 390 20250303_17_54_26_all_samples_stats.tsv >stats_head.tsv

# Sample is column 1, bin is column 2, GC% is column 5, AT% is column 6
# The gsub(/"/, "") function used to remove all double quotes (") from each 
# line before processing the data because TSV and CSV files often enclose 
# their columns with double quotes. Removing double quotes ensures correct
# processing.
# /"/: The regular expression that matches double quotes (").

# I create sample_bin_GC and sample_bin_AT, which are associative arrays 
# (also called hash tables or dictionaries in other programming languages) in awk.
# In awk, arrays can have string keys, not just numeric indices.
# The key in both sample_bin_GC and sample_bin_AT is key = $1 "_" $2 
# (which combines Sample and bin into a single string).
# The values stored in these arrays (count) are the cumulative sums
#  of GC% ($5) and AT% ($6), respectively.
# In the end, we print the key (e.g. SampleA_bin1), the mean GC%, and the mean AT%
awk 'BEGIN {FS=OFS="\t"} 
    FNR==1 && NR==1 {printf  "%s\t%s\t%s", "Sample_Bin", "MeanGC", "MeanAT"; 
      print "";} 

      FNR >1 {
        gsub(/"/, "")
        key = $1 "_" $2
        sample_bin_GC[key] += $5
        sample_bin_AT[key] += $6
        ++count[key]
      }
      END {
        for (key in sample_bin_GC)
          print key, sample_bin_GC[key]/count[key], sample_bin_AT[key]/count[key]
      }' "${seqkit_output_dir}"/"${seqkit_stats_date_time}"_all_samples_stats.tsv > \
        "${seqkit_output_dir}"/"${seqkit_stats_date_time}"_average_gc_at_bins.tsv

# Filter to keep rows with low GC (below 30) (column 2)
awk -v threshold=30 'NR==1 || $2 < threshold'  "${seqkit_output_dir}"/"${seqkit_stats_date_time}"_average_gc_at_bins.tsv \
        > "${seqkit_output_dir}"/"${seqkit_stats_date_time}"_low_gc_bins.tsv


# Extract certain contig by its id
# seqkit grep -n output/mag_assembly/metabat2_output/20250211_22_55_27_2D14/20250211_22_55_27_2D14_bin.13.fa \
#  -p "k141_109001" |grep -v "^>"| tr -d '\n' > output/mag_assembly/2D14_13_k141_109001.fasta

# seqkit grep -n output/mag_assembly/metabat2_output/20250211_22_55_27_2D14/20250211_22_55_27_2D14_bin.10.fa \
#  -p "k141_237915" | seqkit subseq -r 1:10000 > output/mag_assembly/metabat2_output/2D14_bin.10.k141_237915.fasta

# seqkit grep -n output/mag_assembly/megahit_output/20250211_22_55_27_2D14.megahit_asm/20250211_22_55_27_2D14_final.contigs.fa \
#  -p "k141_0 flag=1 multi=1.0000 len=364"  > output/mag_assembly/metabat2_output/2D14_k141_0.fa