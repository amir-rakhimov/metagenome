#!/usr/bin/env bash
#SBATCH -N 4
#SBATCH -n 49
#SBATCH --mem-per-cpu 8g
#SBATCH -J 20250911_17-50-kraken2-pipeline
#SBATCH --output jobreports/20250911_17-50-kraken2-pipeline-out.txt
#SBATCH --error jobreports/20250911_17-50-kraken2-pipeline-out.txt
#I am requesting 4 nodes containing 49 CPUs, with 8 GB memory per CPU. Total: 392 GB

source ~/miniconda3/etc/profile.d/conda.sh
MANIFEST=$1
source "${MANIFEST}"
source config/config.sh
echo $RUN_ID
shopt -s nullglob
# When nullglob is enabled, if a glob pattern does not match any files,
# it expands to nothing (an empty string) instead of returning the pattern itself.
# So, if no matches are found, the script will skip the file

######
# This script performs the following:

# 1. Classifies decontaminated samples using Kraken2.

# Input: decontaminated paired-end FASTQ files in the `data/bowtie2_decontam_fastq` directory,
# Kraken2 reference database path.

# Output: 
# - Kraken2 output files (.kraken2), without minimizers, saved in the 
# `output/kraken2_pipeline/kraken2_output` directory. However,
# it's possible to add minimizers if needed. The output file name format is 
# `${sample_name}_no_minimizer_data.kraken2`. Gzip-compressed.
# - Kraken2 reports (.k2report) with taxonomic classification saved in the 
# `output/kraken2_pipeline/kraken2_reports` directory. The output file name format is 
# `${sample_name}_no_minimizer_data.k2report`. Gzip-compressed.
# - FASTQ files with classified sequences (read 1 and read 2 are in separate files) saved in
# the `output/kraken2_pipeline/kraken2_classified_reads` directory.
# The output file name format is `${sample_name}_classified__1.fq` for read 1
# and  `${sample_name}_classified__2.fq` for read 2.
# - Stdout file saved in the general directory (`~/projects/metagenome`) as 
# `kraken2_stdout_no_minimizer_data.txt`.

# 2. Runs Bracken for abundance estimation of microbiome samples

# Input: Kraken2 reports in the `output/kraken2_pipeline/kraken2_reports` directory, 
# Kraken2 reference database path.

# Output: 
# - Bracken output files (.bracken) saved in the 
# `output/kraken2_pipeline/bracken_output` directory. The output file name format is
# `${sample_name}_no_minimizer_data.bracken`.
# - Bracken reports (.breport) saved in the 
# `output/kraken2_pipeline/bracken_reports` directory. The output file name format is
# `${sample_name}_no_minimizer_data.breport`.
# - Stdout file saved in the general directory (`~/projects/metagenome`) as
# `bracken_stdout_no_minimizer_data.txt`.

# 3. Generates TXT files for Krona plots

# Input: Bracken reports in the `output/kraken2_pipeline/bracken_reports` directory.

# Output: TXT files for Krona generated from Bracken reports (one TXT file per sample) 
# and saved in the `output/kraken2_pipeline/bracken_krona_txt` directory.

# 4. Generates Krona plots

# Input: TXT files for Krona saved in `output/kraken2_pipeline/bracken_krona_txt` directory.

# Output: HTML files with Krona plots saved in `output/kraken2_pipeline/krona_html` directory.

# 5. Extract reads that were mapped to selected taxa: Cryptomeria japonica (taxonomy id 3369),
# Glycine soja (taxid 3848), Ipomoea triloba (taxid 35885), and Macadamia integrifolia (taxid 60698).

# Input: FASTQ files with classified sequences (read 1 and read 2 are in separate files) in 
# the `output/kraken2_pipeline/kraken2_classified_reads` directory.

# Output: FASTA files with reads that were mapped to selected taxa, saved in the 
# `output/kraken2_pipeline/kraken2_filtered_reads` directory. The output file name format is
# `"${date_time}"_"${sample_name}"_taxid_"${taxid}"_R"${readnum}".fasta` where
# `taxid` is the number corresponding to NCBI taxonomic ID of a selected species, and 
# `readnum` is read file number (`1` or `2`) where the sequence was found.

# 6. Classifies reads from selected taxa with nucleotide BLAST.

# Input: FASTA files with reads mapped to selected taxa in the 
# # TODO `output/kraken2_pipeline/kraken2_filtered_reads` directory, reference database for BLAST saved as
# `~/common_data/blast_databases/16S_ribosomal_RNA`.

# Output: TXT files with taxonomic classification from BLAST on selected species, saved
# in the `output/kraken2_pipeline/blast_output` directory. The output file name format is
# `"${date_time}"_"${sample_name}"_R"${readnum}"_blast_out.txt`

######


# Input
# FASTQ files with decontaminated reads (either from pipeline #1 or supplied by user). 
# Please adjust Kraken2 and Bracken parameters like kmer length, minimizer 
# length, etc. according to your needs

# Output
# * Kraken2 output files (.kraken2) without minimizers. However, it's possible 
# to add minimizers if needed
# * Kraken2 reports (.k2report) with taxonomic classification
# * FASTQ files with classified sequences (read 1 and read 2 are in 
# separate files), gzip-compressed
# * Bracken output files (.bracken)
# * Bracken reports (.breport)
# * TXT files for Krona generated from Bracken reports (one TXT file per sample)
# * HTML files with Krona plots generated from TXT files (one HTML file per sample)
# * FASTA files with reads that were mapped to Cryptomeria japonica (taxonomy id 3369),
# Glycine soja (taxid 3848), Ipomoea triloba (taxid 35885), and Macadamia 
# integrifolia (taxid 60698).
# * TXT files with taxonomic classification from BLAST on selected species 
# (using FASTA files)


### Get the current date in yyyy-mm-dd format with date -I
### Remove the hyphen with sed (use global search with /g to remove all hyphens)
### Store as a variable date_var
### directory with database and taxonomy
# kraken2_db_dir=data/kraken2_db/k2_large_${kraken2_db_date} 

conda activate kraken2-tools-2.1.3
echo "${start_date_time}"



# 1. Remove host DNA with bowtie2 (done)

# 2. Classify microbiome samples using Kraken2

# We need the database and FASTA file of sequences
# kraken2 --db $DBNAME seqs.fq
# Run with no minimizer data
# If you want to run with minimizers, the output will be 
#  --output ${kraken2_output_dir}/${base_name}.kraken2 instead of 
#  --output ${kraken2_output_dir}/${base_name}_no_minimizer_data.kraken2
# And the --report will be 
#  --report ${kraken2_reports_dir}/${base_name}.k2report
# instead of --report ${kraken2_reports_dir}/${base_name}_no_minimizer_data.k2report
# And you add  --report-minimizer-data 
# to kraken2 command
# And gzip command will be gzip -9 --best ${kraken2_output_dir}/${base_name}.kraken2 
# instead of gzip -9 --best ${kraken2_output_dir}/${base_name}_no_minimizer_data.kraken2
# And tee command will be tee kraken2_stdout.txt
# instead of tee kraken2_stdout_no_minimizer_data.txt
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
echo "Running kraken2"
for FILE in ${bowtie2_decontam_fastq_dir}/*decontam_R1.fastq.gz
do 
 SAMPLE=$(echo ${FILE} | sed "s/_trim_decontam_R1\.fastq\.gz//")
 base_name=$(basename "$SAMPLE" )
 intermediate_date_time=$(date +"%F %H:%M:%S")
 echo "${intermediate_date_time}"
 echo "Running kraken2 on ${SAMPLE}"
 kraken2 --paired \
	--db ${kraken2_db_dir} \
	--threads ${nthreads} \
	--minimum-hit-groups 3 \
	--gzip-compressed \
	${bowtie2_decontam_fastq_dir}/${base_name}_trim_decontam_R1.fastq.gz \
	${bowtie2_decontam_fastq_dir}/${base_name}_trim_decontam_R2.fastq.gz \
	--output ${kraken2_output_dir}/${base_name}_no_minimizer_data.kraken2 \
	--classified-out ${kraken2_classified_reads_dir}/${base_name}_classified_#.fq \
	--report ${kraken2_reports_dir}/${base_name}_no_minimizer_data.k2report \
	--confidence ${confidence_kraken2}
 ls -alh ${kraken2_output_dir}/${base_name}*
 gzip -9 --best ${kraken2_output_dir}/${base_name}_no_minimizer_data.kraken2
 gzip -9 --best ${kraken2_classified_reads_dir}/${base_name}_classified__1.fq
 gzip -9 --best ${kraken2_classified_reads_dir}/${base_name}_classified__2.fq
 ls -alh ${kraken2_output_dir}/${base_name}*
 echo "running rclone on ${SAMPLE}"
 ~/rclone-v1.70.2-linux-amd64/rclone copy -v --create-empty-src-dirs \
 	${kraken2_output_dir}/${base_name}* \
	gdrive:projects/metagenome/supercomputer/kraken2_pipeline/kraken2_output/
 rm ${kraken2_output_dir}/${base_name}_no_minimizer_data.kraken2.gz
 rm ${kraken2_classified_reads_dir}/${base_name}_classified__*.fq.gz;
done  2>&1 |tee kraken2_stdout_no_minimizer_data.txt   
### --db
### --threads
### -minimum-hit groups
### --gzip-compressed \
### ${bowtie2_decontam_fastq_dir}/${base_name}_decontam_R1.fastq.gz \
### ${bowtie2_decontam_fastq_dir}/${base_name}_decontam_R2.fastq.gz \
### --output ${kraken2_output_dir}/${base_name}.kraken2 \
### --classified-out ${kraken2_classified_reads_dir}/${base_name}_classified_#.fq \

### --classified-out will print classified sequences to filename
### --report-minimizer-data adds two columns to the report file (#4 and #5) , which
### represent the number of minimizers found to be associated with a taxon
### in the read sequences (#4), and the estimate of the number of distinct
### minimizers associated with a taxon in the read sequence data (#5)

# 3. Run bracken for abundance estimation of microbiome samples
# Run with no minimizer data
# If you want to run with minimizers, the FOR command will be 
# for FILE in ${kraken2_reports_dir}/*.k2report
# instead of ${kraken2_reports_dir}/*_no_minimizer_data.k2report
# And sed will be SAMPLE=$(echo ${FILE} | sed "s/\.k2report//"|sed "s///")
# instead of SAMPLE=$(echo ${FILE} | sed "s/\.k2report//"|sed "s/_no_minimizer_data//")
# And input will be -i ${kraken2_reports_dir}/${base_name}.k2report 
# instead of -i ${kraken2_reports_dir}/${base_name}_no_minimizer_data.k2report
# And output will be -o ${bracken_output_dir}/${base_name}.bracken 
# instead of -o ${bracken_output_dir}/${base_name}_no_minimizer_data.bracken 
# And the report will be -w ${bracken_reports_dir}/${base_name}.breport
# instead of -w ${bracken_reports_dir}/${base_name}_no_minimizer_data.breport
# And gzip will be gzip -9 --best ${kraken2_reports_dir}/${base_name}.k2report;
# instead of gzip -9 --best ${kraken2_reports_dir}/${base_name}_no_minimizer_data.k2report
# And tee command will be tee bracken_stdout.txt
# instead of tee bracken_stdout_no_minimizer_data.txt
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
echo "Running bracken"
for FILE in ${kraken2_reports_dir}/*_no_minimizer_data.k2report
do 
  SAMPLE=$(echo ${FILE} | sed "s/\.k2report//"| sed "s/_no_minimizer_data//")
  base_name=$(basename "$SAMPLE" )
  intermediate_date_time=$(date +"%F %H:%M:%S")
  echo "${intermediate_date_time}"
  echo "Running bracken on ${SAMPLE}"
  bracken -d ${kraken2_db_dir} \
	-i ${kraken2_reports_dir}/${base_name}_no_minimizer_data.k2report \
	-r ${read_len} \
	-l S \
	-t ${bracken_threshold} \
	-o ${bracken_output_dir}/${base_name}_no_minimizer_data.bracken \
	-w ${bracken_reports_dir}/${base_name}_no_minimizer_data.breport
  gzip -9 --best ${kraken2_reports_dir}/${base_name}_no_minimizer_data.k2report;
done 2>&1 |tee bracken_stdout_no_minimizer_data.txt
### -d
### -i
### -r is read length
### -l is classification level
### -t is threshold: specifies the minimum number of reads required for 
### a classification at the specified rank. Any classifications with less than 
### the specified threshold will not receive additional reads from higher taxonomy 
### levels when distributing reads for abundance estimation.
### -o
### -w
#FILE: testdir/20230104_f1.txt
#SAMPLE: testdir/f1
#base_name: f1

#FILE:testdir/20230104_f2.txt
#SAMPLE:testdir/f2
#base_name:f2

#FILE: testdir/20230104_file.txt
#SAMPLE:testdir/file
#base_name:file

# 4. Generate Krona plots
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
echo "Generating Krona plots"
for FILE in ${bracken_reports_dir}/*.breport
do 
  SAMPLE=$(echo ${FILE} | sed "s/\.breport//" )
  base_name=$(basename "$SAMPLE" )
  kreport2krona.py -r ${bracken_reports_dir}/${base_name}.breport \
    -o ${bracken_krona_txt_dir}/${base_name}.b.krona.txt \
	  --no-intermediate-ranks
  ktImportText ${bracken_krona_txt_dir}/${base_name}.b.krona.txt \
	-o ${krona_html_dir}/${base_name}.krona.html;
done
### For kreport2krona.py:
### -r
### -o
### --no-intermediate-ranks
### For ktImportText:
### input is
### -o

# Tidy up the bracken report
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"
for FILE in "${bracken_krona_txt_dir}"/"${date_time}"_*_trim_no_minimizer_data.b.krona.txt
 do
 # Extract the sample name
 SAMPLE=$(basename "${FILE}" | sed "s/\.b.krona.txt//"  | sed "s/_trim_no_minimizer_data//")
 # Process the file using awk
 awk -v sample_name="${SAMPLE}" '
	BEGIN {FS=OFS="\t";
		print "Sample", "Abundance", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"
	}
	{
		abundance = $1
		kingdom = ""
		phylum = ""
		class = ""
		order = ""
		family = ""
		genus = ""
		species = ""

		for (i = 2; i <= NF; i++) {
			if ($i ~ /^k__/) kingdom = $i
			else if ($i ~ /^p__/) phylum = $i
			else if ($i ~ /^c__/) class = $i
			else if ($i ~ /^o__/) order = $i
			else if ($i ~ /^f__/) family = $i
			else if ($i ~ /^g__/) genus = $i
			else if ($i ~ /^s__/) species = $i
		}

		print sample_name, abundance,  kingdom, phylum, class, order, family, genus, species
	}
 ' "${bracken_krona_txt_dir}"/"${date_time}"_"${SAMPLE}"_trim_no_minimizer_data.b.krona.txt > \
	"${bracken_krona_txt_dir}"/"${date_time}"_"${SAMPLE}"_trim_no_minimizer_data.tsv
done

awk 'FNR>1' "${bracken_krona_txt_dir}"/"${date_time}"_*_trim_no_minimizer_data.tsv | \
	awk 'BEGIN {FS=OFS="\t";
		print "Sample", "Abundance", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"}
		{print}' > output/rtables/"${date_time}"_combined_report.tsv
intermediate_date_time=$(date +"%F %H:%M:%S")
echo "${intermediate_date_time}"

# Analyse specific sequences
# Cryptomeria japonica: taxid 3369
# Glycine soja: taxid 3848
# Ipomoea triloba: taxid 35885 
# Macadamia integrifolia: taxid 60698

# # 1. Extract reads that were mapped to selected taxid
# for FASTQ_FILE in "${kraken2_classified_reads_dir}"/"${date_time}"_*_classified__1.fq.gz
# do 
#   SAMPLE=$(echo "${FASTQ_FILE}" | sed "s/_classified__1\.fq\.gz//" )
#   base_name=$(basename "$SAMPLE" )
#   for taxid in 3369 3848 35885
# #   for taxid in 338188 10181 823
#   do
# 	for readnum in 1 2
# 	do
#       zcat "${kraken2_classified_reads_dir}"/"${date_time}"_"${base_name}"_classified__"${readnum}".fq.gz |\
# 	    grep --no-group-separator -A 1 -E "@.*kraken:taxid\|$taxid$" | sed "s/@/>/" > \
# 	    "${kraken2_filtered_reads_dir}"/"${date_time}"_"${base_name}"_taxid_"${taxid}"_R"${readnum}".fasta
#     #   gzip -9 --best \
# 	#     "${kraken2_filtered_reads_dir}"/"${date_time}"_"${base_name}"_taxid_"${taxid}"_R"${readnum}".fasta;
# 	done
#   done
# done

# conda activate qc-tools
# mkdir -p ~/common_data/blast_databases
# mkdir -p output/kraken2_pipeline/blast_output
# blast_output_dir=output/kraken2_pipeline/blast_output
# blast_database_path=~/common_data/blast_databases/16S_ribosomal_RNA
# # To download a preformatted NCBI BLAST database, run update_blastdb.pl program followed by any relevant
# # options and the name(s) of the BLAST databases to download
# # update_blastdb.pl --decompress nr [*]
# cd ~/common_data/blast_databases
# update_blastdb.pl --decompress 16S_ribosomal_RNA 
# cd "${project_home_dir}"

# #  *** Formatting options
# #  -outfmt <String>
# #    alignment view options:
# #      0 = Pairwise,
# #      1 = Query-anchored showing identities,
# #      2 = Query-anchored no identities,
# #      3 = Flat query-anchored showing identities,
# #      4 = Flat query-anchored no identities,
# #      5 = BLAST XML,
# #      6 = Tabular,
# #      7 = Tabular with comment lines,
# #      8 = Seqalign (Text ASN.1),
# #      9 = Seqalign (Binary ASN.1),
# #     10 = Comma-separated values,
# #     11 = BLAST archive (ASN.1),
# #     12 = Seqalign (JSON),
# #     13 = Multiple-file BLAST JSON,
# #     14 = Multiple-file BLAST XML2,
# #     15 = Single-file BLAST JSON,
# #     16 = Single-file BLAST XML2,
# #     17 = Sequence Alignment/Map (SAM),
# #     18 = Organism Report


# #    The supported format specifiers for options 6, 7 and 10 are:
# #             qseqid means Query Seq-id
# #                qgi means Query GI
# #            qaccver means Query accession.version
# #               qlen means Query sequence length
# #             sseqid means Subject Seq-id
# #          sallseqid means All subject Seq-id(s), separated by a ';'
# #                sgi means Subject GI
# #             sallgi means All subject GIs
# #            saccver means Subject accession.version
# #            sallacc means All subject accessions
# #               slen means Subject sequence length
# #               qseq means Aligned part of query sequence
# #               sseq means Aligned part of subject sequence
# #              score means Raw score
# #             length means Alignment length
# #             nident means Number of identical matches
# #           positive means Number of positive-scoring matches
# #               gaps means Total number of gaps
# #               ppos means Percentage of positive-scoring matches
# #             frames means Query and subject frames separated by a '/'
# #             qframe means Query frame
# #             sframe means Subject frame
# #               btop means Blast traceback operations (BTOP)
# #             staxid means Subject Taxonomy ID
# #           ssciname means Subject Scientific Name
# #           scomname means Subject Common Name
# #         sblastname means Subject Blast Name
# #          sskingdom means Subject Super Kingdom
# #            staxids means unique Subject Taxonomy ID(s), separated by a ';'
# #                          (in numerical order)
# #          sscinames means unique Subject Scientific Name(s), separated by a ';'
# #          scomnames means unique Subject Common Name(s), separated by a ';'
# #         sblastnames means unique Subject Blast Name(s), separated by a ';'
# #                          (in alphabetical order)
# #         sskingdoms means unique Subject Super Kingdom(s), separated by a ';'
# #                          (in alphabetical order)
# #             stitle means Subject Title
# #         salltitles means All Subject Title(s), separated by a '<>'
# #              qcovs means Query Coverage Per Subject
# #            qcovhsp means Query Coverage Per HSP
# #             qcovus means Query Coverage Per Unique Subject (blastn only)


# # qacc means Query accession
# # sacc means Subject accession
# # pident means Percentage of identical matches
# # mismatch means Number of mismatches
# # gapopen means Number of gap openings
# # qstart means Start of alignment in query
# # qend means End of alignment in query
# # sstart means Start of alignment in subject
# # send means End of alignment in subject
# # evalue means Expect value
# # bitscore means Bit score
# # sstrand means Subject Strand
# # btop means Blast traceback operations (BTOP) (Similar to CIGAR format in SAM, but more flexible)

# # Don't use btop because the output file size will be too big
# for FASTA_FILE in "${kraken2_filtered_reads_dir}"/"${date_time}"_*_taxid_*_R1.fasta
# do 
#   SAMPLE=$(echo "${FASTA_FILE}" | sed "s/_R1\.fasta//")
#   base_name=$(basename "$SAMPLE" )
#   for readnum in 1 2
#     do
# 	  echo "Running BLAST on ${SAMPLE}_R${readnum}"
# 	  blastn -db "${blast_database_path}" \
# 	    -query "${kraken2_filtered_reads_dir}"/"${date_time}"_"${base_name}"_R"${readnum}".fasta \
# 		-out "${blast_output_dir}"/"${date_time}"_"${base_name}"_R"${readnum}"_blast_out.txt \
# 		-outfmt "7 qacc sacc pident mismatch gapopen qstart qend sstart send evalue bitscore";
# 	done
# done

# blastn -db "${blast_database_path}" \
#   -query "${kraken2_filtered_reads_dir}"/Y51b_taxid_823_R1.fasta \
#   -out "${blast_output_dir}"/Y51b_blast_823_16s_local-nobtop.txt \
#   -outfmt "7 qacc sacc pident mismatch gapopen qstart qend sstart send evalue bitscore"

# # Which matches to select? There's many, so need to filter out


# end_date_time=$(date +"%F %H:%M:%S")
# echo "${end_date_time}"