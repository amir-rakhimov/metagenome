#!/usr/bin/env bash
date_var=$(date -I|sed 's/-//g')
kmer_len=35 
read_len=150
fastq_dir=data/fastq/yasuda-fastq
kraken2_db_dir=data/kraken2_db/k2_standard_20231225
kraken2_reports_dir=output/kraken2_reports
kraken2_output_dir=output/kraken2_output
bracken_reports_dir=output/bracken_reports
bracken_output_dir=output/bracken_output
bracken_krona_txt_dir=output/bracken_krona_txt
krona_html_dir=output/krona_html
source ~/miniconda3/etc/profile.d/conda.sh
conda activate kraken2-tools-2.1.3
## Download databases
#kraken2-build --standard --db ${kraken2_db_dir}
## Build the database: uses taxonomy and library
#kraken2-build --build --threads 10 --db ${kraken2_db_dir} --kmer-len 35 --minimizer-len 31
bracken-build -d ${kraken2_db_dir} -t 35 -k ${kmer_len} -l ${read_len}
# 2. Classify microbiome samples using Kraken
for FILE in ${fastq_dir}/*L3_1.fq.gz; 
do 
  SAMPLE=$(echo ${FILE} | sed "s/_L3_1\.fq\.gz//")
  base_name=$(basename "$SAMPLE" )
  kraken2 --paired \
	--db ${kraken2_db_dir} \
	--threads 10 \
	--minimum-hit-groups 3 \
	--gzip-compressed \
	${fastq_dir}/${base_name}_L3_1.fq.gz \
	${fastq_dir}/${base_name}_L3_2.fq.gz \
	--output ${kraken2_output_dir}/${date_var}_${base_name}.kraken2 \
	--classified-out ${kraken2_output_dir}/${date_var}_${base_name}_classified_#.fq \
	--report ${kraken2_reports_dir}/${date_var}_${base_name}.k2report \
	--report-minimizer-data 
done
# 3. Run bracken for abundance estimation of microbiome samples
for FILE in ${kraken2_reports_dir}/*.k2report; 
do 
  SAMPLE=$(echo ${FILE} | sed "s/\.k2report//")
  base_name=$(basename "$SAMPLE" )
  bracken -d ${kraken2_db_dir} \
	-i ${kraken2_reports_dir}/${date_var}_${base_name}.k2report \
	-r ${read_len} \
	-l S \
	-t 10 \
	-o ${bracken_output_dir}/${date_var}_${base_name}.bracken \
	-w ${bracken_reports_dir}/${date_var}_${base_name}.breport;
done
# Generate Krona plots
kreport2krona.py -r ${bracken_reports_dir}/${date_var}_2D10_wms.breport \
	-o ${bracken_krona_txt_dir}/${date_var}_2D10_wms.b.krona.txt \
	--no-intermediate-ranks
ktImportText ${bracken_krona_txt_dir}/${date_var}_2D10_wms.b.krona.txt \
	-o ${krona_html_dir}/${date_var}_2D10_wms.krona.html
kreport2krona.py -r ${bracken_reports_dir}/${date_var}_H4_wms.breport \
	-o ${bracken_krona_txt_dir}/${date_var}_H4_wms.b.krona.txt \
	--no-intermediate-ranks
ktImportText ${bracken_krona_txt_dir}/${date_var}_H4_wms.b.krona.txt \
	-o ${krona_html_dir}/${date_var}_H4_wms.krona.html
kreport2krona.py -r ${bracken_reports_dir}/${date_var}_Y51b_wms.breport \
	-o ${bracken_krona_txt_dir}/${date_var}_Y51b_wms.b.krona.txt \
	--no-intermediate-ranks
ktImportText ${bracken_krona_txt_dir}/${date_var}_Y51b_wms.b.krona.txt \
	-o ${krona_html_dir}/${date_var}_Y51b_wms.krona.html
