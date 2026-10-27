#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob
source config/bash/config.sh
source ../utils/Bash/prefetch_download_and_rename.sh 
author_sra_files_dir="${sra_files_dir}"/tanes
author_fastq_dir="${fastq_dir}"/tanes
author_metadata_dir=data/metadata/tanes
mkdir -p "${author_sra_files_dir}"
mkdir -p "${author_metadata_dir}"
mkdir -p "${author_fastq_dir}"

mamba create -yqn metagenome-reproduce-tools -c conda-forge -c bioconda \
    bioconda::fastp  bioconda::trimmomatic comebin \
    bioconda::panaroo bioconda::diamond \
    mafft  bioconda::fastani

mamba create -yqn bakta-env bioconda::bakta==1.11.4
# Install bakta separately  bioconda::kraken2==2.17.1    krakentools bracken
# Install metaspades and bwa
# wget https://www.cell.com/cms/10.1016/j.celrep.2026.117265/attachment/35992bf2-2290-49e7-81e9-1f3d8b2933b7/mmc3.xlsx \
#     -O "${author_metadata_dir}"/mmc3.xlsx

archaea_mags=("c.013.380.331" "c.013.519.256" "c.107.560.558" "c.107.599.878"
                "c.107.603.871" "c.107.779.295" "NMR.107779295.1"
                "NMR.13357528.1"  "NMR.57.804.626.1" "NMR.603001892.1")

# First, download the SRA Run table from NCBI
# This will print SRA, BioSample, SampleName
# We also filter to keep NMR samples
awk -F',' 'BEGIN {OFS=","} 
    NR ==1 {print "SRA", "BioSample", "SampleName"} 
    NR>1 && $40 ~ /Heterocephalus/{print $1,$7,$55}'  \
        "${author_metadata_dir}"/SraRunTable.csv > \
        "${author_metadata_dir}"/tanes-sra-runinfo-nmr.csv

# Filter to keep only archaea samples
if [[ ! -f "${author_metadata_dir}"/tanes-sra-runinfo-nmr-archaea.csv ]]; then
    awk -F',' -v mags="$(IFS="|"; echo "${archaea_mags[*]}")" '
        BEGIN {OFS="\t";
            split(mags, arr, "|")
            for (i in arr) ids[arr[i]]
        } 
            NR==1 {print $0}; $3 in ids' \
            "${author_metadata_dir}"/tanes-sra-runinfo-nmr.csv > \
            "${author_metadata_dir}"/tanes-sra-runinfo-nmr-archaea.csv
fi

# Next, create metadata (age, sex, etc)
# We filter to keep NMR samples here, too
awk -F',' 'BEGIN {OFS="\t"} 
    NR ==1 {print "sample_id", "host", "host_age", "breeding_status", "colony_id", "host_sex" } 
    NR>1 && $40 ~ /Heterocephalus/ {print $55, $40, $38, $10, $14, $41}'  \
        "${author_metadata_dir}"/SraRunTable.csv > \
        "${author_metadata_dir}"/tanes-metadata.tsv


# Finally, just sample IDs that will be used
awk -F',' 'BEGIN {OFS="\t"} 
    NR ==1 {print "sample_id"}
	NR >1 {print $3}' "${author_metadata_dir}"/tanes-sra-runinfo-nmr-archaea.csv >\
        "${SAMPLES_TSV}"

# # For this one, get the mmc3.tsv file from the paper
# # Next, search selected sample ids in the tsv file (archaea_mags)
# # NMR.13357528 is missing, but NMR.13357528.1 was found
# if [[ -f "${author_metadata_dir}"/tanes-mag.tsv ]]; then
#     awk -F'\t' -v mags="$(IFS="|"; echo "${archaea_mags[*]}")" '
#         BEGIN {OFS=',';
#             split(mags, arr, "|")
#             for (i in arr) ids[arr[i]]
#         } 
#             FNR==1 {print $0}; $2 in ids' "${author_metadata_dir}"/mmc3.tsv > "${author_metadata_dir}"/tanes-mag.tsv
# fi

# # Get SRA IDs for selected samples based on the Bioproject ID and Biosample IDs
# if [[ -f "${author_metadata_dir}"/tanes-sra-runinfo-filtered.csv ]]; then
#     bioproject_id=PRJNA1247447
#     esearch -db sra -query "${bioproject_id}" |efetch -format runinfo > "${author_metadata_dir}"/tanes-sra-runinfo.csv
#     awk -F',' -v biosamples="$(cut -f1 "${author_metadata_dir}"/tanes-mag.tsv | tail +2 | tr '\n' '|')" '
#         BEGIN {OFS=",";
#             split(biosamples, arr, "|")
#             for (i in arr) ids[arr[i]]
#         } 
#             FNR==1 {print $0}; $26 in ids'  "${author_metadata_dir}"/tanes-sra-runinfo.csv |\
#             awk -F',' 'BEGIN{OFS=",";} NR==1 {print "SRA,BioSample,SampleName"} NR>1 {print $1,$26,$12}'  > \
#                 "${author_metadata_dir}"/tanes-sra-runinfo-filtered.csv
# fi

prefetch_download_and_rename "${author_sra_files_dir}" "${author_fastq_dir}" \
    "${author_metadata_dir}"/tanes-sra-runinfo-filtered.csv \
    "${sratoolkit_scripts_dir}" 

