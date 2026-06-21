#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob
RUN_ID=$(date +%Y%m%d_%H%M%S)
# Load config
source config/config.sh

mkdir -p data/4-mag/dbcan_db
mkdir -p data/4-mag/checkm2_db
mkdir -p "${silva_db_dir}"
mkdir -p "${blastdb_path}"
mkdir -p "${blastdb_path}"
mkdir -p "${cazy_db_taxonomy_dir}"
mkdir -p "${cazy_db_fasta_dir}"

mkdir -p "${megahit_output_dir}"
mkdir -p "${megahit_logs_dir}"
mkdir -p "${megahit_aligned_reads_dir}"
mkdir -p "${bbwrap_refs_dir}"
mkdir -p "${contig_coverages_dir}"
mkdir -p "${metaquast_megahit_output_dir}"

mkdir -p "${bam_contig_depths_dir}"
mkdir -p "${metabat2_output_dir}"
mkdir -p "${metaquast_metabat2_output_dir}"
mkdir -p "${metabat2_logs_dir}"

mkdir -p "${checkm2_output_dir}"
mkdir -p "${drep_output_dir}"
mkdir -p "${mag_qc_logs_dir}"

mkdir -p "${gtdbtk_output_dir}"
mkdir -p "${comparem_output_dir}"/aai_output_sa95perc
mkdir -p "${comparem_output_dir}"/aai_output_sa99perc


mkdir -p "${prokka_output_dir}"
mkdir -p "${prokka_logs_dir}"
mkdir -p "${prokka_gtf_dir}"
mkdir -p "${prokka_contig_gene_maps_dir}"
mkdir -p "${mmseqs_easy_cluster_output_dir}"
mkdir -p "${dbcan_logs_dir}"
mkdir -p "${kofam_scan_output_dir}"
mkdir -p "${annotation_fasta_dir}"

mkdir -p "${picard_output_dir}"
mkdir -p "${picard_logs_dir}"
mkdir -p "${htseq_gene_count_dir}"
mkdir -p "${htseq_logs_dir}"
mkdir -p "${gene_lengths_dir}"
mkdir -p "${tpm_files_dir}"

mkdir -p "${coverm_logs_dir}"
mkdir -p "${coverm_drep_95_dir}"
mkdir -p "${coverm_drep_99_dir}"
mkdir -p "${coverm_metabat2_dir}"
mkdir -p "${coverm_megahit_dir}"

mkdir -p "${seqkit_output_dir}"
mkdir -p "${seqkit_megahit_dir}"
mkdir -p "${seqkit_drep_95_dir}"
mkdir -p "${seqkit_metabat2_dir}"

mkdir -p "${blastn_output_dir}"

# mkdir -p "${mag_OUTDIR}"/09-cazyme_blast/fasta ? Not needed?
mkdir -p "${blastp_output_dir}"
mkdir -p "${blastp_logs_dir}"

mkdir -p "${metathermo_out_dir}"

job_id1=$(sbatch --parsable --export=ALL code/slurm-scripts/megahit.slurm)
job_id2=$(sbatch --parsable --dependency=afterok:$job_id1 --export=ALL code/slurm-scripts/metabat2.slurm)
job_id3=$(sbatch --parsable --dependency=afterok:$job_id2 --export=ALL code/slurm-scripts/mag-qc.slurm)
job_id4=$(sbatch --parsable --dependency=afterok:$job_id3 --export=ALL code/slurm-scripts/taxonomic-classification.slurm)
job_id5=$(sbatch --parsable --dependency=afterok:$job_id4 --export=ALL code/slurm-scripts/genome-annotation.slurm)
job_id6=$(sbatch --parsable --dependency=afterok:$job_id5 --export=ALL code/slurm-scripts/quantify-genes.slurm)
job_id7=$(sbatch --parsable --dependency=afterok:$job_id6 --export=ALL code/slurm-scripts/mag-stats.slurm)
job_id8=$(sbatch --parsable --dependency=afterok:$job_id7 --export=ALL code/slurm-scripts/blast-contigs.slurm)
job_id9=$(sbatch --parsable --dependency=afterok:$job_id8 --export=ALL code/slurm-scripts/search-cazymes-blast.slurm)
job_id10=$(sbatch --parsable --dependency=afterok:$job_id9 --export=ALL code/slurm-scripts/meta-thermo.slurm)

echo "Pipeline submitted"
echo " step 1: ${job_id1}"
echo " step 2: ${job_id2}"
echo " step 2: ${job_id3}"
echo " step 2: ${job_id4}"
echo " step 2: ${job_id5}"
echo " step 2: ${job_id6}"
echo " step 2: ${job_id7}"
echo " step 2: ${job_id8}"
echo " step 2: ${job_id9}"
echo " step 2: ${job_id10}"
