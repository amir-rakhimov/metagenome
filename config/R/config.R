project.root <- here::here()

kraken2.run_id<- "20240409_17_32_40"
mag.run_id <- "20250303_17_54_26"


results.kraken2<- function(...) file.path(project.root, "results","1-kraken2", 
                                  kraken2.run_id, ...)
kraken2.analysis.out <- function(...) file.path(project.root, "analysis", 
                                                "1-kraken2", kraken2.run_id, ...)
results.mag <- function(...) file.path(project.root, "results","2-mag", 
                                       mag.run_id, ...)
mag.analysis.out <- function(...)  file.path(project.root, "analysis", "2-mag",
                                             mag.run_id, ...)
maaslin2.out<- file.path(project.root, "analysis", "3-maaslin2")

kraken2.manuscript.figures <- kraken2.analysis.out( "manuscript-figures")
mag.manuscript.figures <- mag.analysis.out( "manuscript-figures")

# Directories with output from R
community.composition.tables <- kraken2.analysis.out("01-community_composition",
                                                            "tables")
community.composition.rdafiles <-kraken2.analysis.out("01-community_composition",
                                                              "rdafiles")
community.composition.figures <-kraken2.analysis.out("01-community_composition",
                                                             "figures")
diffabund.tables <-kraken2.analysis.out("02-differential_abundance", "tables")
diffabund.rdafiles <-kraken2.analysis.out("02-differential_abundance", "rdafiles") 
diffabund.figures <-kraken2.analysis.out("02-differential_abundance", "figures")  



#' Classification statistics (symlink)
kraken2.unclassified.stats.fname <-results.kraken2(paste("unclassified_reads.tsv",
                sep="_"))
#' Directory with metadata:
metadatadir<-paste0("../amplicon_nmr/data/metadata/shared") 
# Kraken2 report (symlink)
combined.report.filename <- results.kraken2("combined_report.tsv")

# Directory with necessary functions
util.functions.r <- "../utils/R"
# Phyloseq object
ps.q.raw.fname <- file.path(community.composition.rdafiles, 
                            paste("phyloseq-kraken2","ps.q.raw.rds",sep="-"))
ps.q.rel.raw.fname <- file.path(community.composition.rdafiles,
                                paste("phyloseq-kraken2","ps.q.rel.raw.rds",sep="-"))
ps.q.raw.otu_table.fname <- file.path(community.composition.rdafiles,
                                      paste("phyloseq-kraken2","ps.q.raw.otu_table.rds",sep="-"))
ps.q.raw.tax_table.fname <- file.path(community.composition.rdafiles,
                                      paste("phyloseq-kraken2","ps.q.raw.tax_table.rds",sep="-"))
rare.num_samples <- 100

# For saving phyloseq objects as dataframes
ps.q.agg.asv.fname.no_ext <- paste("phyloseq-kraken2","OTU", 
                                  "table",sep = "-")
ps.q.agg.species.fname.no_ext <- paste("phyloseq-kraken2","Species", 
                                       "table",sep = "-")
ps.q.agg.genus.fname.no_ext <-paste("phyloseq-kraken2","Genus",
                                    "table",sep = "-")
ps.q.agg.family.fname.no_ext <- paste("phyloseq-kraken2","Family",
                                    "table",sep = "-")
ps.q.agg.phylum.fname.no_ext <- paste("phyloseq-kraken2","Phylum",
                                      "table",sep = "-")

############# TO BE REMOVED #################
resolve_wsl_symlink <- function(wsl_unc_path) {
  # Strip UNC prefix to get Linux path
  linux_path <- gsub("\\\\", "/", wsl_unc_path)
  linux_path <- sub("//wsl.localhost/Ubuntu-22.04", "", linux_path)
  
  # Resolve symlink inside WSL
  resolved_linux <- system(paste("wsl readlink -f", shQuote(linux_path)), intern = TRUE)
  
  # Convert back to UNC
  gsub("^/", "//wsl.localhost/Ubuntu-22.04/", resolved_linux)
}

combined.report.filename <- resolve_wsl_symlink(combined.report.filename)
kraken2.unclassified.stats.fname <- resolve_wsl_symlink(kraken2.unclassified.stats.fname)
######################

#' Import metadata:
custom.md.path<-file.path(metadatadir,"custom.md.ages.rds")

image.formats<-c("png","tiff")



# MAG
mag.tables <- mag.analysis.out("tables")
mag.rdafiles <- mag.analysis.out("rdafiles")
mag.figures <- mag.analysis.out("figures")
#' First, statistics from seqkit.
#' Seqkit stats on MEGAHIT data
seqkit.megahit.date_time<-"20250712_18_07_37"
#' Seqkit stats on dereplicated MAGs
seqkit.drep.date_time<-"20250712_18_07_37"
#' Seqkit stats on MetaBAT2 data
seqkit.metabat2.date_time<-"20250712_18_07_37"
#' Next, taxonomy tables.
#' BLASTN output
blastn.date_time<-"20250701_05_24_10"
#' GTDB-tk output
gtdbtk.date_time<-"20250705_19_24_02"
# MetaBAT2 output to clean up GTDB-tk genome names:
metabat2.date_time<-"20250612_13_37_47"
#' MAG and contig coverages from CoverM 
coverm.drep.date_time<-"20250803_05_07_53"
coverm.metabat2.date_time<-"20250803_05_07_53"
coverm.megahit.date_time<-"20250803_05_07_53"
#' Directories with input files
seqkit.output.dir <- results.mag("07-mag_stats/seqkit")
blastn.output.dir <- results.mag("08-contig_blast/",blastn.date_time)
gtdbtk.output.dir <- results.mag("04-gtdbtk")
coverm.output.dir <- results.mag("07-mag_stats/coverm")
checkm2.output.dir <- results.mag("03-qc/checkm2")
# Specify input file names
seqkit.megahit.fname <- file.path(seqkit.output.dir, "megahit",
                                paste(seqkit.megahit.date_time,
                                      "megahit_combined_stats.tsv",sep="_"))
seqkit.drep.fname<-file.path(seqkit.output.dir, "drep_95",
                             paste(seqkit.drep.date_time,
                                   "sa_95perc_combined_stats.tsv",sep="_"))
seqkit.metabat2.fname<-file.path(seqkit.output.dir,"metabat2",
                                 paste(seqkit.metabat2.date_time,
                                       "metabat2_combined_stats.tsv",sep="_"))
blastn.taxonomy.fname<-file.path(blastn.output.dir,
                                 paste(blastn.date_time,
                                       "blastn_all_samples.tsv",sep="_"))
gtdbtk.taxonomy.fname<-file.path(gtdbtk.output.dir,
                                 paste(gtdbtk.date_time,
                                       "sa_95perc_classification_combined.tsv",
                                       sep="_"))
coverm.drep.fname<-file.path(coverm.output.dir, "drep_95",
                             paste(coverm.drep.date_time,
                                   "all_samples_sa_95perc_coverm.tsv",
                                   sep="_"))
coverm.metabat2.fname<-file.path(coverm.output.dir, #"metabat2", # outside metabat2 folder???
                                 paste(coverm.metabat2.date_time,
                                       "all_samples_bins_coverm.tsv",
                                       sep="_"))
coverm.megahit.fname<-file.path(coverm.output.dir, # "megahit", # outside the megahit folder??
                                paste(coverm.megahit.date_time,
                                      "all_samples_megahit_coverm.tsv",
                                      sep="_"))
drep.clusters.fname <- results.mag("03-qc/drep/20250627_19_56_42_sa_95perc/data_tables/Cdb.csv") 

drep.winning_genomes.fname<-results.mag("03-qc/drep/20250627_19_56_42_sa_95perc/data_tables/Wdb.csv") 

checkm2.high_quality_mags.fname <- file.path(checkm2.output.dir, "20250619_05_47_09_high_quality_mags.tsv") 
checkm2.all_mags.fname <- file.path(checkm2.output.dir, "20250619_05_47_09_quality_reports_merged.tsv") 

#' Gene detection data:
prokka.date_time <- "20250626_22_11_43"
#' MMseqs clusters (mapping between all proteins and representative sequences):
mmseqs.date_time <- "20250626_22_11_43"
#' Htseq gene counts:
htseq.date_time <- "20250710_19_36_58"
# KofamScan gene annotation:
kofamscan.date_time <- "20250718_09_21_41"
# dbCAN results (on MMseqs proteins):
dbcan.date_time <- "20250626_22_11_43"
#' Directories with datasets:
prokka.gtf.dir <- results.mag("05-annotation/prokka/gtf") 
mmseqs.output.dir <- results.mag("05-annotation/mmseqs2/easy_cluster")
# HTseq gene counts from PROKKA (not normalised):
htseq.gene_count.dir <- results.mag("06-gene_quantification/htseq") 
kofamscan.output.dir <- results.mag("05-annotation/kofam_scan") 
dbcan.output.dir <- results.mag("05-annotation/dbcan",
                                paste(dbcan.date_time,"nr_dbcan_proteins",
                                      sep = "_"))
#' Specify input file names:
prokka.contig.gene.map.fname <- file.path(prokka.gtf.dir,
                                        paste(prokka.date_time,
                                              "all_contig_gene_maps.tsv",
                                              sep="_"))
mmseqs.clusters.fname <- file.path(mmseqs.output.dir,
                                 paste(mmseqs.date_time,
                                       "prokka_nr_prot_cluster.tsv",sep = "_"))
kofamscan.fname <- file.path(kofamscan.output.dir,
                           paste(kofamscan.date_time,
                                 "nr_prot_kofam_scan_top_hits.txt",sep="_"))
dbcan.fname <- file.path(dbcan.output.dir,
                       paste(dbcan.date_time,
                             "nr_dbcan_overview.txt",sep = "_"))
htseq.gene_counts.fname <- file.path(htseq.gene_count.dir,
                                   paste(htseq.date_time,"all_prokka_counts.tsv",
                                         sep="_"))
# BLASTP on cazymes
cazymes.blast.df.fname <- results.mag("09-cazyme_blast/blastp_results/20250802_12_29_03_blastp_cazy_db_cazymes_top_hits_with_taxa.tsv") 

imgmc.ko.fname <- "./data/iMGMC_map_functionality/iMGMC-map-GeneID-KeggKO.tab"
imgmc.cazy.fname <- "./data/iMGMC_map_functionality/iMGMC-map-GeneID-CAZy.tab"

# Names of files that we create in this analysis
seqkit.all.df.fname.rds <- file.path(mag.rdafiles,"seqkit.all.df.rds")
seqkit.all.df.fname.tsv <- file.path(mag.tables,"seqkit.all.df.tsv")

gtdbtk.coverm.drep.fname.rds <- file.path(mag.rdafiles,"gtdbtk.coverm.drep.rds")
gtdbtk.coverm.drep.fname.tsv <- file.path(mag.tables,"gtdbtk.coverm.drep.tsv")

gtdbtk.coverm.metabat2.fname.rds <- file.path(mag.rdafiles,"gtdbtk.coverm.metabat2.rds")
gtdbtk.coverm.metabat2.fname.tsv <- file.path(mag.tables,"gtdbtk.coverm.metabat2.tsv")

blastn.coverm.megahit.fname.rds <- file.path(mag.rdafiles,"blastn.coverm.megahit.rds")
blastn.coverm.megahit.fname.tsv <- file.path(mag.tables,"blastn.coverm.megahit.tsv")

gtdbtk.taxonomy.clean.fname.rds <- file.path(mag.rdafiles, "gtdbtk.taxonomy.clean.rds")
gtdbtk.taxonomy.clean.fname.tsv <- file.path(mag.tables, "gtdbtk.taxonomy.clean.tsv")

gene.annotation.df.fname.rds <- file.path(mag.rdafiles, "gene.annotation.df.rds")
gene.annotation.df.fname.tsv <- file.path(mag.tables, "gene.annotation.df.tsv")

dbcan.df.fname.rds <- file.path(mag.rdafiles, "dbcan.df.rds")
dbcan.df.fname.tsv <- file.path(mag.tables, "dbcan.df.tsv")

kofamscan.df.fname.rds <- file.path(mag.rdafiles, "kofamscan.df.rds")
kofamscan.df.fname.tsv <- file.path(mag.tables, "kofamscan.df.tsv")

representative.gene_lengths.fname.rds <- file.path(mag.rdafiles, "representative.gene_lengths.rds")
representative.gene_lengths.fname.tsv <- file.path(mag.tables, "representative.gene_lengths.tsv")

tpm.df.fname.rds <- file.path(mag.rdafiles, "tpm.df.rds")
tpm.df.fname.tsv <- file.path(mag.tables, "tpm.df.tsv")
# seqkit.all.df.fname.no_ext <- "seqkit.all.df" 
# gtdbtk.coverm.drep.fname.no_ext <- "gtdbtk.coverm.drep"
# gtdbtk.coverm.metabat2.fname.no_ext <- "gtdbtk.coverm.metabat2"
# blastn.coverm.megahit.fname.no_ext <- "blastn.coverm.megahit"
# gtdbtk.taxonomy.clean.fname.no_ext <- "gtdbtk.taxonomy.clean"
# gene.annotation.df.fname.no_ext <- "gene.annotation.df.rds"
# dbcan.df.fname.no_ext <- "dbcan.df"
# kofamscan.df.fname.no_ext <- "kofamscan.df"
# representative.gene_lengths.fname.no_ext <- "representative.gene_lengths"
