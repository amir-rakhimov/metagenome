project.root <- here::here()

kraken2.date_time<-"20240409_17_32_40"
kraken2.unclassified.stats<-file.path("./output/kraken2_pipeline",
          paste(kraken2.date_time,"unclassified_reads.tsv",
                sep="_"))

#' Directories with input files:
# rdafiles.directory<-file.path(project.root,"results","rdafiles")
# rtables.directory<-file.path(project.root,"results","rtables")
#' Directory with metadata:
metadatadir<-paste0("../amplicon_nmr/data/metadata/shared") 
# Kraken2 report
combined.report.date_time<-"20240515_10_04_04"
combined.report.filename<-file.path(rtables.directory,
                                    paste(combined.report.date_time,"combined_report.tsv",sep = "_"))


kraken2.analysis.out <- function(...) file.path(project.root, "analysis", "1-kraken2")
mag.analysis.out <- function(...)  file.path(project.root, "analysis", "2-mag")
manuscript.figures <- file.path(project.root,"analysis", "manuscript-figures")


kraken2.tables <- kraken2.analysis.out("tables")
kraken2.rdafiles <- kraken2.analysis.out("rdafiles")
kraken2.figures <- kraken2.analysis.out("figures")

mag.tables <- mag.analysis.out("tables")
mag.rdafiles <- mag.analysis.out("rdafiles")
mag.figures <- mag.analysis.out("figures")



# 002
util.functions.r <- "../utils/R"
## 2. Specifying parameters and directory/file names. #### 
# Import datasets as rds files.
ps.q.raw.fname <- file.path(community.composition.rdafiles,paste(
  "phyloseq-qiime",authorname,truncationlvl,
  "ps.q.raw.rds",sep="-"))
ps.q.rel.raw.fname <- file.path(community.composition.rdafiles,paste(
  "phyloseq-qiime",authorname,truncationlvl,
  "ps.q.rel.raw.rds",sep="-"))
ps.q.raw.otu_table.fname <- file.path(community.composition.rdafiles,paste(
  "phyloseq-qiime",authorname,truncationlvl,
  "ps.q.raw.otu_table.rds",sep="-"))
ps.q.raw.tax_table.fname <- file.path(community.composition.rdafiles,paste(
  "phyloseq-qiime",authorname,truncationlvl,
  "ps.q.raw.tax_table.rds",sep="-"))
ps.q.raw.tree.fname <- file.path(community.composition.rdafiles,paste(
  "phyloseq-qiime",authorname,truncationlvl,
  "ps.q.raw.tree.rds",sep="-"))
rare.num_samples <- 100

ps.q.agg.asv.fname.no_ext <- paste("phyloseq-qiime",authorname,"OTU",read.end.type, # new 
                                  truncationlvl,"table",sep = "-")
ps.q.agg.genus.fname.no_ext <-paste("phyloseq-qiime",authorname,"Genus",read.end.type,
                                    truncationlvl,"table",sep = "-")
ps.q.agg.family.fname.no_ext <- paste("phyloseq-qiime",authorname,"Family",read.end.type,
                                      truncationlvl,"table",sep = "-")
ps.q.agg.phylum.fname.no_ext <- paste("phyloseq-qiime",authorname,"Phylum",read.end.type,
                                      truncationlvl,"table",sep = "-")

ps.q.filtered.fname <- gsub("raw", "filtered", ps.q.raw.fname)
ps.q.rel.filtered.fname <- gsub("raw", "filtered", ps.q.rel.raw.fname)
ps.q.filtered.otu_table.fname <- gsub("raw", "filtered", ps.q.raw.otu_table.fname)
ps.q.filtered.tax_table.fname <- gsub("raw", "filtered", ps.q.raw.tax_table.fname)
ps.q.filtered.tree.fname <- gsub("raw", "filtered", ps.q.raw.tree.fname)


#' Import metadata:
custom.md.path<-file.path(new.metadata.dir,"custom.md.rds")

## Setup plots.
#' Set "pretty" labels for barplot facets that correspond to animal hosts. 
#' Here, the left side of the vector (values) is taken from the metadata, 
#' while the right side (names) are the pretty labels that will be shown on 
#' the final barplot.
pretty.level.names<-c("NMR" = "*H. glaber*", # better labels for facets
                      "DMR" = "*F. damarensis*",
                      "B6mouse" = "B6 mouse",
                      "MSMmouse" = "MSM/Ms mouse",
                      "FVBNmouse" = "FVB/N mouse",
                      "spalax" = "*N. leucodon*",
                      # "pvo" = "*P. v. orii*",
                      "hare" = "*L. europaeus*",
                      "rabbit" = "*O. cuniculus*")
image.formats<-c("png","tiff")



