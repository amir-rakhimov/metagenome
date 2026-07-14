#' ---
#' output: 
#'   bookdown::html_document2:
#'      toc: true
#' ---
#' 
#' ```{r, setup 002-summary-stats-mag.R, include=FALSE}
#' knitr::opts_knit$set(root.dir = '/home/rakhimov/projects/metagenome')
#' ```
#+ echo=FALSE
# Analysing MAG data. ####
#' # Analysing MAG data.
#' 
#+ echo=FALSE
## Introduction ####
#'
#' ## Introduction
#' In this script, we will explore the MAG dataset using tidyverse and phyloseq.
#' 
#' 
#+ echo=FALSE
## 1. Load necessary libraries. ####
#'
#' ## Load necessary libraries.
# install.packages(c("tidyverse"))
# BiocManager::install("phyloseq")
# install.packages(
#   "microViz",
#   repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
# )
library(tidyverse)
library(phyloseq)
library(microViz)
#+ echo=FALSE
## 2. Import the datasets. ####
#'
#' ## Import the datasets.
source(here::here("config/R/config.R"))# config file with global variables
source(here::here("config/R/themes.R"))# config file with themes
#' Seqkit information on every contig and bin:
seqkit.all.df <- readRDS(file = seqkit.all.df.fname.rds)
#' Taxonomy classification and relative abundance of each representative bin from dRep
gtdbtk.coverm.drep <- readRDS(file = gtdbtk.coverm.drep.fname.rds)
#' Taxonomy classification and relative abundance of each bin from MetaBAT2
gtdbtk.coverm.metabat2 <- readRDS(file = gtdbtk.coverm.metabat2.fname.rds)
#' Taxonomy classification and relative abundance of each contig from BLASTN
blastn.coverm.megahit <- readRDS(file = blastn.coverm.megahit.fname.rds)

#+ echo=FALSE
## 3. Basic stats ####
#'
#' ## Basic stats
#' On MEGAHIT data (contigs)
megahit.summary.stats <- seqkit.all.df %>% #seqkit.with_taxonomy%>%
  summarise(total_len=sum(contig_length),
         longest_contig_len=max(contig_length),
         shortest_contig_len=min(contig_length),
         average_len=mean(contig_length),
         n_contig=n_distinct(sample_id,contig_id))
megahit.summary.stats %>%
  knitr::kable(format = "simple")%>%
  print()

#' On MetaBAT2 data (bins)
metabat2.summary.stats <- seqkit.all.df %>% # seqkit.with_taxonomy%>%
  filter(!secondary_cluster %in% c("unassembled"))%>%
  summarise(total_len=sum(contig_length),
         longest_contig_len=max(contig_length),
         shortest_contig_len=min(contig_length),
         average_len=mean(contig_length),
         n_contig=n_distinct(sample_id,contig_id)) 
metabat2.summary.stats%>%
  knitr::kable(format = "simple")%>%
  print()

#' On dRep data (high-quality bins)
drep.summary.stats <- seqkit.all.df %>% # seqkit.with_taxonomy%>%
  filter(is_representative_bin == TRUE)%>%
  summarise(total_len=sum(contig_length),
         longest_contig_len=max(contig_length),
         shortest_contig_len=min(contig_length),
         average_len=mean(contig_length),
         n_contig=n_distinct(sample_id,contig_id))
drep.summary.stats%>%
  knitr::kable(format = "simple")%>%
  print()

#' Combine stats into one Summary statistics from assembly.
final.summary.stats<-rbind(megahit.summary.stats,
                           metabat2.summary.stats,
                           drep.summary.stats)%>%
  mutate(source=c("megahit","metabat2","drep"))%>%
  mutate(average_len=round(average_len))
final.summary.stats%>%
  knitr::kable(format = "simple")%>%
  print()

write.table(final.summary.stats,
            file = "./output/rtables/mag-summary-stats.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = F)
rm(final.summary.stats)
rm(drep.summary.stats)
rm(metabat2.summary.stats)
rm(megahit.summary.stats)

#' Number of representative high-quality bins (MAGs): 128.
seqkit.all.df %>%
  filter(is_representative_bin == TRUE)%>%
  distinct(sample_id,bin_id)%>%
  nrow

#' Length distribution:
seqkit.all.df %>%
  filter(is_representative_bin == TRUE)%>%
  ggplot(aes(x=contig_length))+
  geom_histogram()

#' Number of contigs per bin:
seqkit.all.df %>%
  filter(is_representative_bin == TRUE)%>%
  group_by(sample_id,bin_id)%>%
  summarise(n_contig_per_bin=n_distinct(contig_id))%>%
  group_by(sample_id)%>%
  summarise(mean_n_contig=mean(n_contig_per_bin))

#+ echo=FALSE
## 4. Fix the unclassified cells using tax_fix from microViz. ####
#'
#' ## Fix the unclassified cells using tax_fix from microViz.
#' Need to convert the taxonomy table into a phyloseq tax_table
gtdbtk.drep.tax_fix <- gtdbtk.coverm.drep%>%
  mutate(sample_bin_cluster=paste(sample_id,"bin",bin_id,"cluster",secondary_cluster,sep = "_"))%>%
  select(classification,sample_bin_cluster)%>%
  rename(OTU=sample_bin_cluster)%>%
  separate_wider_delim(delim = ";",cols = classification,
                       names =  c("Kingdom","Phylum","Class","Order",
                                  "Family","Genus","Species"))%>%
  as.data.frame()%>%
  column_to_rownames("OTU")%>%
  as.matrix()
head(gtdbtk.drep.tax_fix)%>%
  knitr::kable(format = "simple")%>%
  print()
gtdbtk.drep.tax_fix <- tax_table(gtdbtk.drep.tax_fix)
gtdbtk.drep.tax_fix <- tax_fix(gtdbtk.drep.tax_fix,
                                    unknowns = "unclassified")
gtdbtk.drep.tax_fix <- gtdbtk.drep.tax_fix%>%
  as.data.frame()%>%
  rownames_to_column(var="sample_bin_cluster")%>%
  as_tibble()%>%
  separate_wider_delim(sample_bin_cluster,delim = "_bin_",names = c("sample_id","bin_cluster"))%>%
  separate_wider_delim(bin_cluster,delim = "_cluster_",names = c("bin_id","secondary_cluster"))%>%
  mutate(bin_id=as.integer(bin_id))
head(gtdbtk.drep.tax_fix)%>%
  knitr::kable(format = "simple")%>%
  print()

if(!file.exists(gtdbtk.taxonomy.clean.fname.rds) & !file.exists(gtdbtk.taxonomy.clean.fname.tsv)){
  saveRDS(gtdbtk.drep.tax_fix,
          file = gtdbtk.taxonomy.clean.fname.rds)
  write.table(file = gtdbtk.taxonomy.clean.fname.tsv,
              gtdbtk.drep.tax_fix,
              sep = '\t',
              row.names = F,
              col.names = T)
}

#+ echo=FALSE
## 5. Classification rate in each sample. ####
#'
#' ## Classification rate in each sample.
gtdbtk.drep.classification.rate <- gtdbtk.coverm.drep%>%
  group_by(sample_id)%>%
  summarise(sum_relative_abundance=sum(relative_abundance))%>%
  mutate(mean_classified=mean(sum_relative_abundance),
         sd_classified=sd(sum_relative_abundance))
gtdbtk.drep.classification.rate %>%
  knitr::kable(format = "simple")%>%
  print()

gtdbtk.drep.classification.rate.fname <- file.path(mag.tables, "gtdbtk-drep-classification-rate.tsv")
write.table(gtdbtk.drep.classification.rate,
            file = gtdbtk.drep.classification.rate.fname,
            sep = "\t",
            col.names = T,
            row.names = F)

#+ echo=FALSE
## 6. Number of phyla, families, genera, species. ####
#'
#' ## Number of phyla, families, genera, species.
#' In GTDB-tk:
gtdbtk.drep.n_taxa <- gtdbtk.drep.tax_fix%>%
  summarise(PhylaPerHost=n_distinct(Phylum),
            FamiliesPerHost=n_distinct(Family),
            GeneraPerHost=n_distinct(Genus),
            SpeciesPerHost=n_distinct(secondary_cluster))
gtdbtk.drep.n_taxa %>%
  knitr::kable(format = "simple")%>%
  print()

gtdbtk.drep.n_taxa.fname <- file.path(mag.tables, "mag-n_taxa-gtdbtk-drep.tsv")
if(!file.exists(gtdbtk.drep.n_taxa.fname)){
  write.table(gtdbtk.drep.n_taxa,
              file = gtdbtk.drep.n_taxa.fname,
              sep = "\t",
              col.names = T,
              row.names = F)
}

#' Recalculate relative abundances as relative abundance/sum(relative 
#' abundances with sample).
gtdbtk.coverm.drep.rescaled <- gtdbtk.coverm.drep%>%
  full_join(gtdbtk.drep.tax_fix,
            by = join_by(sample_id,bin_id,secondary_cluster))%>%
  group_by(sample_id)%>%
  mutate(old_relative_abundance=relative_abundance,
         relative_abundance=old_relative_abundance/sum(old_relative_abundance)*100)%>%
  ungroup

#' Dominant phyla (average relative abundance) and families within those phyla.
#' Rescaled is new, raw is old.
get_dominant_taxa_mag <- function(df, tax.rank){
  if(tax.rank == "Phylum"){
    tax.ranks.to_keep <- "Phylum"
  }else if(tax.rank == "Family"){
    tax.ranks.to_keep <- c( "Phylum", "Family")
  }else if(tax.rank == "Genus"){
    tax.ranks.to_keep <- c("Phylum", "Family", "Genus")
  }
  cols.to_keep <-c(tax.ranks.to_keep,
                   "min_old","max_old","mean_old","sd_old","n_mags",
                   "min_new","max_new","mean_new","sd_new")
  
  dominant.taxa <- df %>%
    group_by_at(tax.rank)%>%
    mutate(min_old=min(old_relative_abundance),
           max_old=max(old_relative_abundance),
           mean_old=mean(old_relative_abundance),
           sd_old=sd(old_relative_abundance),
           n_mags = n(),
           min_new=min(relative_abundance),
           max_new=max(relative_abundance),
           mean_new=mean(relative_abundance),
           sd_new=sd(relative_abundance))%>%
    ungroup%>%
    distinct(dplyr::across(all_of(cols.to_keep)))%>%
    arrange(-n_mags)
  return(dominant.taxa)
}

dominant.phyla.gtdbtk.drep <- 
  get_dominant_taxa_mag(gtdbtk.coverm.drep.rescaled,"Phylum")
  
dominant.phyla.gtdbtk.drep %>%
  knitr::kable(format = "simple")%>%
  print()

dominant.families.gtdbtk.drep<-
  get_dominant_taxa_mag(gtdbtk.coverm.drep.rescaled,"Family")

head(dominant.families.gtdbtk.drep) %>%
  knitr::kable(format = "simple")%>%
  print()

dominant.genera.gtdbtk.drep<-
  get_dominant_taxa_mag(gtdbtk.coverm.drep.rescaled,"Genus")

head(dominant.genera.gtdbtk.drep)%>%
  knitr::kable(format = "simple")%>%
  print()

dominant.phyla.gtdbtk.drep.fname <- file.path(mag.tables,"mag-dominant-phyla-gtdbtk-drep.tsv")
dominant.families.gtdbtk.drep.fname <- file.path(mag.tables,"mag-dominant-families-gtdbtk-drep.tsv")
dominant.genera.gtdbtk.drep.fname <- file.path(mag.tables,"mag-dominant-genera-gtdbtk-drep.tsv")
if(!file.exists(dominant.phyla.gtdbtk.drep.fname)){
  write.table(dominant.phyla.gtdbtk.drep,
              file = dominant.phyla.gtdbtk.drep.fname,
              sep = "\t",
              col.names = T,
              row.names = F)
  
}
if(!file.exists(dominant.families.gtdbtk.drep.fname)){
  write.table(dominant.families.gtdbtk.drep,
              file=dominant.families.gtdbtk.drep.fname,
              sep = "\t",
              col.names = T,
              row.names = F)
  
}

if(!file.exists(dominant.genera.gtdbtk.drep.fname)){
  write.table(dominant.genera.gtdbtk.drep,
              file = dominant.genera.gtdbtk.drep.fname,
              sep = "\t",
              col.names = T,
              row.names = F)  
}

#+ echo=FALSE
## 7. Abundance of specific taxa. ####
#'
#' ## Abundance of specific taxa.
dominant.phyla.gtdbtk.drep%>%
  filter(grepl("Bacteroidota|Spirochaetota|Desulfobacterota",Phylum))%>%
  knitr::kable(format = "simple")%>%
  print()
dominant.families.gtdbtk.drep%>%
  filter(grepl("Bacteroidaceae",Family))%>%
  knitr::kable(format = "simple")%>%
  print()
dominant.genera.gtdbtk.drep%>%
  filter(grepl("Treponema",Genus))%>%
  knitr::kable(format = "simple")%>%
  print()

#+ echo=FALSE
## 8. Count eukaryotic contigs. ####
#'
#' ## Count eukaryotic contigs.
#' Count genera (after the second semicolon from the end):
blastn.coverm.megahit%>%
  filter(grepl("Eukaryota",classification))%>%
  mutate(classification=sub(".*;[^;]*;([^;]*;[^;]*)$", "\\1", classification))%>%
  separate_wider_delim(delim = ";",cols=classification,
                       names=c("Genus","Species"))%>%
  count(Genus,name = "n_contigs")%>%
  arrange(-n_contigs)%>%
  knitr::kable(format = "simple")%>%
  print()
sessionInfo()
rm(list = ls(all=TRUE))
gc()
