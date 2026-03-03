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
#' Seqkit statistics on every contig and bin with taxonomy classification and
#' relative abudnance from CoverM:
seqkit.with_taxonomy <- readRDS(file = "./output/rdafiles/seqkit-with_taxonomy.rds")
#' Taxonomy classification of each bin from dRep
gtdbtk.coverm.drep<-readRDS(file = "./output/rdafiles/gtdbtk-coverm-drep.rds")
#' Taxonomy classification of each bin from MetaBAT2
gtdbtk.coverm.metabat2<-readRDS(file = "./output/rdafiles/gtdbtk-coverm-metabat2.rds")
#' Taxonomy classification of each contig from BLASTN
blastn.coverm.megahit<-readRDS(file = "./output/rdafiles/blastn-coverm-megahit.rds")

#+ echo=FALSE
## 3. Basic stats ####
#'
#' ## Basic stats
#' On MEGAHIT data (contigs)
megahit.summary.stats<-seqkit.with_taxonomy%>%
  mutate(total_len=sum(contig_length),
         longest_contig_len=max(contig_length),
         shortest_contig_len=min(contig_length),
         average_len=mean(contig_length),
         n_contig=n_distinct(sample,contig_id))%>%
  # select(-sample,-contig_id,-contig_length,-GC,-AT)%>%
  select(total_len,longest_contig_len,shortest_contig_len,average_len,
         n_contig)%>%
  distinct()
megahit.summary.stats

#' On MetaBAT2 data (bins)
metabat2.summary.stats<-seqkit.with_taxonomy%>%
  filter(!is.na(bin_id))%>%
  mutate(total_len=sum(contig_length),
         longest_contig_len=max(contig_length),
         shortest_contig_len=min(contig_length),
         average_len=mean(contig_length),
         n_contig=n_distinct(sample,contig_id))%>%
  # select(-sample,-contig_id,-contig_length,-GC,-AT)%>%
  select(total_len,longest_contig_len,shortest_contig_len,average_len,
         n_contig)%>%
  distinct()
metabat2.summary.stats

#' On dRep data (high-quality bins)
drep.summary.stats<-seqkit.with_taxonomy%>%
  filter(is_drep_bin)%>%
  mutate(total_len=sum(contig_length),
         longest_contig_len=max(contig_length),
         shortest_contig_len=min(contig_length),
         average_len=mean(contig_length),
         n_contig=n_distinct(sample,contig_id))%>%
  # select(-sample,-contig_id,-contig_length,-GC,-AT)%>%
  select(total_len,longest_contig_len,shortest_contig_len,average_len,
         n_contig)%>%
  distinct()
drep.summary.stats

#' Combine stats into one Summary statistics from assembly.
final.summary.stats<-rbind(megahit.summary.stats,
                           metabat2.summary.stats,
                           drep.summary.stats)%>%
  mutate(source=c("megahit","metabat2","drep"))%>%
  mutate(average_len=round(average_len))
final.summary.stats
# write.table(final.summary.stats,
#             file = "./output/rtables/mag-summary-stats.tsv",
#             sep = "\t",
#             quote = FALSE,
#             row.names = F)
rm(final.summary.stats)
rm(drep.summary.stats)
rm(metabat2.summary.stats)
rm(megahit.summary.stats)

#' Number of total bins (MAGs): 128.
seqkit.with_taxonomy%>%
  filter(is_drep_bin)%>%
  distinct(sample,bin_id)%>%
  nrow

#' Length distribution:
seqkit.with_taxonomy%>%
  filter(is_drep_bin)%>%
  ggplot(aes(x=contig_length))+
  geom_histogram()

#' Number of contigs per bin:
seqkit.with_taxonomy%>%
  filter(is_drep_bin)%>%
  group_by(sample,bin_id)%>%
  summarise(n_contig_per_bin=n_distinct(contig_id))%>%
  group_by(sample)%>%
  summarise(mean_n_contig=mean(n_contig_per_bin))

#+ echo=FALSE
## 4. Fix the unclassified cells using tax_fix from microViz. ####
#'
#' ## Fix the unclassified cells using tax_fix from microViz.
#' Need to convert the taxonomy table into a phyloseq tax_table
gtdbtk.coverm.drep.tax_fix<-gtdbtk.coverm.drep%>%
  mutate(sample_bin_cluster=paste(sample,"bin",bin_id,"cluster",secondary_cluster,sep = "_"))%>%
  select(classification,sample_bin_cluster)%>%
  rename(OTU=sample_bin_cluster)%>%
  separate_wider_delim(delim = ";",cols = classification,
                       names =  c("Kingdom","Phylum","Class","Order",
                                  "Family","Genus","Species"))%>%
  as.data.frame()%>%
  column_to_rownames("OTU")%>%
  as.matrix()
head(gtdbtk.coverm.drep.tax_fix)
gtdbtk.coverm.drep.tax_fix<-tax_table(gtdbtk.coverm.drep.tax_fix)
gtdbtk.coverm.drep.tax_fix<-tax_fix(gtdbtk.coverm.drep.tax_fix,
                                    unknowns = "unclassified")
gtdbtk.coverm.drep.tax_fix<-gtdbtk.coverm.drep.tax_fix%>%
  as.data.frame()%>%
  rownames_to_column(var="sample_bin_cluster")%>%
  as_tibble()%>%
  separate_wider_delim(sample_bin_cluster,delim = "_bin_",names = c("sample","bin_cluster"))%>%
  separate_wider_delim(bin_cluster,delim = "_cluster_",names = c("bin_id","secondary_cluster"))%>%
  mutate(bin_id=as.integer(bin_id))
head(gtdbtk.coverm.drep.tax_fix)
# saveRDS(gtdbtk.coverm.drep.tax_fix,
#         file="./output/rdafiles/gtdbtk-taxonomy.rds")
# write.table(file="./output/rtables/gtdbtk-taxonomy.tsv",
#             gtdbtk.coverm.drep.tax_fix,sep = '\t',
#             row.names = F,col.names = T)

#+ echo=FALSE
## 5. Classification rate in each sample. ####
#'
#' ## Classification rate in each sample.
gtdbtk.drep.classification.rate<-gtdbtk.coverm.drep%>%
  group_by(sample)%>%
  summarise(sum_relative_abundance=sum(relative_abundance))%>%
  mutate(mean_classified=mean(sum_relative_abundance),
         sd_classified=sd(sum_relative_abundance))
gtdbtk.drep.classification.rate
# write.table(gtdbtk.drep.classification.rate,
#             file="output/rtables/gtdbtk-drep-classification-rate.tsv",
#             sep = "\t",
#             col.names = T,
#             row.names = F)

#+ echo=FALSE
## 6. Number of phyla, families, genera, species. ####
#'
#' ## Number of phyla, families, genera, species.
#' In GTDB-tk:
gtdbtk.drep.n_taxa<-gtdbtk.coverm.drep.tax_fix%>%
  summarise(PhylaPerHost=n_distinct(Phylum),
            FamiliesPerHost=n_distinct(Family),
            GeneraPerHost=n_distinct(Genus),
            SpeciesPerHost=n_distinct(secondary_cluster))
gtdbtk.drep.n_taxa

# write.table(gtdbtk.drep.n_taxa,
#             file="output/rtables/mag-n_taxa-gtdbtk-drep.tsv",
#             sep = "\t",
#             col.names = T,
#             row.names = F)

#' Recalculate relative abundances as relative abundance/sum(relative 
#' abundances with sample).
gtdbtk.coverm.drep.rescaled<-gtdbtk.coverm.drep%>%
  full_join(gtdbtk.coverm.drep.tax_fix,
            by = join_by(sample,bin_id,secondary_cluster))%>%
  group_by(sample)%>%
  mutate(old_relative_abundance=relative_abundance,
         relative_abundance=old_relative_abundance/sum(old_relative_abundance)*100)%>%
  ungroup

# gtdbtk.coverm.drep%>%
#   separate_wider_delim(delim = ";",cols = classification,
#                        names =  c("Kingdom","Phylum","Class","Order",
#                                      "Family","Genus","Species"))%>%
#   group_by(sample,Phylum)%>%
#   summarise(relative_abundance=sum(relative_abundance))%>%
#   group_by(sample)%>%
#   mutate(relative_abundance=relative_abundance/sum(relative_abundance)*100)%>%
#   group_by(Phylum)%>%
#   mutate(MeanRelativeAbundance=mean(relative_abundance))%>%
#   distinct(Phylum,MeanRelativeAbundance)

#' Dominant phyla (average relative abundance) and families within those phyla.
#' Rescaled is new, raw is old.
dominant.phyla.gtdbtk.drep<-gtdbtk.coverm.drep.rescaled%>%
  group_by(Phylum)%>%
  mutate(min_old=min(old_relative_abundance),
         max_old=max(old_relative_abundance),
         mean_old=mean(old_relative_abundance),
         sd_old=sd(old_relative_abundance),
         n=n(),
         min_new=min(relative_abundance),
         max_new=max(relative_abundance),
         mean_new=mean(relative_abundance),
         sd_new=sd(relative_abundance))%>%
  ungroup%>%
  distinct(Phylum,min_old,max_old,mean_old,sd_old,n,
           min_new,max_new,mean_new,sd_new)%>%
  arrange(-n)
dominant.phyla.gtdbtk.drep

dominant.families.gtdbtk.drep<-gtdbtk.coverm.drep.rescaled%>%
  group_by(Family)%>%
  mutate(min_old=min(old_relative_abundance),
         max_old=max(old_relative_abundance),
         mean_old=mean(old_relative_abundance),
         sd_old=sd(old_relative_abundance),
         n=n(),
         min_new=min(relative_abundance),
         max_new=max(relative_abundance),
         mean_new=mean(relative_abundance),
         sd_new=sd(relative_abundance))%>%
  ungroup%>%
  distinct(Phylum,Family,min_old,max_old,mean_old,sd_old,n,
           min_new,max_new,mean_new,sd_new)%>%
  arrange(-n)
head(dominant.families.gtdbtk.drep)

dominant.genera.gtdbtk.drep<-gtdbtk.coverm.drep.rescaled%>%
  group_by(Genus)%>%
  mutate(min_old=min(old_relative_abundance),
         max_old=max(old_relative_abundance),
         mean_old=mean(old_relative_abundance),
         sd_old=sd(old_relative_abundance),
         n=n(),
         min_new=min(relative_abundance),
         max_new=max(relative_abundance),
         mean_new=mean(relative_abundance),
         sd_new=sd(relative_abundance))%>%
  ungroup%>%
  distinct(Phylum,Family,Genus,min_old,max_old,mean_old,sd_old,n,
           min_new,max_new,mean_new,sd_new)%>%
  arrange(-n)
head(dominant.genera.gtdbtk.drep)

# write.table(dominant.phyla.gtdbtk.drep,
#             file="output/rtables/mag-dominant-phyla-gtdbtk-drep.tsv",
#             sep = "\t",
#             col.names = T,
#             row.names = F)
# 
# write.table(dominant.families.gtdbtk.drep,
#             file="output/rtables/mag-dominant-families-gtdbtk-drep.tsv",
#             sep = "\t",
#             col.names = T,
#             row.names = F)
# 
# write.table(dominant.genera.gtdbtk.drep,
#             file="output/rtables/mag-dominant-genera-gtdbtk-drep.tsv",
#             sep = "\t",
#             col.names = T,
#             row.names = F)

#+ echo=FALSE
## 7. Abundance of specific taxa. ####
#'
#' ## Abundance of specific taxa.
dominant.phyla.gtdbtk.drep%>%
  filter(grepl("Bacteroidota|Spirochaetota|Desulfobacterota",Phylum))
dominant.families.gtdbtk.drep%>%
  filter(grepl("Bacteroidaceae",Family))
dominant.genera.gtdbtk.drep%>%
  filter(grepl("Treponema",Genus))

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
  arrange(-n_contigs)
sessionInfo()
rm(list = ls(all=TRUE))
gc()