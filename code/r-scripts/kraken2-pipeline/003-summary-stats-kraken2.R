#' ---
#' output: 
#'   bookdown::html_document2:
#'      toc: true
#' ---
#' 
#' ```{r, setup 003-summary-stats-kraken2.R, include=FALSE}
#' knitr::opts_knit$set(root.dir = '/home/rakhimov/projects/amplicon_nmr')
#' ```
#' ```{r, echo = FALSE}
#' # For showing images, tables, etc: Use global path
#' #knitr::spin("code/r-scripts/003-summary-stats-kraken2.R", knit = FALSE)
#' #file.rename("code/r-scripts/003-summary-stats-kraken2.Rmd", "markdown/003-summary-stats-kraken2.Rmd")
#' #rmarkdown::render('./markdown/003-summary-stats-kraken2.Rmd', 'html_document',
#' # knit_root_dir="/home/rakhimov/projects/amplicon_nmr/")
#' ```
#' 
#+ echo=FALSE
# Analysing phyloseq data (Kraken2) ####
#' # Analysing phyloseq data (Kraken2)
#' 
#+ echo=FALSE
## Introduction ####
#'
#' ## Introduction
#' In this script, we will explore the imported dataset from QIIME2 (using 
#' phyloseq).
#' 
#' We will use the data from 002-phyloseq-kraken2.R script (ps.q.agg
#' agglomerated tables at phylum, family, genus, and Species level).

#+ echo=FALSE
## 1. Load necessary libraries and scripts. ####
#'
#' ## Load necessary libraries and scripts.
# install.packages(c("tidyverse","ggtext"))
library(tidyverse)
library(ggtext)
library(vegan)

#' Load necessary scripts.
source("../amplicon_nmr/code/r-scripts/get_n_uniq_taxa_per_host.R")
source("../amplicon_nmr/code/r-scripts/create_summary_stats_table.R")
source("../amplicon_nmr/code/r-scripts/add_relab_to_tax_df.R")
source("../amplicon_nmr/code/r-scripts/add_agegroup_to_tax_df.R")
source("../amplicon_nmr/code/r-scripts/get_unclassified_summary_stats.R")
source("../amplicon_nmr/code/r-scripts/ggplot_species.R")
###
source("../amplicon_nmr/code/r-scripts/get_dominant_taxa_in_host.R")
#+ echo=FALSE
## 2. Specifying parameters and directory/file names. #### 
#'
#' ## Specifying parameters and directory/file names. 
#' Directories with input files:
rdafiles.directory<-"./output/rdafiles"
rtables.directory<-"./output/rtables"

#' Import datasets as rds files
ps.q.agg.date_time<-"20241003_13_52_43"
ps.q.agg.genus.date_time<-"20241003_13_56_28"
ps.q.agg.family.date_time<-"20241003_13_57_30"
ps.q.agg.phylum.date_time<-"20241003_13_57_59"

ps.q.agg<-readRDS(file=file.path(
  rdafiles.directory,
  paste(ps.q.agg.date_time,"phyloseq","kraken2","Species",
        "table.rds",sep = "-")))
ps.q.agg.genus<-readRDS(file=file.path(
  rdafiles.directory,
  paste(ps.q.agg.genus.date_time,"phyloseq","kraken2","Genus",
        "table.rds",sep = "-")))

ps.q.agg.family<-readRDS(file=file.path(
  rdafiles.directory,
  paste(ps.q.agg.family.date_time,"phyloseq","kraken2","Family",
        "table.rds",sep = "-")))
ps.q.agg.phylum<-readRDS(file=file.path(
  rdafiles.directory,
  paste(ps.q.agg.phylum.date_time,"phyloseq","kraken2","Phylum",
        "table.rds",sep = "-")))
# Import metadata
custom.md<-readRDS("../amplicon_nmr/output/rdafiles/custom.md.ages.rds")
custom.md<-custom.md%>%
  filter(Sample%in%ps.q.agg$Sample, 
         sequencing_type =="Naked mole-rat whole metagenome sequencing")

#+ echo=FALSE
## 3. Setup plots. ####
#'
#' ## Setup plots.
sample.levels<-custom.md%>%
  ungroup()%>%
  select(Sample,age)%>%
  arrange(age)%>%
  distinct()%>%
  mutate(NewSample=paste0(Sample," (",age,")"))%>%
  mutate(Sample=factor(Sample,levels=Sample),
         NewSample=factor(NewSample,levels=NewSample))

pretty.level.names<-names(table(custom.md$old_agegroup))
names(pretty.level.names)<-names(table(custom.md$agegroup))
custom.levels<-names(pretty.level.names)
gg.labs.name<-"Age group"
gg.title.groups<-"age groups"

#+ echo=FALSE
## 4. Calculating summary statistics. ####
#'
#' ## Calculating summary statistics.
#' 
#+ echo=FALSE
### 4.1 Check the total number of unique Species/phyla/families/genera per class. ####
#' 
#' ### Check the total number of unique Species/phyla/families/genera per class.
#' Here, Species are called ASV because the same function is used on 16S rRNA
#' gene sequencing analysis
#' We will use it for the summary table in the next section.
n.species.per.host<-get_n_uniq_taxa_per_host(ps.q.agg,"Species")
n.phylum.per.host<-get_n_uniq_taxa_per_host(ps.q.agg.phylum,"Phylum")
n.family.per.host<-get_n_uniq_taxa_per_host(ps.q.agg.family,"Family")
n.genus.per.host<-get_n_uniq_taxa_per_host(ps.q.agg.genus,"Genus")
#+ echo=FALSE
### 4.2 Create a summary table. ####
#'
#' ### Create a summary table.
#' The columns are:  
#' * Total reads  
#' * Library size (mean Abundance ± SD)  
#' * Num of Species per host  
#' * Num of phyla per host  
#' * Num of families per host  
#' * Num of genera per host  
summary.stats.table<-create_summary_stats_table(ps.q.agg,
                                          n.asv.table = n.species.per.host,
                                          n.phylum.table = n.phylum.per.host,
                                          n.family.table = n.family.per.host,
                                          n.genus.table = n.genus.per.host)
summary.stats.table
# write.table(summary.stats.table,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "kraken2","summary-table.tsv",sep="-")),
#             row.names = F,sep = "\t")

#+ echo=FALSE
### 4.2 Create a summary table for each sample. ####
#'
#' ### Create a summary table for each sample.
n.taxa.per_sample<-ps.q.agg%>%
  group_by(Sample)%>%
  summarise(n_species=n_distinct(Species),
            n_genera=n_distinct(Genus),
            n_families=n_distinct(Family),
            n_phyla=n_distinct(Phylum))
n.taxa.per_sample
# write.table(n.taxa.per_sample,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "kraken2","n_taxa-per-sample-table.tsv",sep="-")),
#             row.names = F,sep = "\t")

#+ echo=FALSE
## 5. Add relative abundance and average relative abundance columns. ####
#'
#' ## Add relative abundance and average relative abundance columns.
ps.q.agg.phylum.relab<-add_relab_to_tax_df(ps.q.agg.phylum,"Phylum")
ps.q.agg.family.relab<-add_relab_to_tax_df(ps.q.agg.family,"Family")
ps.q.agg.genus.relab<-add_relab_to_tax_df(ps.q.agg.genus,"Genus")
ps.q.agg.relab<-add_relab_to_tax_df(ps.q.agg,"Species")
head(ps.q.agg.phylum.relab)
head(ps.q.agg.family.relab)
head(ps.q.agg.genus.relab)
head(ps.q.agg.relab)

#+ echo=FALSE
## 6. Add agegroup variable to NMR data (must run for plotting). ####
#'
#' ## Add agegroup variable to NMR data (must run for plotting).
ps.q.agg.relab<-add_agegroup_to_tax_df(ps.q.agg.relab,"Species",custom.md)
ps.q.agg.genus.relab<-add_agegroup_to_tax_df(ps.q.agg.genus.relab,"Genus",custom.md)
ps.q.agg.family.relab<-add_agegroup_to_tax_df(ps.q.agg.family.relab,"Family",custom.md)

#+ echo=FALSE
## 8. Rarefy the table and check the percentage of unclassified taxa. ####
#'
#' ## Rarefy the table and check the percentage of unclassified taxa. ####
#' Convert the data frame into wide format: rows are samples and columns
#' are taxa
get_rarefied_table<-function(tax.df,tax.rank,host.classes){
  tax.df.wide<-tax.df%>%
    filter(class %in% host.classes,Abundance!=0)%>%
    dplyr::select(Sample,Abundance,class,all_of(tax.rank))%>%
    filter(Abundance!=0)%>%
    pivot_wider(names_from = all_of(tax.rank),
                values_from = "Abundance",
                values_fill = 0)%>%
    as.data.frame()%>%
    column_to_rownames("Sample")%>% # Set sample names as row names
    dplyr::select(-class)
  # Find the smallest sample size
  min.n_seqs.all<-tax.df%>%
    filter(class %in% host.classes,Abundance!=0)%>%
    dplyr::select(Sample,all_of(tax.rank),Abundance)%>%
    group_by(Sample)%>%
    summarize(n_seqs=sum(Abundance))%>%
    summarize(min=min(n_seqs))%>%
    pull(min)
  print(paste("Smallest sample size:", min.n_seqs.all))
  
  ### Rarefied asv table with vegan ####
  set.seed(1)
  tax.df.rare<-rrarefy(tax.df.wide,sample=min.n_seqs.all)
  tax.df.rare<-tax.df.rare%>%
    as_tibble(rownames="Sample")%>%
    pivot_longer(-Sample)%>%
    as.data.frame()%>%
    left_join(unique(tax.df[,c("Sample","class")]),
              by="Sample")
  if(tax.rank=="OTU"){
    tax.df.rare<-tax.df.rare%>%
      rename(OTU=name,Abundance=value)%>%
      filter(Abundance!=0)  
  }else{
    # rename the 'name' column corresponding to the tax.rank
    tax.df.rare[,paste(tax.rank)]<-tax.df.rare$name
    tax.df.rare<-tax.df.rare%>%
      dplyr::select(-name)%>%
      rename(Abundance=value)%>%
      filter(Abundance!=0)
  }
  write.table(tax.df.rare,
              file = file.path(rtables.directory,paste0(
                paste(
                  paste(format(Sys.time(),format="%Y%m%d"),
                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                  "kraken2-ps.q.df.rare-nonfiltered",tax.rank,
                  paste(host.classes,collapse = '-'),sep = "-"),
                ".tsv")),
              row.names = F,
              sep = "\t")
  saveRDS(tax.df.rare,
          file = file.path(rdafiles.directory,paste0(
            paste(
              paste(format(Sys.time(),format="%Y%m%d"),
                    format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
              "kraken2-ps.q.df.rare-nonfiltered",tax.rank,
              paste(host.classes,collapse = '-'),sep = "-"),
            ".rds")))
  return(tax.df.rare)
}

ps.q.agg.rare<-get_rarefied_table(ps.q.agg,"Species","NMR")
head(ps.q.agg.rare)

#+ echo=FALSE
### 8.1 Add relative abundances and taxonomic information to the rarefied dataframe. ####
#' 
#' ### Add relative abundances and taxonomic information to the rarefied dataframe. 
#' NMR (Species)
ps.q.agg.rare.relab<-add_relab_to_tax_df(ps.q.agg.rare,"Species")
head(ps.q.agg.rare.relab)

#' Add other taxonomic ranks to the dataframe
ps.q.agg.rare.relab<-ps.q.agg.rare.relab%>%
  left_join(unique(ps.q.agg[,c("Kingdom","Phylum","Class","Order","Family",
                               "Genus","Species","OTU")]))
head(ps.q.agg.rare.relab)

#+ echo=FALSE
### 8.3 Create a summary stats table for the rarefied dataframe. ####
#' 
#' ### Create a summary stats table for the rarefied dataframe.
n.asv.per.host.rare<-get_n_uniq_taxa_per_host(ps.q.agg.rare.relab,"Species")
n.phylum.per.host.rare<-get_n_uniq_taxa_per_host(ps.q.agg.rare.relab,"Phylum")
n.family.per.host.rare<-get_n_uniq_taxa_per_host(ps.q.agg.rare.relab,"Family")
n.genus.per.host.rare<-get_n_uniq_taxa_per_host(ps.q.agg.rare.relab,"Genus")

summary.stats.table.rare<-create_summary_stats_table(ps.q.agg.rare,
                                                     n.asv.per.host.rare,
                                                     n.phylum.per.host.rare,
                                                     n.family.per.host.rare,
                                                     n.genus.per.host.rare)
summary.stats.table.rare
summary.stats.table

#+ echo=FALSE
## 10. Check the most abundant phyla, families, genera in NMR. ####
#'
#' ## Check the most abundant phyla, families, genera in NMR.
#' Phyla:
ps.q.agg.dominant.phyla<-ps.q.agg.phylum.relab%>%
  group_by(class,Phylum)%>%
  mutate(min=min(RelativeAbundance),
         max=max(RelativeAbundance),
         n=n())%>%
  distinct(class,Kingdom,Phylum, MeanRelativeAbundance,sdRelativeAbundance,
           min, max, n)%>%
  group_by(class,Kingdom,Phylum)%>%
  arrange(class,desc(MeanRelativeAbundance))%>%
  ungroup()
head(ps.q.agg.dominant.phyla)

#' Families
ps.q.agg.dominant.families<-ps.q.agg.family.relab%>%
  group_by(class,Family)%>%
  mutate(min=min(RelativeAbundance),
         max=max(RelativeAbundance),
         n=n())%>%
  distinct(class, Kingdom,Phylum, Family, MeanRelativeAbundance, sdRelativeAbundance, 
           min,max,n)%>%
  group_by(class, Kingdom,Phylum,Family)%>%
  arrange(class,desc(MeanRelativeAbundance))%>%
  ungroup()
head(ps.q.agg.dominant.families)

#' Genera:
ps.q.agg.dominant.genera<-ps.q.agg.genus.relab%>%
  group_by(class, Genus)%>%
  mutate(min=min(RelativeAbundance),
         max=max(RelativeAbundance),
         n=n())%>%
  distinct(class, Kingdom,Phylum,Family,Genus, MeanRelativeAbundance,sdRelativeAbundance,
           min, max, n)%>%
  group_by(class, Kingdom, Phylum,Family,Genus)%>%
  arrange(class,desc(MeanRelativeAbundance))%>%
  ungroup()
head(ps.q.agg.dominant.genera)

#' Species:
ps.q.agg.dominant.species<-ps.q.agg.relab%>%
  group_by(class, Species)%>%
  mutate(min=min(RelativeAbundance),
         max=max(RelativeAbundance),
         n=n())%>%
  distinct(class, Kingdom,Phylum,Family,Genus,Species, 
           MeanRelativeAbundance,sdRelativeAbundance,
           min, max, n)%>%
  group_by(class, Kingdom, Phylum,Family,Genus,Species)%>%
  arrange(class,desc(MeanRelativeAbundance))%>%
  ungroup()
head(ps.q.agg.dominant.species)

#' Non-bacterial phyla:
ps.q.agg.dominant.phyla.nonbact<-ps.q.agg.dominant.phyla%>%
  filter(Kingdom != "Bacteria")
#' Non-bacterial families:
ps.q.agg.dominant.families.nonbact<-ps.q.agg.dominant.families%>%
  filter(Kingdom != "Bacteria")
#' Non-bacterial genera:
ps.q.agg.dominant.genera.nonbact<-ps.q.agg.dominant.genera%>%
  filter(Kingdom != "Bacteria")
#' Non-bacterial species:
ps.q.agg.dominant.species.nonbact<-ps.q.agg.dominant.species%>%
  filter(Kingdom != "Bacteria")

write.table(ps.q.agg.dominant.phyla,
            file=file.path(rtables.directory,
                           paste("20241003_15_45_51",#paste(format(Sys.time(),format="%Y%m%d"),
                                 # format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "kraken2","dominant-phyla.tsv",sep="-")),
            row.names = F,sep = "\t")

write.table(ps.q.agg.dominant.phyla.nonbact,
            file=file.path(rtables.directory,
                           paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "kraken2","dominant-phyla-nonbact.tsv",sep="-")),
            row.names = F,sep = "\t")

write.table(ps.q.agg.dominant.families,
            file=file.path(rtables.directory,
                           paste("20241003_15_47_55",#paste(format(Sys.time(),format="%Y%m%d"),
                           #             format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "kraken2","dominant-families.tsv",sep="-")),
            row.names = F,sep = "\t")

write.table(ps.q.agg.dominant.families.nonbact,
            file=file.path(rtables.directory,
                           paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "kraken2","dominant-families-nonbact.tsv",sep="-")),
            row.names = F,sep = "\t")
write.table(ps.q.agg.dominant.genera,
            file=file.path(rtables.directory,
                           paste("20241003_15_48_15",#paste(format(Sys.time(),format="%Y%m%d"),
                                 # format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "kraken2","dominant-genera.tsv",sep="-")),
            row.names = F,sep = "\t")
write.table(ps.q.agg.dominant.genera.nonbact,
            file=file.path(rtables.directory,
                           paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "kraken2","dominant-genera-nonbact.tsv",sep="-")),
            row.names = F,sep = "\t")
write.table(ps.q.agg.dominant.species,
            file=file.path(rtables.directory,
                           paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "kraken2","dominant-species-all.tsv",sep="-")),
            row.names = F,sep = "\t")
write.table(ps.q.agg.dominant.species.nonbact,
            file=file.path(rtables.directory,
                           paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "kraken2","dominant-species-nonbact.tsv",sep="-")),
            row.names = F,sep = "\t")




#+ echo=FALSE
## 11. Check how much Bacteroidaceae are in NMR.  ####
#'
#' ## Check how much Bacteroidaceae are in NMR.
bacteroidaceae.nmr<-ps.q.agg.dominant.families%>%
  filter(Family=="Bacteroidaceae")
bacteroidaceae.nmr
# write.table(bacteroidaceae.nmr,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "kraken2","bacteroidaceae-nmr-table.tsv",sep="-")),
#             row.names = F,sep = "\t")

#+ echo=FALSE
### 11.1 Check the most dominant Bacteroidota families in NMR. ####
#'
#'### Check the most dominant Bacteroidota families in NMR.
bacteroidota.nmr<-ps.q.agg.dominant.families%>%
  filter(Phylum=="Bacteroidota")
bacteroidota.nmr
# write.table(bacteroidota.nmr,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "kraken2","bacteroidota-nmr-table.tsv",sep="-")),
#             row.names = F,sep = "\t")

#+ echo=FALSE
## 12. Check Spirochaetaceae, Spirochaetota, and Treponema in NMR. ####
#'
#' ## Check Spirochaetaceae, Spirochaetota, and Treponema in NMR. 
spirochaetaceae.nmr<-ps.q.agg.dominant.families%>%
  filter(Family=="Spirochaetaceae")
spirochaetaceae.nmr

spirochaetota.nmr<-ps.q.agg.dominant.phyla%>%
  filter(Phylum=="Spirochaetota")
spirochaetota.nmr
# write.table(spirochaetaceae.nmr,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "kraken2", "spirochaetaceae-nmr-table.tsv",sep="-")),
#             row.names = F,sep = "\t")
# write.table(spirochaetota.nmr,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "kraken2", "spirochaetota-nmr-table.tsv",sep="-")),
#             row.names = F,sep = "\t")

treponema.nmr<-ps.q.agg.dominant.genera%>%
  filter(Genus=="Treponema")
treponema.nmr

treponema.nmr.names<-ps.q.agg.dominant.species%>%
  filter(Genus=="Treponema")
treponema.nmr.names
# write.table(treponema.nmr,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "kraken2","treponema-nmr-table.tsv",sep="-")),
#             row.names = F,sep = "\t")
# write.table(treponema.nmr.names,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "kraken2","treponema-nmr-names.tsv",sep="-")),
#             row.names = F,sep = "\t")

#+ echo=FALSE
### 12.1 Check the number of Species in Treponema from NMR. ####
#' 
#' ### Check the number of Species in Treponema from NMR.
ps.q.agg%>%
  ungroup()%>%
  filter(Genus=="Treponema",class=="NMR")%>%
  distinct(OTU)%>%
  tally

#+ echo=FALSE
## 13. Check Mogibacteriaceae (renamed to Anaerovoracaceae). ####
#' 
#' ## Check Mogibacteriaceae (renamed to Anaerovoracaceae).
mogibacteriaceae_anaerovoracaceae.all<-ps.q.agg.dominant.families%>%
  filter(Family=="Anaerovoracaceae")
# write.table(mogibacteriaceae_anaerovoracaceae.all,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                 "kraken2", "mogibacteriaceae_anaerovoracaceae-all-table.tsv",sep="-")),
#             row.names = F,sep = "\t")

#+ echo=FALSE
## 14. Analyse sulfur-metabolising bacteria in NMR. ####
#'
#' ## Analyse sulfur-metabolising bacteria in NMR. 
#' We don't know where the "sulf" pattern is (Phylum, Order, Genus, or anywhere 
#' else), so we use sapply with grepl.
#' https://stackoverflow.com/questions/47941680/grepl-across-multiple-specified-columns
all.ranks<-c("Kingdom", "Phylum", "Class", "Order", "Family","Genus","Species")
sulfur.bact.nmr<-
  ps.q.agg[!!rowSums(sapply(ps.q.agg.relab[,all.ranks], grepl, pattern = "sulf|thio") ),]%>%
  distinct(OTU,.keep_all = T)%>%
  select(all_of(all.ranks))%>%
  inner_join(ps.q.agg.dominant.species)
# write.table(sulfur.bact.nmr,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "kraken2","sulfur.bact.nmr-nmr-table.tsv",sep="-")),
#             row.names = F,sep = "\t")

#+ echo=FALSE
### 14.2 Plot sulfur metabolising bacteria. ####
#'
#' ### Plot sulfur metabolising bacteria.
# Add zero rows
# ps.q.agg.relab<-ps.q.agg.relab%>%ungroup()
sulfur.bact.plot<-ps.q.agg.relab%>%
  filter(Species %in% sulfur.bact.nmr$Species)%>%
  mutate(Species = factor(Species, levels = sulfur.bact.nmr$Species),
         Species = gsub("_"," ",Species),
         Species = paste0("<i>",Species,"</i>")
  )%>%
  left_join(custom.md[c("Sample","agegroup")])%>%
  mutate(Sample=factor(Sample,levels=sample.levels$Sample))%>%
  group_by_at(c("class","Species"))%>%
  ggplot(aes(x=Sample,
             y=RelativeAbundance,
             fill=factor(agegroup)))+
  geom_bar(stat="identity")+
  facet_wrap(~Species,
             scales = "free",
             ncol = 3)+
  theme_bw()+
  labs(x="",
       y="Relative abundance (%)",
       fill=gg.labs.name)+
  coord_cartesian(expand = c("bottom" = F))+
  scale_color_manual(breaks = unname(pretty.level.names),
                     labels=unname(pretty.level.names))+
  scale_x_discrete(labels=sample.levels$Sample,
                   limits=sample.levels$Sample)+
  # scale_fill_manual(labels=pretty.level.names)+
  scale_fill_viridis_d(option = "C")+
  theme(axis.title.y = element_text(size = 15),
        axis.title = element_text(size = 10),
        axis.text.y = ggtext::element_markdown(size=10),
        axis.text.x = element_text(size=6),
        strip.text.x = ggtext::element_markdown(size=10),
        plot.title = element_text(size = 17),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 15),
        legend.position = "right",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
print(sulfur.bact.plot)
  # ggtitle(paste0("Relative abundance of sulfur-utilizing bacteria in different naked mole-rat age groups"))
# for (image.format in c("png","tiff")){
#   ggsave(paste0("./images/barplots/",
#                 paste(paste(format(Sys.time(),format="%Y%m%d"),
#                             format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                       "sulfur-bacteria-nmr",
#                       sep = "-"),".",image.format),
#          plot=sulfur.bact.plot,
#          width=11, height=8,units="in",
#          # width = 7000,height = 9000, units = "px",
#          dpi=300,device = image.format)
# }

treponema.sp.nmr<-ps.q.agg.relab%>%
  filter(Genus=="Treponema",class=="NMR")%>%
  group_by(Species)%>%
  select(Species,MeanRelativeAbundance)%>%
  arrange(-MeanRelativeAbundance)%>%
  ungroup()%>%
  distinct()%>%
  pull(Species)

treponema.plot<-ps.q.agg.relab%>%
  filter(Species %in% treponema.sp.nmr)%>%
  mutate(Species = factor(Species, levels = treponema.sp.nmr),
         Species = gsub("_"," ",Species),
         Species = paste0("<i>",Species,"</i>")
  )%>%
  left_join(custom.md[c("Sample","agegroup")])%>%
  mutate(Sample=factor(Sample,levels=sample.levels$Sample))%>%
  group_by_at(c("class","Species"))%>%
  ggplot(aes(x=Sample,
             y=RelativeAbundance,
             fill=factor(agegroup)))+
  geom_bar(stat="identity")+
  facet_wrap(~Species,
             scales = "free",
             ncol = 3)+
  theme_bw()+
  labs(x="",
       y="Relative abundance (%)",
       fill=gg.labs.name)+
  coord_cartesian(expand = c("bottom" = F))+
  scale_color_manual(breaks = unname(pretty.level.names),
                     labels=unname(pretty.level.names))+
  scale_x_discrete(labels=sample.levels$Sample,
                   limits=sample.levels$Sample)+
  scale_fill_viridis_d(option = "C")+
  theme(axis.title.y = element_text(size = 15),
        axis.title = element_text(size = 10),
        axis.text.y = ggtext::element_markdown(size=10),
        axis.text.x = element_text(size=6),
        strip.text.x = ggtext::element_markdown(size=10),
        plot.title = element_text(size = 17),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 15),
        legend.position = "right",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

# for (image.format in c("png","tiff")){
#   ggsave(paste0("./images/barplots/",
#                 paste(paste(format(Sys.time(),format="%Y%m%d"),
#                             format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                       "treponema.sp-nmr",
#                       sep = "-"),".",image.format),
#          plot=treponema.plot,
#          width=11, height=8,units="in",
#          # width = 7000,height = 9000, units = "px",
#          dpi=300,device = image.format)
# }




## Groups 1, 2, and 3 (**Supplementary figure**) ####
group1.genera<-c("Faecalibacterium",
                 "Roseburia",
                 "Coprococcus",
                 "Prevotella")
group1.species<-c("Eubacterium_rectale") # Not in our data
group2.increased.genera<-c("Eggerthella",
                           "Bilophila",
                           "Desulfovibrio",
                           "Fusobacterium",
                           "Anaerotruncus",
                           "Streptococcus",
                           "Escherichia")
group2.unhealthy.genera<-c("Eggerthella",
                           "Coprobacillus",
                           "Streptococcus",
                           "Bilophila",
                           "Actinomyces",
                           "Desulfovibrio",
                           "Campylobacter",
                           "Veillonella",
                           "Enterococcus")
group2.unhealthy.families<-c("Atopobiaceae",
                             "Enterobacteriaceae")
group2.unhealthy.species<-c("Bacteroides_fragilis",
                            "Clostridium_hathewayi",
                            "Enterocloster_clostridioformis", #Clostridium_clostridioforme
                            "Enterocloster_bolteae", # Clostridium_bolteae
                            "Clostridium_cindens",
                            "[Ruminococcus]_torques", #Ruminococcus_torques
                            "Mediterraneibacter_gnavus", # Ruminococcus_gnavus
                            "Clostridioides_difficile")

group3.genera<-c("Akkermansia",
                 "Odoribacter",
                 "Butyricimonas",
                 "Butyrivibrio",
                 "Barnesiella",
                 "Oscillospira")
group3.families<-c("Christensenellaceae")

group1.genera
group1.genera.plot<-ggplot_species(group1.genera,
                                   ps.q.agg.genus.relab,
                                   tax.rank="Genus",
                                   sample.order=sample.levels,
                                   group.names = pretty.level.names,
                                   grouping.variable = "agegroup",
                                   metadata.df = custom.md,
                                   ggplot.fill.name = gg.labs.name)
#+ plot.width=10, plot.height=7
print(group1.genera.plot+
        ggtitle(paste0("Relative abundance of Group 1 members"))) # 4 plots

group1.species%in%ps.q.agg.relab$Species # no need to plot

#' Group 2 increased genera
group2.increased.genera
group2.increased.genera.plot<-ggplot_species(taxa.to.plot =group2.increased.genera,
                                             tax.df =ps.q.agg.genus.relab,
                                             tax.rank="Genus",
                                             sample.order=sample.levels,
                                             group.names = pretty.level.names,
                                             grouping.variable = "agegroup",
                                             metadata.df = custom.md,
                                             ggplot.fill.name = gg.labs.name)
#+ plot.width=10, plot.height=7
print(group2.increased.genera.plot+
        ggtitle(paste0("Relative abundance of Group 2 members increased with age")))

#' Group 2 unhealthy genera
group2.unhealthy.genera
group2.unhealthy.genera.plot<-ggplot_species(taxa.to.plot=group2.unhealthy.genera,
                                              tax.df=ps.q.agg.genus.relab,
                                              tax.rank="Genus",
                                              sample.order=sample.levels,
                                              group.names = pretty.level.names,
                                              grouping.variable = "agegroup",
                                              metadata.df = custom.md,
                                              ggplot.fill.name = gg.labs.name)
#+ plot.width=10, plot.height=7
print(group2.unhealthy.genera.plot+
  ggtitle(paste0("Relative abundance of Group 2 members associated with unhealthy aging")))

#' Group 2 unhealthy families
group2.unhealthy.families
group2.unhealthy.families.plot<-ggplot_species(taxa.to.plot=group2.unhealthy.families,
                                              tax.df=ps.q.agg.family.relab,
                                              tax.rank="Family",
                                              sample.order=sample.levels,
                                              group.names = pretty.level.names,
                                              grouping.variable = "agegroup",
                                              metadata.df = custom.md,
                                              ggplot.fill.name = gg.labs.name)
#+ plot.width=10, plot.height=7
print(group2.unhealthy.families.plot+
  ggtitle(paste0("Relative abundance of Group 2 members associated with unhealthy aging")))

#' Group 2 unhealthy species
group2.unhealthy.species
group2.unhealthy.species.plot<-ggplot_species(taxa.to.plot =group2.unhealthy.species,
                                             tax.df=ps.q.agg.relab,
                                              tax.rank="Species",
                                             sample.order=sample.levels,
                                             group.names = pretty.level.names,
                                             grouping.variable = "agegroup",
                                             metadata.df = custom.md,
                                             ggplot.fill.name = gg.labs.name)

#+ plot.width=10, plot.height=7
print(group2.unhealthy.species.plot+
  ggtitle(paste0("Relative abundance of Group 2 members associated with unhealthy aging")))

#' Group 3 genera
group3.genera
group3.genera.plot<-ggplot_species(taxa.to.plot =group3.genera,
                                    tax.df = ps.q.agg.genus.relab,
                                    tax.rank = "Genus",
                                    sample.order=sample.levels,
                                    group.names = pretty.level.names,
                                    grouping.variable = "agegroup",
                                    metadata.df = custom.md,
                                    ggplot.fill.name = gg.labs.name)
#+ plot.width=10, plot.height=7
print(group3.genera.plot+
  ggtitle(paste0("Relative abundance of Group 3 members")))

#' Group 3 genera
group3.families.plot<-ggplot_species(taxa.to.plot = group3.families,
                                      tax.df = ps.q.agg.family.relab,
                                    tax.rank = "Family",
                                    sample.order=sample.levels,
                                    group.names = pretty.level.names,
                                    grouping.variable = "agegroup",
                                    metadata.df = custom.md,
                                    ggplot.fill.name = gg.labs.name)
#+ plot.width=10, plot.height=7
print(group3.families.plot+
  ggtitle(paste0("Relative abundance of Group 3 members")))

for (image.format in c("png","tiff")){
  ggsave(paste0("./images/barplots/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "group1.genera-nmr",
                      sep = "-"),".",image.format),
         plot=group1.genera.plot,
         width=10, height=7,units="in",
         dpi=300,device = image.format) # 4 plots
  ggsave(paste0("./images/barplots/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "group2.increased.genera-nmr",
                      sep = "-"),".",image.format),
         plot=group2.increased.genera.plot,
         width=10, height=7,units="in",
         dpi=300,device = image.format) # 4 plots
  ggsave(paste0("./images/barplots/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "group2.unhealthy.genera-nmr",
                      sep = "-"),".",image.format),
         plot=group2.unhealthy.genera.plot,
         width=10, height=7,units="in",
         dpi=300,device = image.format) # 7 plots
  ggsave(paste0("./images/barplots/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "group2.unhealthy.families-nmr",
                      sep = "-"),".",image.format),
         plot=group2.unhealthy.families.plot,
         width=10, height=5,units="in",
         dpi=300,device = image.format)    # 2 plots
  ggsave(paste0("./images/barplots/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "group2.unhealthy.species-nmr",
                      sep = "-"),".",image.format),
         plot=group2.unhealthy.species.plot,
         width=10, height=8,units="in",
         dpi=300,device = image.format)    # 6 plots
  ggsave(paste0("./images/barplots/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "group3.genera-nmr",
                      sep = "-"),".",image.format),
         plot= group3.genera.plot ,
         width=10, height=8,units="in",
         dpi=300,device = image.format)    # 5 plots
  ggsave(paste0("./images/barplots/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "group3.families-nmr",
                      sep = "-"),".",image.format),
         plot= group3.families.plot ,
         width=7, height=5,units="in",
         dpi=300,device = image.format)    # 1 plot
  
}

## Analyse each kingdom ####
#' Number of bacterial species
ps.q.agg.dominant.species%>%
  filter(Kingdom=="Bacteria")%>%
  nrow
# Save the table of bacterial abundances
ps.q.agg.bacteria<-ps.q.agg.dominant.species%>%
  ungroup%>%
  filter(Kingdom=="Bacteria")
# write.table(ps.q.agg.bacteria,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "kraken2","bacteria.tsv",sep="-")),
#             row.names = F,sep = "\t")

# Non-bacterial species ####
### 1. Eukaryotes ####
#### 1.1. How many eukaryotes do we have ####
ps.q.agg.dominant.species.nonbact%>%
  filter(Kingdom=="Eukaryota")%>%
  distinct(Species)%>%
  tally

#### 1.2 Total abundance of eukaryotes in classified data ####
ps.q.agg.relab%>%
  ungroup%>%
  filter(Kingdom=="Eukaryota")%>%
  reframe(TotalKingdom=sum(Abundance)/TotalClass*100)%>%
  distinct

#### 1.3 Eukaryotic phyla ####
ps.q.agg.dominant.species.nonbact%>%
  filter(Kingdom=="Eukaryota")%>%
  group_by(Phylum)%>%
  distinct(Species,.keep_all = T)%>%
  tally

#### 1.4 Eukaryotic species that aren't Streptophyta ####
ps.q.agg.dominant.species.nonbact%>%
  filter(Kingdom=="Eukaryota",Phylum!="Streptophyta")

# Save the table of eukaryotic abundances
ps.q.agg.eukaryotes<-ps.q.agg.dominant.species.nonbact%>%
  filter(Kingdom=="Eukaryota")
# write.table(ps.q.agg.eukaryotes,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "kraken2","eukaryotes.tsv",sep="-")),
#             row.names = F,sep = "\t")

### 2. Archaea ####
#### 2.1. How many archaea do we have ####
ps.q.agg.dominant.species%>%
  filter(Kingdom=="Archaea")%>%
  distinct(Species)%>%
  tally

#### 2.2 Total abundance of archaea in classified data ####
ps.q.agg.relab%>%
  ungroup%>%
  filter(Kingdom=="Archaea")%>%
  reframe(TotalKingdom=sum(Abundance)/TotalClass*100)%>%
  distinct

#### 2.3 Archaea phyla ####
ps.q.agg.dominant.species.nonbact%>%
  filter(Kingdom=="Archaea")%>%
  group_by(Phylum)%>%
  distinct(Species,.keep_all = T)%>%
  tally

#### 2.4 Archaea species ###
ps.q.agg.dominant.species.nonbact%>%
  filter(Kingdom=="Archaea")

# Save the table of archaea abundances
ps.q.agg.archaea<-ps.q.agg.dominant.species.nonbact%>%
  filter(Kingdom=="Archaea")
# write.table(ps.q.agg.archaea,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "kraken2","archaea.tsv",sep="-")),
#             row.names = F,sep = "\t")

### 3. Viruses ####
#### 2.1. How many archaea do we have ####
ps.q.agg.dominant.species.nonbact%>%
  filter(Kingdom=="Viruses")%>%
  distinct(Species)%>%
  tally

#### 2.2 Total abundance of viruses in classified data ####
ps.q.agg.relab%>%
  ungroup%>%
  filter(Kingdom=="Viruses")%>%
  reframe(TotalKingdom=sum(Abundance)/TotalClass*100)%>%
  distinct

#### 2.3 Viral phyla ####
ps.q.agg.dominant.species.nonbact%>%
  filter(Kingdom=="Viruses")%>%
  group_by(Phylum)%>%
  distinct(Species,.keep_all = T)%>%
  tally

#### 2.4 Viral species ###
ps.q.agg.dominant.species.nonbact%>%
  filter(Kingdom=="Viruses")%>%
  group_by(Phylum)%>%
  distinct(Species,.keep_all = T)

# Save the table of viral abundances
ps.q.agg.viruses<-ps.q.agg.dominant.species.nonbact%>%
  filter(Kingdom=="Viruses")
# write.table(ps.q.agg.viruses,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "kraken2","viruses.tsv",sep="-")),
#             row.names = F,sep = "\t")
sessionInfo()
rm(list = ls(all=TRUE))
gc()