#' ---
#' output: 
#'   bookdown::html_document2:
#'      toc: true
#' ---
#' 
#' ```{r, setup 003-summary-stats-kraken2.R, include=FALSE}
#' knitr::opts_knit$set(root.dir = '/home/rakhimov/projects/metagenome')
#' ```
#' 
#+ echo=FALSE
# Analysing Kraken2 data with phyloseq. ####
#' # Analysing Kraken2 data with phyloseq. 
#' 
#+ echo=FALSE
## Introduction ####
#'
#' ## Introduction
#' In this script, we will explore the imported dataset from Kraken2 using 
#' phyloseq.
#' 
#' We will use the data from 002-phyloseq-kraken2.R script (ps.q.agg
#' agglomerated tables at phylum, family, genus, and Species level).

#+ echo=FALSE
## 1. Load necessary libraries and scripts. ####
#'
#' ## Load necessary libraries and scripts.
# install.packages(c("tidyverse","ggtext"))
library(phyloseq)
library(tidyverse)
library(ggtext)
library(vegan)
#+ echo=FALSE
## 2. Import the config file. ####
#' 
#' ## Import the config file.
source(here::here("config/R/config.R"))# config file with global variables
source(here::here("config/R/themes.R"))# config file with themes

#' Load necessary scripts.
source(file.path(util.functions.r,"get_mean_and_sd_relab_from_phyloseq.R"))
source(file.path(util.functions.r,"create_summary_stats.R"))
source(file.path(util.functions.r,"multi_rarefy_phyloseq.R"))
source(file.path(util.functions.r,"get_average_summary_stats.R"))

source(file.path(util.functions.r,"get_unclassified_summary_stats.R"))
source(file.path(util.functions.r,"get_average_unclassified_summary_stats.R"))
source(file.path(util.functions.r,"ggplot_species.R"))

#+ echo=FALSE
## 3. Import raw datasets as rds files. ####
#'
#' Import raw datasets as rds files.
ps.q <- readRDS(file = ps.q.raw.fname)
ps.q.rel <- readRDS(file = ps.q.rel.raw.fname)
custom.md <- readRDS(custom.md.path)%>%
  filter(sequencing_type =="Naked mole-rat whole metagenome sequencing")
gg.labs.name <- "Age group"
#+ echo=FALSE
## 4. Calculating summary statistics. ####
#'
#' ## Calculating summary statistics.
#' 
#' Load the NMR age metadata.
custom.md.ages<-readRDS(file.path(metadatadir,"custom.md.ages.rds"))%>%
  filter(sequencing_type == "Naked mole-rat whole metagenome sequencing")
pretty.level.names <-names(table(custom.md$agegroup))
names(pretty.level.names)<-names(table(custom.md$agegroup))

#' Setup sample levels for NMR for barplots.
sample.levels<-custom.md.ages%>%
  ungroup()%>%
  select(Sample,age)%>%
  arrange(age)%>%
  distinct()%>%
  mutate(NewSample=paste0(Sample," (",age,")"))%>%
  mutate(Sample=factor(Sample,levels=Sample),
         NewSample=factor(NewSample,levels=NewSample))


#' Agglomerate the `ps.q.filtered` object and convert it into a dataframe with `psmelt()` 
#' function. Then, save it in tab-separated format and as an rds object. 
melt_phyloseq_and_save <-function (ps, ps.index){
  ps.indexes <- c("Kingdom","Phylum","Class","Order", "Family", "Genus","Species")
  
  if(ps.index == "OTU"){
    ps.df <- ps %>%
      psmelt()
  }else{
    ps.df <- ps %>%
      tax_glom(ps.index, NArm=FALSE)%>%
      psmelt()
  }
  
  if(ps.index == "OTU"){
    ps.df <- ps.df %>%
      dplyr::select(all_of(c("OTU","Sample", "Abundance", "class", ps.indexes)))
  }else{
    selected.indexes <- ps.indexes[1:match(ps.index, ps.indexes)] # vector of taxonomic ranks to select
    ps.df <- ps.df %>%
      dplyr::select(all_of(c("Sample", "Abundance", "class", selected.indexes)))
  }
  # Save the tables in TSV format and as an RDS object
  ps.df.fname <- case_when(ps.index =="OTU" ~ ps.q.agg.asv.fname.no_ext,
                           ps.index =="Species" ~ ps.q.agg.species.fname.no_ext,
                           ps.index =="Genus" ~ ps.q.agg.genus.fname.no_ext,
                           ps.index =="Family" ~ps.q.agg.family.fname.no_ext,
                           ps.index =="Phylum" ~ ps.q.agg.phylum.fname.no_ext
  ) 
  tsv.fname <- file.path(community.composition.tables, paste0(ps.df.fname,".tsv"))
  rda.fname <- file.path(community.composition.rdafiles, paste0(ps.df.fname,".rds"))
  if(!file.exists(tsv.fname)){
    write.table(ps.df,
                file = tsv.fname,
                row.names = F,sep = "\t")
  }
  if(!file.exists(rda.fname)){
    saveRDS(ps.df,
            file = rda.fname)
  }
  return(ps.df)
}

ps.q.agg.asv<-melt_phyloseq_and_save(ps.q, ps.index = "OTU")%>%
  filter(Abundance!=0)%>%
  as_tibble()
ps.q.rel.agg.asv<- ps.q.rel %>%
  psmelt() %>%
  dplyr::select(all_of(c("OTU","Sample", "Abundance", "class", 
                         c("Kingdom","Phylum","Class","Order", "Family", "Genus","Species"))))%>%
  filter(Abundance!=0)%>%
  as_tibble()%>%
  rename("RelativeAbundance" = "Abundance")
ps.q.agg.species<-melt_phyloseq_and_save(ps.q, ps.index = "Species")%>%
  filter(Abundance!=0)%>%
  as_tibble()
ps.q.agg.genus <- melt_phyloseq_and_save(ps.q, ps.index = "Genus")%>%
  filter(Abundance!=0)%>%
  as_tibble()
ps.q.agg.family <- melt_phyloseq_and_save(ps.q, ps.index = "Family")%>%
  filter(Abundance!=0)%>%
  as_tibble()
ps.q.agg.phylum<-melt_phyloseq_and_save(ps.q, ps.index = "Phylum")%>%
  filter(Abundance!=0)%>%
  as_tibble()


#' Save filtered phyloseq objects and the contents (OTU table, taxonomy, tree).
if(!file.exists(ps.q.raw.fname)){
  saveRDS(ps.q,
          file = ps.q.raw.fname)
}

if(!file.exists(ps.q.rel.raw.fname)){
  saveRDS(ps.q.rel,
          file = ps.q.rel.raw.fname)
}

if(!file.exists(ps.q.raw.otu_table.fname)){
  saveRDS(ps.q@otu_table %>% as.matrix() %>%as.data.frame(),
          file = ps.q.raw.otu_table.fname)
}
if(!file.exists(ps.q.raw.tax_table.fname)){
  saveRDS(ps.q@tax_table%>% as.matrix() %>%as.data.frame(),
          file = ps.q.raw.tax_table.fname)
}
ps.q.raw.otu_table.fname.tsv <- 
  file.path(community.composition.tables, 
            gsub("rds", "tsv",basename(ps.q.raw.otu_table.fname))
)
if(!file.exists(ps.q.raw.otu_table.fname.tsv)){
  write.table(ps.q@otu_table%>% as.matrix() %>%as.data.frame() %>% rownames_to_column("OTU"),
              file = ps.q.raw.otu_table.fname.tsv,
              sep = "\t", 
              row.names = FALSE)
  
}
ps.q.raw.tax_table.fname.tsv <- 
  file.path(community.composition.tables, 
            gsub("rds", "tsv",basename(ps.q.raw.tax_table.fname)))
if(!file.exists(ps.q.raw.tax_table.fname.tsv)){
  write.table(ps.q@tax_table%>% as.matrix() %>%as.data.frame()%>% rownames_to_column("OTU"),
              file = ps.q.raw.tax_table.fname.tsv,
              sep = "\t", 
              row.names = FALSE)
  
}


#+ echo=FALSE
### 4.1 Calculating summary statistics per host. ####
#'
#' ### Calculating summary statistics per host.
#' 
#+ echo=FALSE
#' The columns are:  
#' - Total reads  
#' - Library size (mean Abundance ± SD)  
#' - Number of ASV per host  
#' - Number of phyla per host  
#' - Number of families per host  
#' - Number of genera per host  
#' Here, Species are called ASV because the same function is used in the 
#' 16S rRNA gene sequencing analysis.  
summary.stats <- create_summary_stats(ps.q.agg.asv)
knitr::kable(as.data.frame(summary.stats),format = "simple")
summary.stats.fname<-file.path(community.composition.tables,
                               "summary-table.tsv")
if(!file.exists(summary.stats.fname)){
  write.table(summary.stats,
              file=summary.stats.fname,
              row.names = F,sep = "\t")
}

#+ echo=FALSE
### 4.2 Create a summary table for each sample. ####
#'
#' ### Create a summary table for each sample.
n.taxa.per_sample<-ps.q.agg.asv%>%
  group_by(Sample)%>%
  summarise(n_species=n_distinct(Species),
            n_genera=n_distinct(Genus),
            n_families=n_distinct(Family),
            n_phyla=n_distinct(Phylum))
knitr::kable(n.taxa.per_sample,format = "simple")
n.taxa.per_sample.fname<-file.path(community.composition.tables,
                               "kraken2-n_taxa-per-sample-table.tsv")
if(!file.exists(n.taxa.per_sample.fname)){
  write.table(n.taxa.per_sample,
              file=n.taxa.per_sample.fname,
              row.names = F,sep = "\t")
}

#+ echo=FALSE
## 5. Calculate average, min, max, and SD relative abundances from phyloseq object. ####
#'
#' ## Calculate average, min, max, and SD relative abundances from phyloseq object.
#' For phyla, families, and genera, no grouping. For species 
#' (here, we call them OTU), create a dataframe by age group and total (no grouping).
mean_sd_relab.nmr_no_groups.phyla <- get_mean_and_sd_relab_from_phyloseq(ps.q.rel, "class", "Phylum")
mean_sd_relab.nmr_no_groups.families <- get_mean_and_sd_relab_from_phyloseq(ps.q.rel, "class", "Family")
mean_sd_relab.nmr_no_groups.genera <- get_mean_and_sd_relab_from_phyloseq(ps.q.rel, "class", "Genus")
#' By age group.
mean_sd_relab.nmr_ages.asv <- 
  get_mean_and_sd_relab_from_phyloseq(ps.q.rel, "agegroup", "OTU")
#' In total (no grouping).
mean_sd_relab.nmr_no_groups.asv <- 
  get_mean_and_sd_relab_from_phyloseq(ps.q.rel, "class", "OTU")

knitr::kable(head(mean_sd_relab.nmr_no_groups.phyla),format = "simple")
knitr::kable(head(mean_sd_relab.nmr_no_groups.families),format = "simple")
knitr::kable(head(mean_sd_relab.nmr_no_groups.genera),format = "simple")
knitr::kable(head(mean_sd_relab.nmr_ages.asv),format = "simple")
knitr::kable(head(mean_sd_relab.nmr_no_groups.asv),format = "simple")

#' Save tables.
mean_sd_relab.nmr_no_groups.phyla.fname <- 
  file.path(community.composition.tables,
            "mean-max-min-sd-relab-nmr_no_groups-Phylum.tsv")
if(!file.exists(mean_sd_relab.nmr_no_groups.phyla.fname)){
  write.table(mean_sd_relab.nmr_no_groups.phyla,
              file = mean_sd_relab.nmr_no_groups.phyla.fname,
              row.names = F,sep = "\t")
}

mean_sd_relab.nmr_no_groups.families.fname <-
  file.path(community.composition.tables,
            "mean-max-min-sd-relab-nmr_no_groups-Family.tsv")
if(!file.exists (mean_sd_relab.nmr_no_groups.families.fname)){
  write.table(mean_sd_relab.nmr_no_groups.families,
              file = mean_sd_relab.nmr_no_groups.families.fname,
              row.names = F,sep = "\t")
}

mean_sd_relab.nmr_no_groups.genera.fname <-
  file.path(community.composition.tables,
            "mean-max-min-sd-relab-nmr_no_groups-Genus.tsv")
if(! file.exists (mean_sd_relab.nmr_no_groups.genera.fname)){
  write.table(mean_sd_relab.nmr_no_groups.genera,
              file = mean_sd_relab.nmr_no_groups.genera.fname,
              row.names = F,sep = "\t")
}

mean_sd_relab.nmr_ages.asv.fname <-
  file.path(community.composition.tables,
            "mean-max-min-sd-relab-nmr_ages-Species.tsv")
if(! file.exists (mean_sd_relab.nmr_ages.asv.fname)){
  write.table(mean_sd_relab.nmr_ages.asv,
              file = mean_sd_relab.nmr_ages.asv.fname,
              row.names = F,sep = "\t")
}
mean_sd_relab.nmr_no_groups.asv.fname <-
  file.path(community.composition.tables,
            "mean-max-min-sd-relab-nmr_no_groups-Species.tsv")
if(! file.exists (mean_sd_relab.nmr_no_groups.asv.fname)){
  write.table(mean_sd_relab.nmr_no_groups.asv,
              file = mean_sd_relab.nmr_no_groups.asv.fname,
              row.names = F,sep = "\t")
}

#+ echo=FALSE
## 6. Rarefy the table and check the percentage of unclassified taxa. ####
#'
#' ## Rarefy the table and check the percentage of unclassified taxa. ####
#' First, perform multiple rarefaction:
# print(paste("Rarefying data with", rare.num_samples, "subsamples"))
# rarefied.fname.all <- file.path(community.composition.rdafiles,paste0(
#   paste("ps.q.rare","OTU",paste0(rare.num_samples,"_iter"),
#         paste("NMR",collapse = '-'),sep = "-"),".rds"))
# if(!file.exists(rarefied.fname.all)){
#   ps.q.agg.asv.nmr.rare<-multi_rarefy_phyloseq(ps.q,
#                                                niter = rare.num_samples,
#                                                output.filename.rds =  rarefied.fname.all
#                                                )
# }else{
#   ps.q.agg.asv.nmr.rare<- readRDS(rarefied.fname.all)
# }

#' #+ echo=FALSE
#' ### 6.3 Create a summary stats table for the rarefied dataframe. ####
#' #' 
#' #' ### Create a summary stats table for the rarefied dataframe.
#' summary.stats.rare<-get_average_summary_stats(rare.list.asv = ps.q.agg.asv.nmr.rare)
#' knitr::kable(summary.stats.rare, format ="simple")
#' summary.stats.rare.fname<-file.path(community.composition.tables,
#'                                     "summary-table-rarefied.tsv")
#' if(!file.exists(summary.stats.rare.fname)){
#'   write.table(summary.stats.rare,
#'               file=summary.stats.rare.fname,
#'               row.names = F,sep = "\t")
#' }

#+ echo=FALSE
## 8. Check the most abundant bacterial phyla, families, genera in NMR. ####
#'
## Check the most abundant bacterial phyla, families, genera in NMR.
#' Phyla:
mean_sd_relab.nmr_no_groups.phyla %>%
  filter(Kingdom == "Bacteria")%>%
  head()%>%
  knitr::kable(format = "simple")

#' Families:
mean_sd_relab.nmr_no_groups.families %>%
  filter(Kingdom == "Bacteria")%>%
  head()%>%
  knitr::kable(format = "simple")

#' Genera:
mean_sd_relab.nmr_no_groups.genera %>%
  filter(Kingdom == "Bacteria")%>%
  head()%>%
  knitr::kable(format = "simple")


#' Species:
mean_sd_relab.nmr_no_groups.asv %>%
  filter(Kingdom == "Bacteria")%>%
  head()%>%
  knitr::kable(format = "simple")

#+ echo=FALSE
## 8.1 Check the most abundant non-bacterial phyla, families, genera in NMR. ####
#'
#' ## Check the most abundant non-bacterial phyla, families, genera in NMR.
#' Non-bacterial phyla:
mean_sd_relab.nmr_no_groups.phyla.nonbact<-mean_sd_relab.nmr_no_groups.phyla%>%
  filter(Kingdom != "Bacteria")
head(mean_sd_relab.nmr_no_groups.phyla.nonbact)%>%
  knitr::kable(format = "simple")
#' Non-bacterial families:
mean_sd_relab.nmr_no_groups.families.nonbact<-mean_sd_relab.nmr_no_groups.families%>%
  filter(Kingdom != "Bacteria")
head(mean_sd_relab.nmr_no_groups.families.nonbact)%>%
  knitr::kable(format = "simple")
#' Non-bacterial genera:
mean_sd_relab.nmr_no_groups.genera.nonbact<-mean_sd_relab.nmr_no_groups.genera%>%
  filter(Kingdom != "Bacteria")
head(mean_sd_relab.nmr_no_groups.genera.nonbact)%>%
  knitr::kable(format = "simple")
#' Non-bacterial species:
mean_sd_relab.nmr_no_groups.asv.nonbact<-mean_sd_relab.nmr_no_groups.asv%>%
  filter(Kingdom != "Bacteria")
head(mean_sd_relab.nmr_no_groups.asv.nonbact)%>%
  knitr::kable(format = "simple")

mean_sd_relab.nmr_no_groups.phyla.nonbact.fname <-file.path(community.composition.tables,
                                                            "mean-max-min-sd-relab-nmr_no_groups-Phylum-nonbact.tsv")

mean_sd_relab.nmr_no_groups.families.nonbact.fname <-file.path(community.composition.tables,
                                                               "mean-max-min-sd-relab-nmr_no_groups-Family-nonbact.tsv")

mean_sd_relab.nmr_no_groups.genera.nonbact.fname <- file.path(community.composition.tables,
                                                              "mean-max-min-sd-relab-nmr_no_groups-Genus-nonbact.tsv")

mean_sd_relab.nmr_no_groups.asv.nonbact.fname <- file.path(community.composition.tables,
                                                           "mean-max-min-sd-relab-nmr_no_groups-OTU-nonbact.tsv")

if(!file.exists(mean_sd_relab.nmr_no_groups.phyla.nonbact.fname)){
  write.table(mean_sd_relab.nmr_no_groups.phyla.nonbact,
              file = mean_sd_relab.nmr_no_groups.phyla.nonbact.fname,
              row.names = F,sep = "\t")
}
if(!file.exists(mean_sd_relab.nmr_no_groups.families.nonbact.fname)){
  write.table(mean_sd_relab.nmr_no_groups.families.nonbact,
              file = mean_sd_relab.nmr_no_groups.families.nonbact.fname,
              row.names = F,sep = "\t")
}
if(!file.exists(mean_sd_relab.nmr_no_groups.genera.nonbact.fname)){
  write.table(mean_sd_relab.nmr_no_groups.genera.nonbact,
              file = mean_sd_relab.nmr_no_groups.genera.nonbact.fname,
              row.names = F,sep = "\t")
}
if(!file.exists(mean_sd_relab.nmr_no_groups.asv.nonbact.fname)){
  write.table(mean_sd_relab.nmr_no_groups.asv.nonbact,
              file = mean_sd_relab.nmr_no_groups.asv.nonbact.fname,
              row.names = F,sep = "\t")
}

#+ echo=FALSE
## 9. Check how much Bacteroidaceae are in NMR.  ####
#'
#' ## Check how much Bacteroidaceae are in NMR.
bacteroidaceae.nmr<-mean_sd_relab.nmr_no_groups.families%>%
  filter(Family=="Bacteroidaceae")
knitr::kable(bacteroidaceae.nmr, format ="simple")
bacteroidaceae.nmr.fname <-file.path(community.composition.tables,
                                  paste("kraken2", 
                                        "bacteroidaceae-nmr-table.tsv",sep="-"))
if(! file.exists (bacteroidaceae.nmr.fname)){
  write.table(bacteroidaceae.nmr,
              file = bacteroidaceae.nmr.fname,
              row.names = F,sep = "\t")
}

#+ echo=FALSE
### 9.1 Check the most dominant Bacteroidota families in NMR. ####
#'
#'### Check the most dominant Bacteroidota families in NMR.
bacteroidota.nmr<-mean_sd_relab.nmr_no_groups.families%>%
  filter(Phylum=="Bacteroidota")
knitr::kable(bacteroidota.nmr, format ="simple")
bacteroidota.nmr.fname <-file.path(community.composition.tables,
                                     paste("kraken2", 
                                           "bacteroidota-nmr-table.tsv",sep="-"))
if(! file.exists (bacteroidota.nmr.fname)){
  write.table(bacteroidota.nmr,
              file = bacteroidota.nmr.fname,
              row.names = F,sep = "\t")
}


#+ echo=FALSE
## 10. Check Spirochaetaceae, Spirochaetota, and Treponema in NMR. ####
#'
#' ## Check Spirochaetaceae, Spirochaetota, and Treponema in NMR. 
spirochaetaceae.nmr<-mean_sd_relab.nmr_no_groups.families%>%
  filter(Family=="Spirochaetaceae")
knitr::kable(spirochaetaceae.nmr, format ="simple")

spirochaetota.nmr<-mean_sd_relab.nmr_no_groups.phyla%>%
  filter(Phylum=="Spirochaetota")
knitr::kable(spirochaetota.nmr, format ="simple")

spirochaetaceae.nmr.fname <- file.path(community.composition.tables,
                                       paste("kraken2", 
                                             "spirochaetaceae-nmr-table.tsv",sep="-") )
spirochaetota.nmr.fname <- file.path(community.composition.tables,
                                       paste("kraken2", 
                                             "spirochaetota-nmr-table.tsv",sep="-") )

if(!file.exists(spirochaetaceae.nmr.fname)){
  write.table(spirochaetaceae.nmr,
              file = spirochaetaceae.nmr.fname,
              row.names = F,sep = "\t")
}
if(!file.exists(spirochaetota.nmr.fname)){
  write.table(spirochaetota.nmr,
              file = spirochaetota.nmr.fname,
              row.names = F,sep = "\t")
}

#' Check the number of Species in Treponema from NMR.
mean_sd_relab.nmr_no_groups.asv%>%
  ungroup()%>%
  filter(Genus=="Treponema",class=="NMR")%>%
  distinct(OTU)%>%
  tally%>%
  knitr::kable(format = "simple")


treponema.nmr<-mean_sd_relab.nmr_no_groups.genera%>%
  filter(Genus=="Treponema")
knitr::kable(treponema.nmr, format ="simple")

treponema.nmr.names<-mean_sd_relab.nmr_no_groups.asv%>%
  filter(Genus=="Treponema")
knitr::kable(treponema.nmr.names, format ="simple")

treponema.nmr.fname <- file.path(community.composition.tables,
                                 paste("kraken2", 
                                       "treponema-nmr-table.tsv",sep="-") )

treponema.nmr.names.fname <- file.path(community.composition.tables,
                                       paste("kraken2", 
                                             "treponema-nmr-names.tsv",sep="-") )
if(!file.exists(treponema.nmr.fname)){
  write.table(treponema.nmr,
              file = treponema.nmr.fname,
              row.names = F,sep = "\t")
}
if(!file.exists(treponema.nmr.names.fname)){
  write.table(treponema.nmr.names,
              file = treponema.nmr.names.fname,
              row.names = F,sep = "\t")
}

#+ echo=FALSE
## 11. Analyse sulfur-metabolising bacteria in NMR. ####
#'
#' ## Analyse sulfur-metabolising bacteria in NMR. 
#' We don't know where the "sulf" pattern is (Phylum, Order, Genus, or anywhere 
#' else), so we use sapply with grepl.
#' https://stackoverflow.com/questions/47941680/grepl-across-multiple-specified-columns
all.ranks<-c("Kingdom", "Phylum", "Class", "Order", "Family","Genus","Species")
sulfur.bact.nmr<-
  ps.q.agg.species[!!rowSums(sapply(ps.q.agg.species[,all.ranks], grepl, 
                            pattern = "sulf|thio") ),]%>%
  distinct(Species,.keep_all = T)%>%
  select(all_of(all.ranks))%>%
  inner_join(mean_sd_relab.nmr_no_groups.asv)
head(sulfur.bact.nmr)%>%
  knitr::kable(format = "simple")
sulfur.bact.nmr.fname <- file.path(community.composition.tables,
                                   paste("kraken2", 
                                         "sulfur.bact.nmr-nmr-table.tsv",sep="-") )
if(!file.exists(sulfur.bact.nmr.fname)){
  write.table(sulfur.bact.nmr,
              file = sulfur.bact.nmr.fname,
              row.names = F,sep = "\t")
}

#+ echo=FALSE
### 11.2 Plot sulfur metabolising bacteria. ####
#'
#' ### Plot sulfur metabolising bacteria.
sulfur.bact.plot <- ps.q.rel.agg.asv%>%
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
#+ plot.width=11, plot.height=8
print(sulfur.bact.plot +
        theme(axis.title.y = element_text(size = 5),
              axis.title = element_text(size = 5),
              axis.text.y = ggtext::element_markdown(size=5),
              axis.text.x = element_text(size=6, angle=45, hjust=1),
              strip.text.x = ggtext::element_markdown(size=6),
              plot.title = element_text(size = 10),
              legend.text = element_text(size = 8),
              legend.title = element_text(size = 10),
              legend.position = "right",
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank())+
        ggtitle(paste0("Relative abundance of sulfur-utilizing bacteria in different naked mole-rat age groups")))
for (image.format in c("png","tiff")){
  ggsave(filename = paste0("sulfur-bacteria-nmr.", image.format),
         path = community.composition.figures,
         plot=sulfur.bact.plot,
         width=11, height=8,units="in",
         dpi=300,device = image.format)
}

treponema.plot<-ps.q.rel.agg.asv%>%
  filter(Genus=="Treponema")%>%
  arrange(Species)%>%
  mutate(Species = factor(Species, levels = unique(Species)),
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
#+ plot.width=11, plot.height=8
print(treponema.plot+
        theme(axis.title.y = element_text(size = 5),
              axis.title = element_text(size = 5),
              axis.text.y = ggtext::element_markdown(size=5),
              axis.text.x = element_text(size=6, angle=45, hjust=1),
              strip.text.x = ggtext::element_markdown(size=8),
              plot.title = element_text(size = 10),
              legend.text = element_text(size = 8),
              legend.title = element_text(size = 10),
              legend.position = "right",
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank())+
        ggtitle(paste0("Relative abundance of Treponemas in naked mole-rats")))
for (image.format in c("png","tiff")){
  ggsave(filename = paste0("treponema.sp-nmr.", image.format),
         path = community.composition.figures,
         plot=treponema.plot,
         width=11, height=8,units="in",
         dpi=300,device = image.format)
}

#+ echo=FALSE
## 12. Groups 1, 2, and 3 (**Supplementary figure**). ####
#'
#' ## Groups 1, 2, and 3 (**Supplementary figure**).
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

#' Abundances of group 1 genera (decreased with age and associated with 
#' healthy aging in humans):
group1.genera
group1.genera.plot<-ggplot_species(group1.genera,
                                   ps.q.rel.agg.asv,
                                   tax.rank="Genus",
                                   sample.order=sample.levels,
                                   group.names = pretty.level.names,
                                   grouping.variable = "agegroup",
                                   metadata.df = custom.md,
                                   ggplot.fill.name = gg.labs.name)
#+ plot.width=10, plot.height=7
print(group1.genera.plot+
        theme(axis.text.x = element_text(size=6, angle=45, hjust=1))+
        ggtitle(paste0("Relative abundance of Group 1 members"))) # 4 plots

group1.species%in%ps.q.rel.agg.asv$Species # no need to plot

#' Abundances of group 2 genera (increased with age):
group2.increased.genera
group2.increased.genera.plot<-ggplot_species(taxa.to.plot =group2.increased.genera,
                                             tax.df =ps.q.rel.agg.asv,
                                             tax.rank="Genus",
                                             sample.order=sample.levels,
                                             group.names = pretty.level.names,
                                             grouping.variable = "agegroup",
                                             metadata.df = custom.md,
                                             ggplot.fill.name = gg.labs.name)
#+ plot.width=10, plot.height=7
print(group2.increased.genera.plot+
        theme(plot.title = element_text(size = 13),
              axis.text.x = element_text(size=6, angle=45, hjust=1))+
        ggtitle(paste0("Relative abundance of Group 2 members increased with age")))

#' Abundances of group 2 genera (associated with unhealthy aging):
group2.unhealthy.genera
group2.unhealthy.genera.plot<-ggplot_species(taxa.to.plot=group2.unhealthy.genera,
                                              tax.df=ps.q.rel.agg.asv,
                                              tax.rank="Genus",
                                              sample.order=sample.levels,
                                              group.names = pretty.level.names,
                                              grouping.variable = "agegroup",
                                              metadata.df = custom.md,
                                              ggplot.fill.name = gg.labs.name)
#+ plot.width=10, plot.height=10
print(group2.unhealthy.genera.plot+
        theme(plot.title = element_text(size = 13),
              axis.text.x = element_text(size=6, angle=45, hjust=1))+
  ggtitle(paste0("Relative abundance of Group 2 members associated with unhealthy aging")))

#' Abundances of group 2 families (associated with unhealthy aging):
group2.unhealthy.families
group2.unhealthy.families.plot<-ggplot_species(taxa.to.plot=group2.unhealthy.families,
                                              tax.df=ps.q.rel.agg.asv,
                                              tax.rank="Family",
                                              sample.order=sample.levels,
                                              group.names = pretty.level.names,
                                              grouping.variable = "agegroup",
                                              metadata.df = custom.md,
                                              ggplot.fill.name = gg.labs.name)
#+ plot.width=10, plot.height=7
print(group2.unhealthy.families.plot+
        theme(plot.title = element_text(size = 13),
              axis.text.x = element_text(size=6, angle=45, hjust=1))+
        ggtitle(paste0("Relative abundance of Group 2 members associated with unhealthy aging")))

#' Abundances of group 2 species (associated with unhealthy aging):
group2.unhealthy.species
group2.unhealthy.species.plot<-ggplot_species(taxa.to.plot =group2.unhealthy.species,
                                             tax.df = ps.q.rel.agg.asv,
                                              tax.rank="Species",
                                             sample.order=sample.levels,
                                             group.names = pretty.level.names,
                                             grouping.variable = "agegroup",
                                             metadata.df = custom.md,
                                             ggplot.fill.name = gg.labs.name)

#+ plot.width=10, plot.height=7
print(group2.unhealthy.species.plot+
        theme(plot.title = element_text(size = 13),
              axis.text.x = element_text(size=6, angle=45, hjust=1))+
        ggtitle(paste0("Relative abundance of Group 2 members associated with unhealthy aging")))

#' Abundances of group 3 genera (increased with age and associated with 
#' healthy aging):
group3.genera
group3.genera.plot<-ggplot_species(taxa.to.plot =group3.genera,
                                    tax.df = ps.q.rel.agg.asv,
                                    tax.rank = "Genus",
                                    sample.order=sample.levels,
                                    group.names = pretty.level.names,
                                    grouping.variable = "agegroup",
                                    metadata.df = custom.md,
                                    ggplot.fill.name = gg.labs.name)
#+ plot.width=10, plot.height=7
print(group3.genera.plot+
        theme(axis.text.x = element_text(size=6, angle=45, hjust=1))+
        ggtitle(paste0("Relative abundance of Group 3 members")))

#' Abundances of group 3 families (increased with age and associated 
#' with healthy aging):
group3.families.plot<-ggplot_species(taxa.to.plot = group3.families,
                                      tax.df = ps.q.rel.agg.asv,
                                    tax.rank = "Family",
                                    sample.order=sample.levels,
                                    group.names = pretty.level.names,
                                    grouping.variable = "agegroup",
                                    metadata.df = custom.md,
                                    ggplot.fill.name = gg.labs.name)
#+ plot.width=10, plot.height=7
print(group3.families.plot+
        theme(axis.text.x = element_text(size=6, angle=45, hjust=1))+
        ggtitle(paste0("Relative abundance of Group 3 members")))

for (image.format in c("png","tiff")){
  ggsave(filename = paste0("group1.genera-nmr.", image.format),
         path = community.composition.figures,
         plot=group1.genera.plot,
         width=10, height=7,units="in",
         dpi=300,device = image.format) # 4 plots
  ggsave(filename = paste0("group2.increased.genera-nmr.",image.format),
         path = community.composition.figures,
         plot=group2.increased.genera.plot,
         width=10, height=7,units="in",
         dpi=300,device = image.format) # 4 plots
  ggsave(filename = paste0("group2.unhealthy.genera-nmr.",image.format),
         path = community.composition.figures,
         plot=group2.unhealthy.genera.plot,
         width=10, height=7,units="in",
         dpi=300,device = image.format) # 7 plots
  ggsave(filename = paste0("group2.unhealthy.families-nmr.",image.format),
         path = community.composition.figures,
         plot = group2.unhealthy.families.plot,
         width=10, height=5,units="in",
         dpi=300,device = image.format)    # 2 plots
  ggsave(filename = paste0("group2.unhealthy.species-nmr.",image.format),
         path = community.composition.figures,
         plot = group2.unhealthy.species.plot,
         width=10, height=8,units="in",
         dpi=300,device = image.format)    # 6 plots
  ggsave(paste0("group3.genera-nmr.", image.format), 
         path = community.composition.figures,
         plot= group3.genera.plot ,
         width=10, height=8,units="in",
         dpi=300,device = image.format)    # 5 plots
  ggsave(paste0("group3.families-nmr.", image.format),
         path = community.composition.figures,
         plot = group3.families.plot ,
         width=7, height=5,units="in",
         dpi=300,device = image.format)    # 1 plot

}

#+ echo=FALSE
## 13. Analyse each kingdom in detail. ####
#'
#' ## Analyse each kingdom in detail. 
#'
#' Number of bacterial species
mean_sd_relab.nmr_no_groups.asv%>%
  filter(Kingdom=="Bacteria")%>%
  nrow

#+ echo=FALSE
### 13.1 Eukaryotes ####
#'
#' ### Eukaryotes
#' How many eukaryote species do we have
mean_sd_relab.nmr_no_groups.asv.nonbact%>%
  filter(Kingdom=="Eukaryota")%>%
  distinct(Species)%>%
  tally()
#' Total abundance of eukaryotes in classified data
mean_sd_relab.nmr_no_groups.asv.nonbact%>%
  ungroup%>%
  filter(Kingdom=="Eukaryota")%>%
  distinct() %>%
  knitr::kable(format = "simple")

#' Eukaryotic phyla
mean_sd_relab.nmr_no_groups.asv.nonbact%>%
  filter(Kingdom=="Eukaryota")%>%
  group_by(Phylum)%>%
  distinct(Species,.keep_all = T)%>%
  tally()%>%
  knitr::kable(format = "simple")
#' Eukaryotic species that aren't Streptophyta
mean_sd_relab.nmr_no_groups.asv.nonbact%>%
  filter(Kingdom=="Eukaryota",Phylum!="Streptophyta")%>%
  knitr::kable(format = "simple")

#+ echo=FALSE
### 13.2 Archaea ####
#'
#' ### 13.2 Archaea
#' How many archaea do we have 
mean_sd_relab.nmr_no_groups.asv%>%
  filter(Kingdom=="Archaea")%>%
  distinct(Species)%>%
  tally()
#' Total abundance of archaea in classified data
mean_sd_relab.nmr_no_groups.asv.nonbact%>%
  ungroup%>%
  filter(Kingdom=="Archaea")%>%
  distinct()%>%
  knitr::kable(format = "simple")
#' Archaea phyla
mean_sd_relab.nmr_no_groups.asv.nonbact%>%
  filter(Kingdom=="Archaea")%>%
  group_by(Phylum)%>%
  distinct(Species,.keep_all = T)%>%
  tally()%>%
  knitr::kable(format = "simple")

#+ echo=FALSE
### 13.3 Viruses ####
#'
#' ### 13.3 Viruses
#' How many viruses do we have
mean_sd_relab.nmr_no_groups.asv.nonbact%>%
  filter(Kingdom=="Viruses")%>%
  distinct(Species)%>%
  tally()
#' Total abundance of viruses in classified data
mean_sd_relab.nmr_no_groups.asv.nonbact%>%
  ungroup%>%
  filter(Kingdom=="Viruses")%>%
  distinct()%>%
  knitr::kable(format = "simple")
#' Viral phyla 
mean_sd_relab.nmr_no_groups.asv.nonbact%>%
  filter(Kingdom=="Viruses")%>%
  group_by(Phylum)%>%
  distinct(Species,.keep_all = T)%>%
  tally()%>%
  knitr::kable(format = "simple")

sessionInfo()
rm(list =setdiff(ls(all.names = TRUE), c("markdown.dir")))
gc()
