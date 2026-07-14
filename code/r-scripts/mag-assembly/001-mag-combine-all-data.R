#' ---
#' output: 
#'   bookdown::html_document2:
#'      toc: true
#' ---
#' 
#' ```{r, setup 001-mag-combine-all-data.R, include=FALSE}
#' knitr::opts_knit$set(root.dir = '/home/rakhimov/projects/metagenome')
#' ```
#' # Metagenome assembly analysis {-}
#+ echo=FALSE
# Combine all the data from MAG assembly ####
#' # Combine all the data from MAG assembly
#' 
#+ echo=FALSE
## Introduction ####
#'
#' ## Introduction
#' Here, we combine BLASTN and GTDB-tk results (on all contigs from MEGAHIT) 
#' with seqkit stats (on contigs and dereplicated bins).
#+ echo=FALSE
## 1. Load necessary libraries. ####
#'
#' ## Load necessary libraries
# install.packages(c("tidyverse", "data.table"))
library(tidyverse)
library(data.table)

# load(file="./output/rdafiles/mag-combine-all-data.Rdata")
#+ echo=FALSE
## 2. Load data. ####
#'
#' ## Load data.
source(here::here("config/R/config.R"))# config file with global variables
source(here::here("config/R/themes.R"))# config file with themes
dir.create(file.path(mag.rdafiles),recursive = TRUE)
dir.create(file.path(mag.tables),recursive = TRUE)
dir.create(file.path(mag.figures),recursive = TRUE)

#+ echo=FALSE
## 3. Read the seqkit stats. ####
#'
#' ## Read the seqkit stats.
#' Seqkit on MEGAHIT [VERY LARGE]
seqkit.megahit.df<-data.table::fread(seqkit.megahit.fname,
                                     header = TRUE,
                              sep="\t")%>%
  as_tibble()%>%
  rename("sample_id" = "sample")
#' Seqkit on MetaBAT2
seqkit.metabat2.df<-read.table(seqkit.metabat2.fname,header = TRUE,
                              sep="\t",comment.char = "")%>%
  as_tibble()%>%
  rename("sample_id" = "sample")
#' Seqkit on dRep
seqkit.drep.df<-read.table(seqkit.drep.fname,header = TRUE,
                              sep="\t",comment.char = "")%>%
  as_tibble()%>%
  rename("sample_id" = "sample")
head(seqkit.megahit.df)%>%
  knitr::kable(format = "simple")%>%
  print()

head(seqkit.metabat2.df)%>%
  knitr::kable(format = "simple")%>%
  print()

head(seqkit.drep.df)%>%
  knitr::kable(format = "simple")%>%
  print()

#' 3118148 contigs from MEGAHIT
nrow(seqkit.megahit.df)
#' 156182 contigs from MetaBAT2
nrow(seqkit.metabat2.df)
#' 11456 contigs from dRep
nrow(seqkit.drep.df)
summary(seqkit.megahit.df)
summary(seqkit.metabat2.df)
summary(seqkit.drep.df)

#+ echo=FALSE
## 4. Add drep clusters information. ####
#'
#' ## Add drep clusters information.
#' We want to know which bins were assigned to which clusters, and which
#' bins weren't because they had low quality.
#' We need it to get a more precise estimate of mapped reads. 
#' Read the data with all high-quality MAGs and their clusters
drep.clusters<-read.csv(drep.clusters.fname)%>%
  as_tibble()%>%
  mutate(genome=gsub(paste0(metabat2.date_time,"_"),"",genome),
         genome=gsub(".fa","",genome))%>%
  separate_wider_delim(genome,delim = "_bin.",names = c("sample_id","bin_id"))%>%
  mutate(bin_id=as.integer(bin_id))
head(drep.clusters)%>%
  knitr::kable(format = "simple")%>%
  print()

#' Read the data with representative genomes and their clusters
drep.winning_genomes<-read.csv(drep.winning_genomes.fname)%>%
  as_tibble()%>%
  mutate(genome=gsub(paste0(metabat2.date_time,"_"),"",genome),
         genome=gsub(".fa","",genome))%>%
  separate_wider_delim(genome,delim = "_bin.",names = c("sample_id","bin_id"))%>%
  mutate(bin_id=as.integer(bin_id))
head(drep.winning_genomes)%>%
  knitr::kable(format = "simple")%>%
  print()

#' Combine drep MAGs with other checkm2 MAGs
drep.clusters.joined<-drep.winning_genomes%>%
  select(sample_id,bin_id,cluster)%>%
  rename("secondary_cluster" = "cluster",
         "drep_sample"= "sample_id",
         "drep_bin_id"= "bin_id")%>%
  full_join(drep.clusters[,c("sample_id","bin_id","secondary_cluster")],
            join_by("secondary_cluster"))%>%
  # Add a column that shows if a bin is representative or not.
  # Then, remove drep_sample and drep_bin_id because it's the same information
  mutate(is_representative_bin=ifelse(sample_id==drep_sample&bin_id==drep_bin_id,
                            TRUE,FALSE))%>%
  select(-drep_sample,-drep_bin_id)
head(drep.clusters.joined)%>%
  knitr::kable(format = "simple")%>%
  print()


rm(drep.clusters)
rm(drep.clusters.fname)
rm(drep.winning_genomes)
rm(drep.winning_genomes.fname)
#+ echo=FALSE
### 4.1 Join cluster information with seqkit stats. ####
#'
#' ### Join cluster information with seqkit stats.
seqkit.drep.df<-seqkit.drep.df%>%
  left_join(drep.clusters.joined,by = join_by(sample_id, bin_id))

seqkit.metabat2.df <- seqkit.metabat2.df%>%
  left_join(drep.clusters.joined,by = join_by(sample_id, bin_id))%>%
  mutate(secondary_cluster = ifelse(is.na(secondary_cluster), "low_quality_mag",
                                    secondary_cluster),
         is_representative_bin = ifelse(is.na(is_representative_bin), "low_quality_mag",
                                        is_representative_bin)
         )
  
head(seqkit.drep.df)%>%
  knitr::kable(format = "simple")%>%
  print()

head(seqkit.metabat2.df)

rm(drep.clusters.joined)
#+ echo=FALSE
## 5. Combine all seqkit statistics into one table. ####
#'
#' ## Combine all seqkit statistics into one table.
#' First, we create a column with each contig + its bin + sample_id.
seqkit.metabat2.df <- seqkit.metabat2.df%>%
  mutate(sample_contig=paste(sample_id,contig_id,sep = "_"),
         sample_bin_contig=paste(sample_id,bin_id,contig_id,sep = "_"))
seqkit.drep.df <- seqkit.drep.df%>%
  mutate(sample_contig=paste(sample_id,contig_id,sep = "_"),
         sample_bin_contig=paste(sample_id,bin_id,contig_id,sep = "_"))
seqkit.megahit.df <- seqkit.megahit.df%>%
  mutate(sample_contig=paste(sample_id,contig_id,sep = "_"))
head(seqkit.drep.df)%>%
  knitr::kable(format = "simple")%>%
  print()

head(seqkit.metabat2.df)%>%
  knitr::kable(format = "simple")%>%
  print()

head(seqkit.megahit.df)%>%
  knitr::kable(format = "simple")%>%
  print()

nrow(seqkit.drep.df)
nrow(seqkit.metabat2.df)
nrow(seqkit.megahit.df)

#' Get the common IDs between datasets to know whether a contig
#' was assembled (metabat2 or metabat2_drep) or not (megahit).
#' If it's metabat2, it is not representative. If it's metabat2_drep,
#' it is representative.
megahit.metabat2.common.ids<-intersect(seqkit.megahit.df$sample_contig,
                                      seqkit.metabat2.df$sample_contig)
metabat2.drep.common.ids<-intersect(seqkit.metabat2.df$sample_bin_contig,
                                   seqkit.drep.df$sample_bin_contig)
megahit.drep.common.ids<-intersect(seqkit.megahit.df$sample_contig,
                                   seqkit.drep.df$sample_contig)
length(megahit.metabat2.common.ids)
length(metabat2.drep.common.ids)
length(megahit.drep.common.ids)
#' All contigs from drep are found in MetaBAT2 (no differences):
length(setdiff(seqkit.drep.df$sample_bin_contig,
        seqkit.metabat2.df$sample_bin_contig))
#' All contigs from MetaBAT2 are found in MEGAHIT (no differences):
length(setdiff(seqkit.metabat2.df$sample_contig,
        seqkit.megahit.df$sample_contig))

#' Finally combine the data:
seqkit.all.df <- seqkit.megahit.df %>%
  full_join(subset(seqkit.metabat2.df,select = -c(AT,GC,sample_contig,
                                                      contig_length)),
            by = join_by(sample_id, contig_id))%>%
  relocate(sample_id,contig_id,bin_id)%>%
  # If bin_id, secondary_cluster, and is_representative_bin are NA, the contig was not assembled
  mutate(bin_id = ifelse(is.na(bin_id), "unassembled", bin_id),
         secondary_cluster = ifelse(is.na(secondary_cluster), "unassembled", secondary_cluster),
         is_representative_bin = ifelse(is.na(is_representative_bin), "unassembled", is_representative_bin))%>%
  select(-sample_contig,-sample_bin_contig)
head(seqkit.all.df)%>%
  knitr::kable(format = "simple")%>%
  print()

nrow(seqkit.all.df)
str(seqkit.all.df)
summary(seqkit.all.df)
rm(megahit.drep.common.ids)
rm(megahit.metabat2.common.ids)
rm(metabat2.drep.common.ids)

#' Barplot to count how many contigs were retained: most contigs weren't assembled
seqkit.all.df%>%
  mutate(assembly_status = ifelse(! secondary_cluster %in% c("low_quality_mag",
                                                               "unassembled"),
                                    "high_quality_mag", secondary_cluster))%>%
  count(sample_id,assembly_status )%>%
  # pivot_wider(names_from = "assembly_status",
  #             values_from="n")%>%
  ggplot(aes(x=sample_id,fill=assembly_status,y=n))+
  geom_bar(stat="identity")
rm(seqkit.drep.df)
rm(seqkit.megahit.df)
rm(seqkit.metabat2.df)
gc()

#+ echo=FALSE
## 6. Add taxonomy to the seqkit stats data. ####
#'
#' ## Add taxonomy to the seqkit stats data.
#+ echo=FALSE
### 6.1 Add GTDB-tk taxonomy (for dereplicated bins). ####
#'
#' ### Add GTDB-tk taxonomy (for dereplicated bins).
gtdbtk.taxonomy<-read.table(gtdbtk.taxonomy.fname,header = T,
                            sep="\t",comment.char = "")%>%
  as_tibble()%>%
  select(user_genome,classification) # remove unnecessary columns
head(gtdbtk.taxonomy)%>%
  knitr::kable(format = "simple")%>%
  print()

# 128 bins
nrow(gtdbtk.taxonomy)

#' Tidy up GTDB-tk taxonomy: substitute s__ and g__ with unclassified
#' Separate sample_bin column into separate columns
gtdbtk.taxonomy<-gtdbtk.taxonomy %>%
  mutate(classification=gsub(";g__;s__",";unclassified;unclassified",
                             classification))%>%
  mutate(classification=gsub(";s__$",";unclassified",classification))%>%
  mutate(classification=gsub("d__|p__|c__|o__|f__|g__|s__","",classification)) %>%
  mutate(user_genome=gsub(paste0(metabat2.date_time,"_"),"",user_genome))%>%
  separate_wider_delim(user_genome,delim="_bin.",names=c("sample_id","bin_id"))%>%
  mutate(bin_id=as.integer(bin_id))

#+ echo=FALSE
### 6.2 Add BLAST taxonomy (for all contigs from MEGAHIT). ####
#'
#' ### Add BLAST taxonomy (for all contigs from MEGAHIT).
blastn.taxonomy<-read.table(blastn.taxonomy.fname,header = F,
                            sep="\t",comment.char = "")%>%
  as_tibble()
head(blastn.taxonomy)%>%
  knitr::kable(format = "simple")%>%
  print()

# 1899 contigs
nrow(blastn.taxonomy)
colnames(blastn.taxonomy)<-c("sample_id","contig_id","silva_ref_id","perc_identity",
                             "alignment_length","bitscore","classification")

#+ echo=FALSE
### 6.3 Join stats with taxonomy. ####
#'
#' ### Join stats with taxonomy.

###NOT NEEDED?
# seqkit.with_taxonomy<-seqkit.all.df%>%
#   full_join(gtdbtk.taxonomy,
#             by = join_by(sample_id, bin_id))%>%
#   full_join(blastn.taxonomy[,c("sample_id","contig_id","classification")],
#             by = join_by(sample_id, contig_id))%>%
#   rename("gtdbtk_result"="classification.x",
#          "blastn_result"="classification.y")%>%
#   ### 6.4 Add column to show if a contig was classified, and which pipeline classified it ####
#   mutate(tax.classification.type=ifelse(!is.na(gtdbtk_result)&!is.na(blastn_result),
#                         "gtdbtk_and_blastn","unclassified"),
#          tax.classification.type=ifelse(!is.na(gtdbtk_result)&is.na(blastn_result),
#                         "gtdbtk_only",tax.classification.type),
#          tax.classification.type=ifelse(is.na(gtdbtk_result)&!is.na(blastn_result),
#                         "blastn_only",tax.classification.type))
# nrow(seqkit.with_taxonomy)
# head(seqkit.with_taxonomy)
# str(seqkit.with_taxonomy)
# summary(seqkit.with_taxonomy)
# rm(seqkit.all.df)

#+ echo=FALSE
## 7. Add bin and contig coverages from CoverM. ####
#'
#' ## Add bin and contig coverages from CoverM.
#+ echo=FALSE
### 7.1 CoverM data from dRep. ####
#'
#' ### CoverM data from dRep.
coverm.drep.df<-read.table(coverm.drep.fname,header = T,
                           sep = "\t")%>%
  as_tibble()%>%
  mutate(genome=gsub(paste0(metabat2.date_time,"_"),"",genome))
head(coverm.drep.df)%>%
  knitr::kable(format = "simple")%>%
  print()

nrow(coverm.drep.df)

#+ echo=FALSE
### 7.1 CoverM data from MetaBAT2. ####
#'
#' ### CoverM data from MetaBAT2.
coverm.metabat2.df<-read.table(coverm.metabat2.fname,header = T,
                               sep = "\t")%>%
  as_tibble()%>%
  mutate(genome=gsub(paste0(metabat2.date_time,"_"),"",genome))
head(coverm.metabat2.df)%>%
  knitr::kable(format = "simple")%>%
  print()

nrow(coverm.metabat2.df)

#+ echo=FALSE
### 7.1 CoverM data from MEGAHIT (use data.table package). ####
#'
#' ### CoverM data from MEGAHIT (use data.table package).
coverm.megahit.df<-data.table::fread(coverm.megahit.fname,header = T,
                                  sep = "\t")%>%
  as_tibble()%>%
  rename(sample_contig=contig)%>%
  separate_wider_delim(sample_contig,delim = "_",
                       names = c("sample_id","contig_id"),
                       too_many="merge" )
head(coverm.megahit.df)%>%
  knitr::kable(format = "simple")%>%
  print()

nrow(coverm.megahit.df)

#' How much relative abundance is classified in drep bins:
coverm.drep.df%>%
  mutate(sample_id=ifelse(grepl("unmapped",genome),
                       gsub("unmapped_","",genome),
                       gsub("_bin.[0-9]*","",genome)))%>%
  filter(!grepl("unmapped",genome))%>%
  group_by(sample_id)%>%
  summarise(sum_relative_ab=sum(relative_abundance))

#' How much relative abundance is classified in MetaBAT2 bins:
coverm.metabat2.df%>%
  mutate(sample_id=ifelse(grepl("unmapped",genome),
                       gsub("unmapped_","",genome),
                       gsub("_bin.[0-9]*","",genome)))%>%
  filter(!grepl("unmapped",genome))%>%
  group_by(sample_id)%>%
  summarise(sum_relative_ab=sum(relative_abundance))

#' How much relative abundance is classified in MetaBAT2 bins that belonged to
#' drep clusters:
# seqkit.with_taxonomy%>%
seqkit.all.df%>%
  filter(!secondary_cluster %in% c("unassembled", "low_quality_mag"))%>%
  mutate(genome=paste(sample_id,bin_id,sep="_bin."))%>%
  distinct(sample_id,bin_id,genome,secondary_cluster)%>%
  left_join(coverm.metabat2.df,by=join_by("genome"))%>%
  group_by(sample_id)%>%
  summarise(sum_relative_ab=sum(relative_abundance))

#+ echo=FALSE
## 8. Join relative abundances with taxonomy. ####
#'
#' ## Join relative abundances with taxonomy.
#+ echo=FALSE
### 8.1 GTDB-tk + dRep ####
#'
#' ### GTDB-tk + dRep
gtdbtk.coverm.drep <- seqkit.all.df%>% # seqkit.with_taxonomy%>%
  filter(!secondary_cluster %in% c("unassembled", "low_quality_mag"))%>%
  select(sample_id,bin_id,secondary_cluster)%>%
  distinct()%>%
  mutate(bin_id = as.integer(bin_id))%>% # conver bin_id into integer because we will join it with the gtdbtk data
  right_join(gtdbtk.taxonomy)%>%
  mutate(genome=paste(sample_id,bin_id,sep="_bin."))%>%
  left_join(coverm.drep.df[,c("genome","relative_abundance")],
            by=join_by("genome"))%>% # here, we add the relative abundances
  select(-genome)
head(gtdbtk.coverm.drep)%>%
  knitr::kable(format = "simple")%>%
  print()

#+ echo=FALSE
### 8.2 GTDB-tk + MetaBAT2 ####
#'
#' ### GTDB-tk + MetaBAT2
gtdbtk.coverm.metabat2 <- seqkit.all.df%>% # seqkit.with_taxonomy%>%
  filter(!secondary_cluster %in% c("unassembled", "low_quality_mag"))%>%
  distinct(sample_id,bin_id,secondary_cluster)%>%
  mutate(bin_id = as.integer(bin_id))%>%
  # full_join(gtdbtk.taxonomy,by = join_by(sample_id, bin_id)) %>%
  full_join(gtdbtk.coverm.drep[,c("secondary_cluster","classification")],
            by = join_by(secondary_cluster)) %>% # add the previous dataset because it has secondary_cluster info
  group_by(secondary_cluster) %>%
  mutate(classification = classification[!is.na(classification)][1L])%>% # remove NAs
  ungroup%>%
  mutate(genome=paste(sample_id,bin_id,sep="_bin."))%>%
  left_join(coverm.metabat2.df[,c("genome","relative_abundance")], 
            by=join_by("genome"))%>% # add relative abundances
  select(-genome)
head(gtdbtk.coverm.metabat2) %>%
  knitr::kable(format = "simple")%>%
  print()

#+ echo=FALSE
### 8.3 BLASTN + MEGAHIT ####
#'
#' ### BLASTN + MEGAHIT
blastn.coverm.megahit <- coverm.megahit.df%>%
  select(sample_id, contig_id,tpm)%>%
  inner_join(blastn.taxonomy,by = join_by(sample_id, contig_id))
head(blastn.coverm.megahit)%>%
  knitr::kable(format = "simple")%>%
  print()

# Remove CoverM data because it's already join with taxonomy
rm(coverm.drep.df)
rm(coverm.metabat2.df)
rm(coverm.megahit.df)
rm(blastn.taxonomy)

#' Save data
if (!file.exists(seqkit.all.df.fname.rds) & !file.exists(seqkit.all.df.fname.tsv)){
  saveRDS(seqkit.all.df, file = seqkit.all.df.fname.rds)
  data.table::fwrite(seqkit.all.df,
              file = seqkit.all.df.fname.tsv,
              col.names = TRUE,
              row.names = FALSE,
              quote = FALSE,
              sep = '\t',
              na =NA)
  
}

if (!file.exists(gtdbtk.coverm.drep.fname.rds) & !file.exists(gtdbtk.coverm.drep.fname.tsv)){
  saveRDS(gtdbtk.coverm.drep,file = gtdbtk.coverm.drep.fname.rds)
  write.table(gtdbtk.coverm.drep,
              file = gtdbtk.coverm.drep.fname.tsv,
              col.names = TRUE,
              row.names = FALSE,
              quote = FALSE,
              sep = '\t')
  
}
if (!file.exists(gtdbtk.coverm.metabat2.fname.rds) & !file.exists(gtdbtk.coverm.metabat2.fname.tsv)){
  saveRDS(gtdbtk.coverm.metabat2,file = gtdbtk.coverm.metabat2.fname.rds)
  write.table(gtdbtk.coverm.metabat2,
              file = gtdbtk.coverm.metabat2.fname.tsv,
              col.names = TRUE,
              row.names = FALSE,
              quote = FALSE,
              sep = '\t')
  
}
if (!file.exists(blastn.coverm.megahit.fname.rds) & !file.exists(blastn.coverm.megahit.fname.tsv)){
  saveRDS(blastn.coverm.megahit,file = blastn.coverm.megahit.fname.rds)
  write.table(blastn.coverm.megahit,
              file = blastn.coverm.megahit.fname.tsv,
              col.names = TRUE,
              row.names = FALSE,
              quote = FALSE,
              sep = '\t')
  
}

sessionInfo()
rm(list = ls(all=TRUE))
gc()
