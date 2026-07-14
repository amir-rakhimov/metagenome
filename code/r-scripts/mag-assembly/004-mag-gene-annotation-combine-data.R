#' ---
#' output: 
#'   bookdown::html_document2:
#'      toc: true
#' ---

#' ```{r, setup 004-mag-gene-annotation-combine-data.R, include=FALSE}
#' knitr::opts_knit$set(root.dir = '/home/rakhimov/projects/metagenome')
#' ```
#+ echo=FALSE
# Gene annotation analysis. ####
#' # Gene annotation analysis.
#'
#+ echo=FALSE
## Introduction ####
#'
#' ## Introduction
#' Import the data from gene annotation and convert into tibble format.
#' 
#+ echo=FALSE
## 1. Load necessary libraries. ####
#'
#' ## Load necessary libraries.
# install.packages(c("tidyverse", "data.table"))
library(tidyverse)
library(data.table)
#+ echo=FALSE
## 2. Import datasets. ####
#'
#' ## Import datasets.
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
#' Tidy taxonomy table
gtdbtk.taxonomy <- readRDS(gtdbtk.taxonomy.clean.fname.rds)

#+ echo=FALSE
## 3. Analyse PROKKA mapping between samples, contigs, and genes. ####
#'
#' ## Analyse PROKKA mapping between samples, contigs, and genes.
prokka.contig.gene.map.df <- data.table::fread(prokka.contig.gene.map.fname,
                                      header = FALSE,sep="\t")%>%
  as_tibble()
colnames(prokka.contig.gene.map.df)<-c("sample_id","contig_id",
                                       "locus_tag","gene_length","gene_type")
#' 3,847,253 genes (including redundancies):
nrow(prokka.contig.gene.map.df)
head(prokka.contig.gene.map.df) %>%
  knitr::kable(format = "simple")%>%
  print()
#' 1,633,598 contigs had some genes from PROKKA (no matter how small or 
#' whether they're not CDS):
prokka.contig.gene.map.df%>%
  distinct(sample_id,contig_id)%>%
  nrow
#' 3771950 CDS in total:
prokka.contig.gene.map.df<-prokka.contig.gene.map.df%>%
  filter(gene_type=="CDS")
prokka.contig.gene.map.df%>%
  distinct(locus_tag)%>%
  nrow
#' 3764752 CDS >=100 bp:
prokka.contig.gene.map.df<-prokka.contig.gene.map.df%>%
  filter(gene_type=="CDS",gene_length>=100)
prokka.contig.gene.map.df%>%
  distinct(locus_tag)%>%
  nrow
#' 1,618,027 contigs had CDS that were >=100 bp long
prokka.contig.gene.map.df%>%
  distinct(sample_id,contig_id)%>%
  nrow

#+ echo=FALSE
## 4. Add gene count data from HTseq. Then, normalise to TPM. ####
#'
#' ## Add gene count data from HTseq. Then, normalise to TPM.
#' Sources: https://metagenomics-workshop.readthedocs.io/en/latest/annotation/quantification.html  
#' https://github.com/EnvGen/metagenomics-workshop/blob/master/in-house/tpm_table.py
htseq.gene_counts.df<-data.table::fread(htseq.gene_counts.fname,
                                 header = FALSE,sep="\t")%>%
  as_tibble()
colnames(htseq.gene_counts.df)<-c("locus_tag","raw_count")
#' 3,847,253 rows, like in prokka.contig.gene.map.df
nrow(htseq.gene_counts.df)
head(htseq.gene_counts.df) %>%
  knitr::kable(format = "simple")%>%
  print()

#+ echo=FALSE
## 5. MMseqs2 data has mapping between all locus_tag and representative clusters. ####
#'
#' ## MMseqs2 data has mapping between all locus_tag and representative clusters.
mmseqs.clusters.df <- data.table::fread(mmseqs.clusters.fname,
                               header = FALSE,sep = "\t")%>%
  as_tibble()
colnames(mmseqs.clusters.df)<-c("locus_tag_cluster","locus_tag")
head(mmseqs.clusters.df) %>%
  knitr::kable(format = "simple")%>%
  print()
#' 3,764,752 rows
nrow(mmseqs.clusters.df)
#' 1,922,857 unique clusters
mmseqs.clusters.df%>%
  distinct(locus_tag_cluster)

#' Map locus_tag_clusters to samples: 3847253 entries (all mappings, 
#' including redundancies).  
#' Add gene counts: raw counts, so we need to normalise.
gene.annotation.df <- mmseqs.clusters.df%>%
  full_join(prokka.contig.gene.map.df,by = join_by(locus_tag))
#' Most genes are short:
gene.annotation.df%>%
  filter(gene_type=="CDS",!is.na(locus_tag_cluster))%>%
  distinct(locus_tag_cluster,.keep_all = T)%>%
  summary()

#+ echo=FALSE
## 6. TPM normalisation. ####
#'
#' ## TPM normalisation.
#' Combine locus_tag with their counts. Since we have locus_tag_cluster info,
#' we know the abundance for clusters is the sum of all locus_tag 
#' counts in the cluster. We're working with CDS >=100 bp long.
gene.annotation.df <- gene.annotation.df%>%
  full_join(htseq.gene_counts.df,by = join_by(locus_tag))%>%
  filter(!is.na(locus_tag_cluster))
head(gene.annotation.df) %>%
  knitr::kable(format = "simple")%>%
  print()
#' Find representative gene lengths (clustered or unclustered gene length):
representative.gene_lengths <- gene.annotation.df%>%
  group_by(locus_tag_cluster)%>%
  select(locus_tag_cluster,locus_tag,gene_length)%>%
  summarise(rep_gene_len=gene_length[locus_tag == locus_tag_cluster][1])
head(representative.gene_lengths) %>%
  knitr::kable(format = "simple")%>%
  print()
#' Sum counts within each sample:
tpm.df <- gene.annotation.df%>%
  group_by(sample_id,locus_tag_cluster)%>%
  summarise(locus_tag_cluster_sum_count=sum(raw_count))%>%
  ungroup%>%
  # Add gene_count*read_length/representative_gene_length
  left_join(representative.gene_lengths)%>%
  mutate(rgxrl_over_gene_len=(locus_tag_cluster_sum_count*150)/rep_gene_len)%>%
  # Get the final TPM
  group_by(sample_id)%>%
  mutate(T_denom=sum(rgxrl_over_gene_len),
         tpm=(locus_tag_cluster_sum_count*150*10^6) / (rep_gene_len * T_denom))%>%
  ungroup
head(tpm.df) %>%
  knitr::kable(format = "simple")%>%
  print()

rm(prokka.contig.gene.map.df)
rm(mmseqs.clusters.df)
rm(htseq.gene_counts.df)
gc()

#+ echo=FALSE
## 7. Add KofamScan annotation. ####
#'
#' ## Add KofamScan annotation.
kofamscan.df <- read.table(kofamscan.fname,header = TRUE,sep = "\t")%>%
  as_tibble()
colnames(kofamscan.df)[1] <- "locus_tag_cluster"
head(kofamscan.df) %>%
  knitr::kable(format = "simple")%>%
  print()
#' Not all mmseqs clusters are in kofamscan.df, but all kofamscan.df are 
#' in clusters (obviously)
table(unique(gene.annotation.df$locus_tag_cluster)%in%kofamscan.df$locus_tag_cluster)
table(kofamscan.df$locus_tag_cluster%in%unique(gene.annotation.df$locus_tag_cluster))

#+ echo=FALSE
## 8. Add dbcan annotation. ####
#'
#' ## Add dbcan annotation.
dbcan.df <- read.table(dbcan.fname,header=T,sep = "\t",
                     stringsAsFactors = FALSE,comment.char = "")%>%
  as_tibble()
colnames(dbcan.df)<-c("locus_tag_cluster","ec","hmmer","dbcan_sub","diamond","n_of_tools")
head(dbcan.df) %>%
  knitr::kable(format = "simple")%>%
  print()

#' Not all mmseqs clusters are in dbcan, but all dbcan are in clusters (obviously)
table(unique(gene.annotation.df$locus_tag_cluster)%in%dbcan.df$locus_tag_cluster)
table(dbcan.df$locus_tag_cluster%in%unique(gene.annotation.df$locus_tag_cluster))
#' 16752 CAZymes were identified by 1 tool
dbcan.df%>%
  count(n_of_tools)
#' Filter for confident CAZyme hits (≥2 tools agree): count how many tools
#' gave a hit per gene:
dbcan.df <- dbcan.df %>%
  # By default, use hmmer. But if hmmer is unclassified, use dbcan_sub
  mutate(caz_with_aa_coord = ifelse(hmmer != "-", hmmer, dbcan_sub),
         caz_with_aa_coord = ifelse(hmmer=="-" & dbcan_sub == "-", diamond, caz_with_aa_coord))%>%
  # Remove empty entries ('-') for tools that didn't identify a locus_tag_cluster
  select(-ec,-n_of_tools,-hmmer,-dbcan_sub,-diamond)%>%
  distinct(locus_tag_cluster, caz_with_aa_coord, .keep_all = T)  %>% # Remove duplicates
  # Remove text in brackets and CBM (e.g. "GH5(1-100)+CBM50(110-180)")
  mutate(caz_subclass_combination = gsub("_[0-9]*|_e[0-9]*","", caz_with_aa_coord),
         caz_subclass_combination = gsub("\\([0-9]+\\-[0-9]+\\)","", caz_subclass_combination),
         caz_subclass_combination = gsub("in[A-Z]+[0-9]+","", caz_subclass_combination)
         )%>%
  # Split multi-domain strings (e.g. "GH5+CBM50")
  mutate(caz_subclass = caz_subclass_combination)%>%
  separate_longer_delim(caz_subclass, delim = "+") %>%
  mutate(caz_class=ifelse(grepl("[A-Za-z]+[0-9]*$",caz_subclass),
                          gsub("[0-9]*$","",caz_subclass),caz_subclass))
#' Extract caz_subclass name: at the start of the string, match any character
#' except ( for one or more times (that's why we add + sign).
#' Matches a sequence of one or more characters until the first () from the
#' beginning of the string.  
#' ```mutate(caz_subclass = str_extract(caz, "^[^\\(]+"),  
#'        caz_subclass=gsub("_[0-9]*|_e[0-9]*","",caz_subclass))```
#' CAZ class is caz subclass without numbers at the end. For example, GH1 becomes GH.
head(dbcan.df)%>%
  knitr::kable(format = "simple")%>%
  print()

dbcan.df%>%
  nrow()

#' 24063 CAZymes
dbcan.df%>%
  summarise(n_caz=n_distinct(caz_with_aa_coord),
            n_caz_subclass_combination=n_distinct(caz_subclass_combination),
            n_caz_subclass=n_distinct(caz_subclass),
            n_caz_class=n_distinct(caz_class))

#+ echo=FALSE
## 9. Map dbCAN and kofamscan locus_tag_clusters with MMseqs locus_tag_clusters and samples. ####
#'
#' ## Map dbCAN and kofamscan locus_tag_clusters with MMseqs locus_tag_clusters and samples.
# Then, add columns that show whether a locus_tag_cluster was classified by dbcan
# or by kofamscan, or both
# gene.annotation.df <- gene.annotation.df %>%
#   full_join(dbcan.df,by = join_by(locus_tag_cluster),
#             relationship = "many-to-many")%>%
#   full_join(kofamscan.df[,c("locus_tag_cluster","ko")],
#             by = join_by(locus_tag_cluster))%>%
#   mutate(annotation_type=ifelse(is.na(caz)&is.na(ko),"None","Annotated"),
#          annotation_type=ifelse(!is.na(caz)&is.na(ko),"dbCAN",annotation_type),
#          annotation_type=ifelse(is.na(caz)&!is.na(ko),"KofamScan",annotation_type),
#          annotation_type=ifelse(!is.na(caz)&!is.na(ko),"dbCAN+KofamScan",annotation_type))
# head(gene.annotation.df)
# rm(dbcan.df)
# rm(kofamscan.df)
gc()

#' Save data:
if(!file.exists(gene.annotation.df.fname.rds) & !file.exists(gene.annotation.df.fname.tsv)){
  saveRDS(gene.annotation.df, file = gene.annotation.df.fname.rds)
  data.table::fwrite(gene.annotation.df, file = gene.annotation.df.fname.tsv,
                     quote = TRUE,sep = '\t',row.names = FALSE, col.names = TRUE)
}

if(!file.exists(dbcan.df.fname.rds) & !file.exists(dbcan.df.fname.tsv)){
  saveRDS(dbcan.df, file = dbcan.df.fname.rds)
  data.table::fwrite(dbcan.df, file = dbcan.df.fname.tsv,
                     quote = TRUE,sep = '\t',row.names =FALSE,col.names = TRUE)
}

if(!file.exists(kofamscan.df.fname.rds) & !file.exists(kofamscan.df.fname.tsv)){
  saveRDS(kofamscan.df, file = kofamscan.df.fname.rds)
  data.table::fwrite(kofamscan.df, file = kofamscan.df.fname.tsv,
                     quote = TRUE,sep = '\t',row.names = FALSE,col.names = TRUE)
}

if(!file.exists(representative.gene_lengths.fname.rds) & !file.exists(representative.gene_lengths.fname.tsv)){
  saveRDS(representative.gene_lengths, file = representative.gene_lengths.fname.rds)
  data.table::fwrite(representative.gene_lengths, file = representative.gene_lengths.fname.tsv,
                     quote = TRUE,sep = '\t',row.names = FALSE,col.names = TRUE)
}
if(!file.exists(tpm.df.fname.rds) & !file.exists(tpm.df.fname.tsv)){
  saveRDS(tpm.df, file = tpm.df.fname.rds)
  data.table::fwrite(tpm.df, file = tpm.df.fname.tsv,
                     quote = TRUE,sep = '\t',row.names = FALSE,col.names = TRUE)
}
sessionInfo()
rm(list = ls(all=TRUE))
gc()
