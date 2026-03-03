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
# load("./output/rdafiles/gene-annotation-workspace.RData")
#' Gene detection data:
prokka.date_time<-"20250626_22_11_43"
#' MMseqs clusters (mapping between all proteins and representative sequences):
mmseqs.date_time<-"20250626_22_11_43"
#' Htseq gene counts:
htseq.date_time<-"20250710_19_36_58"
# KofamScan gene annotation:
kofamscan.date_time<-"20250718_09_21_41"
# dbCAN results (on MMseqs proteins):
dbcan.date_time<-"20250626_22_11_43"
#' Directories with datasets:
prokka.gtf.dir<-"output/mag_assembly/prokka_gtf_files"
mmseqs.output.dir<-file.path("output/mag_assembly/mmseqs_output/easy_cluster")
# HTseq gene counts from PROKKA (not normalised):
htseq.gene_count.dir<- "output/mag_assembly/htseq_gene_count"
kofamscan.output.dir<-file.path("output/mag_assembly/kofam_scan_output")
dbcan.output.dir<-file.path("output/mag_assembly/dbcan_output", 
                            paste(dbcan.date_time,"nr_dbcan_proteins",
                                  sep = "_"))
#' Specify input file names:
prokka.contig.gene.map.fname<-file.path(prokka.gtf.dir,
                                        paste(prokka.date_time,
                                              "all_contig_gene_maps.tsv",
                                              sep="_"))
mmseqs.clusters.fname<-file.path(mmseqs.output.dir,
                                 paste(mmseqs.date_time,
                                       "prokka_nr_prot_cluster.tsv",sep = "_"))
kofamscan.fname<-file.path(kofamscan.output.dir,
                           paste(kofamscan.date_time,
                                 "nr_prot_kofam_scan_top_hits.txt",sep="_"))
dbcan.fname<-file.path(dbcan.output.dir,
                       paste(dbcan.date_time,
                             "nr_dbcan_overview.txt",sep = "_"))
htseq.gene_counts.fname<-file.path(htseq.gene_count.dir,
                                   paste(htseq.date_time,"all_prokka_counts.tsv",
                                         sep="_"))
#+ echo=FALSE
## 3. Analyse PROKKA mapping between samples, contigs, and genes. ####
#'
#' ## Analyse PROKKA mapping between samples, contigs, and genes.
prokka.contig.gene.map.df<-data.table::fread(prokka.contig.gene.map.fname,
                                      header = F,sep="\t")%>%
  as_tibble()
colnames(prokka.contig.gene.map.df)<-c("sample","contig_id",
                                       "locus_tag","gene_length","gene_type")
#' 3,847,253 genes (including redundancies):
nrow(prokka.contig.gene.map.df)
head(prokka.contig.gene.map.df)
#' 1,633,598 contigs had some genes from PROKKA (no matter how small or 
#' whether they're not CDS):
prokka.contig.gene.map.df%>%
  distinct(sample,contig_id)%>%
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
  distinct(sample,contig_id)%>%
  nrow

#+ echo=FALSE
## 4. Add gene count data from HTseq. Then, normalise to TPM. ####
#'
#' ## Add gene count data from HTseq. Then, normalise to TPM.
#' Sources: https://metagenomics-workshop.readthedocs.io/en/latest/annotation/quantification.html
#' https://github.com/EnvGen/metagenomics-workshop/blob/master/in-house/tpm_table.py
htseq.gene_counts.df<-data.table::fread(htseq.gene_counts.fname,
                                 header = F,sep="\t")%>%
  as_tibble()
colnames(htseq.gene_counts.df)<-c("locus_tag","raw_count")
#' 3,847,253 rows, like in prokka.contig.gene.map.df
nrow(htseq.gene_counts.df)
head(htseq.gene_counts.df)

#+ echo=FALSE
## 5. MMseqs2 data has mapping between all locus_tag and representative clusters. ####
#'
#' ## MMseqs2 data has mapping between all locus_tag and representative clusters.
mmseqs.clusters.df<-data.table::fread(mmseqs.clusters.fname,
                               header = F,sep = "\t")%>%
  as_tibble()
colnames(mmseqs.clusters.df)<-c("locus_tag_cluster","locus_tag")
head(mmseqs.clusters.df)
#' 3,764,752 rows
nrow(mmseqs.clusters.df)
#' 1,922,857 unique clusters
mmseqs.clusters.df%>%
  distinct(locus_tag_cluster)

#' Map locus_tag_clusters to samples: 3847253 entries (all mappings, 
#' including redundancies).  
#' Add gene counts: raw counts, so we need to normalise.
gene.annotation.df<-mmseqs.clusters.df%>%
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
gene.annotation.df<-gene.annotation.df%>%
  full_join(htseq.gene_counts.df,by = join_by(locus_tag))%>%
  filter(!is.na(locus_tag_cluster))
#' Find representative gene lengths (clustered or unclustered gene length):
representative.gene_lengths<-gene.annotation.df%>%
  group_by(locus_tag_cluster)%>%
  select(locus_tag_cluster,locus_tag,gene_length)%>%
  summarise(rep_gene_len=gene_length[locus_tag == locus_tag_cluster][1])
#' Sum counts within each sample:
tpm.df<-gene.annotation.df%>%
  group_by(sample,locus_tag_cluster)%>%
  summarise(locus_tag_cluster_count=sum(raw_count))%>%
  ungroup%>%
  # Add gene_count*read_length/representative_gene_length
  left_join(representative.gene_lengths)%>%
  mutate(rgxrl_over_gene_len=(locus_tag_cluster_count*150)/rep_gene_len)%>%
  # Get the final TPM
  group_by(sample)%>%
  mutate(T_denom=sum(rgxrl_over_gene_len),
         tpm=(locus_tag_cluster_count*150*10^6) / (rep_gene_len * T_denom))%>%
  ungroup
#' Combine gene annotation data with TPM:
gene.annotation.df<-gene.annotation.df%>%
  left_join(tpm.df[,c("sample","locus_tag_cluster","tpm")],
            by = join_by(sample, locus_tag_cluster))%>%
  select(-raw_count)

rm(tpm.df)
rm(prokka.contig.gene.map.df)
rm(mmseqs.clusters.df)
rm(htseq.gene_counts.df)
gc()

#+ echo=FALSE
## 7. Add KofamScan annotation. ####
#'
#' ## Add KofamScan annotation.
kofamscan.df<-read.table(kofamscan.fname,header = T,sep = "\t")%>%
  as_tibble()
colnames(kofamscan.df)[1]<-"locus_tag_cluster"
head(kofamscan.df)
#' Not all mmseqs clusters are in kofamscan.df, but all kofamscan.df are 
#' in clusters (obviously)
table(unique(gene.annotation.df$locus_tag_cluster)%in%kofamscan.df$locus_tag_cluster)
table(kofamscan.df$locus_tag_cluster%in%unique(gene.annotation.df$locus_tag_cluster))

#' KO definitions for reference
ko.defs<-kofamscan.df%>%
  select(ko,ko_definition)%>%
  distinct()

#+ echo=FALSE
## 8. Add dbcan annotation. ####
#'
#' ## Add dbcan annotation.
dbcan.df<-read.table(dbcan.fname,header=T,sep = "\t",
                     stringsAsFactors = FALSE,comment.char = "")%>%
  as_tibble()
colnames(dbcan.df)<-c("locus_tag_cluster","ec","hmmer","dbcan_sub","diamond","n_of_tools")
head(dbcan.df)

#' Not all mmseqs clusters are in dbcan, but all dbcan are in clusters (obviously)
table(unique(gene.annotation.df$locus_tag_cluster)%in%dbcan.df$locus_tag_cluster)
table(dbcan.df$locus_tag_cluster%in%unique(gene.annotation.df$locus_tag_cluster))
#' 16752 CAZymes were identified by 1 tool
dbcan.df%>%
  count(n_of_tools)
#' Filter for confident CAZyme hits (≥2 tools agree): count how many tools
#' gave a hit per gene:
dbcan.df<-dbcan.df%>%
  # By default, use hmmer. But if hmmer is unclassified, use dbcan_sub
  mutate(caz=ifelse(hmmer!="-",hmmer,dbcan_sub),
         caz=ifelse(hmmer=="-"&dbcan_sub=="-",diamond,caz))%>%
  # filter(n_of_tools>=2)%>%
  # Pivot to long format
  # pivot_longer(cols = -c("locus_tag_cluster","ec","n_of_tools"),
  #                                   names_to = "tool",
  #                                   values_to = "caz") %>%
  # Remove empty entries ('-') for tools that didn't identify a locus_tag_cluster
  select(-ec,-n_of_tools,-hmmer,-dbcan_sub,-diamond)%>%
  distinct(locus_tag_cluster, caz,.keep_all = T)  %>% # Avoid duplicates
  # Remove text in brackets and CBM (e.g. "GH5(1-100)+CBM50(110-180)")
  mutate(caz_combo=gsub("_[0-9]*|_e[0-9]*","",caz),
         caz_combo=gsub("\\([0-9]+\\-[0-9]+\\)","",caz_combo),
         caz_combo=gsub("in[A-Z]+[0-9]+","",caz_combo)#,
         # caz_combo=gsub("\\+CBM[0-9]+|^CBM[0-9]+\\+","",caz_combo)
         )%>%
  # Split multi-domain strings (e.g. "GH5+CBM50")
  mutate(caz_subclass=caz_combo)%>%
  separate_longer_delim(caz_subclass, delim = "+") %>%
  mutate(caz_class=ifelse(grepl("[A-Za-z]+[0-9]*$",caz_subclass),
                          gsub("[0-9]*$","",caz_subclass),caz_subclass))
#' Extract caz_subclass name: at the start of the string, match any character
#' except ( for one or more times (that's why we add + sign).
#' Matches a sequence of one or more characters until the first () from the
#' beginning of the string
#' mutate(caz_subclass = str_extract(caz, "^[^\\(]+"),
#'        caz_subclass=gsub("_[0-9]*|_e[0-9]*","",caz_subclass))
#' CAZ class is caz subclass without numbers at the end. GH1 becomes GH

dbcan.df%>%
  nrow()

#' 24063 cazymes
dbcan.df%>%
  summarise(n_caz=n_distinct(caz),
            n_caz_combo=n_distinct(caz_combo),
            n_caz_subclass=n_distinct(caz_subclass),
            n_caz_class=n_distinct(caz_class))

#+ echo=FALSE
## 9. Map dbCAN and kofamscan locus_tag_clusters with MMseqs locus_tag_clusters and samples. ####
#'
#' ## Map dbCAN and kofamscan locus_tag_clusters with MMseqs locus_tag_clusters and samples.
# Then, add columns that show whether a locus_tag_cluster was classified by dbcan
# or by kofamscan, or both
gene.annotation.df<-gene.annotation.df%>%
  full_join(dbcan.df,by = join_by(locus_tag_cluster),
            relationship = "many-to-many")%>%
  full_join(kofamscan.df[,c("locus_tag_cluster","ko")],
            by = join_by(locus_tag_cluster))%>%
  mutate(annotation_type=ifelse(is.na(caz)&is.na(ko),"None","Annotated"),
         annotation_type=ifelse(!is.na(caz)&is.na(ko),"dbCAN",annotation_type),
         annotation_type=ifelse(is.na(caz)&!is.na(ko),"KofamScan",annotation_type),
         annotation_type=ifelse(!is.na(caz)&!is.na(ko),"dbCAN+KofamScan",annotation_type))
rm(dbcan.df)
rm(kofamscan.df)
gc()

#' Save data:
# saveRDS(gene.annotation.df,file="./output/rdafiles/gene-annotation-df.rds")
# saveRDS(ko.defs,file="./output/rdafiles/ko-defs.rds")
# saveRDS(representative.gene_lengths,file="./output/rdafiles/representative-gene_lengths.rds")
# data.table::fwrite(gene.annotation.df,file="./output/rtables/gene-annotation-df.tsv",
#             quote = T,sep = '\t',row.names = F,col.names = T)
# data.table::fwrite(ko.defs,file="./output/rtables/ko-defs.tsv",
#             quote = T,sep = '\t',row.names = F,col.names = T)
# data.table::fwrite(representative.gene_lengths,file="./output/rtables/representative-gene_lengths.tsv",
#             quote = T,sep = '\t',row.names = F,col.names = T)
sessionInfo()
rm(list = ls(all=TRUE))
gc()