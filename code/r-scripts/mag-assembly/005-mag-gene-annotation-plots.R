#' ---
#' output: 
#'   bookdown::html_document2:
#'      toc: true
#' ---

#' ```{r, setup 005-mag-gene-annotation-plots.R, include=FALSE}
#' knitr::opts_knit$set(root.dir = '/home/rakhimov/projects/metagenome')
#' ```
#+ echo=FALSE
# Gene annotation plots. ####
#' # Gene annotation plots.
#'
#+ echo=FALSE
## Introduction ####
#'
#' ## Introduction
#' Plot gene annotation data.
#' 
#+ echo=FALSE
## 1. Load necessary libraries. ####
#'
#' ## Load necessary libraries.
# install.packages(c("tidyverse"))
library(tidyverse)
library(gplots)
#+ echo=FALSE
## 2. Import datasets. ####
#'
#' ## Import datasets.
source(here::here("config/R/config.R"))# config file with global variables
source(here::here("config/R/themes.R"))# config file with themes

#' Mapping between PROKKA, MMseqs2, samples, contigs, and
#' HTseq quantification datasets:
gene.annotation.df <- readRDS(file = gene.annotation.df.fname.rds)
#' KOfamScan output (KEGG Ortholog classification):
kofamscan.df <- readRDS(file = kofamscan.df.fname.rds)
#' dbCAN output (CAZyme classification):
dbcan.df <- readRDS(file = dbcan.df.fname.rds)
#' Representative gene lengths:
representative.gene_lengths <- readRDS(file = representative.gene_lengths.fname.rds)
#' BLASTP classification of CAZymes:
cazymes.blast.df <- read.table(file = cazymes.blast.df.fname,
                               sep = "\t",header = F,fill=T)%>%
  as_tibble()

colnames(cazymes.blast.df)<-c("qacc","sacc", "qseqid","sseqid", 
                              "pident", "length", "mismatch", "gapopen", 
                              "qstart", "qend", "sstart", "send",
                              "evalue", "bitscore", "sframe","uniparc",
                              "caz_db","Kingdom","Species")
#' Clean the dataset by removing the "fasta" suffix, removing qseqid and sseqid,
#' renaming qacc, and renaming some taxa:
cazymes.blast.df<-cazymes.blast.df%>%
  mutate(caz_db=gsub(".*fasta\\|","",caz_db))%>%
  select(-qseqid,-sseqid)%>%
  rename("locus_tag_cluster"="qacc")%>%
  mutate(Kingdom=ifelse(Kingdom=="","unmapped",Kingdom),
         Species=ifelse(Kingdom=="unmapped","unmapped",Species))
head(cazymes.blast.df)%>%
  knitr::kable(format = "simple")%>%
  print()
#' 3764752 CDS:
gene.annotation.df%>%
  filter(gene_type=="CDS")%>%
  distinct(locus_tag, .keep_all =T )%>%
  nrow
#' Min, max, mean gene lengths
gene.annotation.df%>%
  filter(gene_type=="CDS")%>%
  distinct(locus_tag_cluster, .keep_all =T )%>%
  select(gene_length)%>%
  summary()%>%
  knitr::kable(format = "simple")%>%
  print()

#+ echo=FALSE
## 3. Number of unique genes by annotation in total (kofamscan, dbcan, all). ####
#'
#' ## Number of unique genes by annotation in total (kofamscan, dbcan, all).
#' Inspired by https://www.nature.com/articles/s42003-021-02827-2
gene.annotation.tool.list <- list("dbCAN" = unique(dbcan.df$locus_tag_cluster),
                                  "KofamScan" = unique(kofamscan.df$locus_tag_cluster),
                                  "all" = unique(gene.annotation.df$locus_tag_cluster))
gene.annotation.venn <- venn(gene.annotation.tool.list, show.plot = FALSE)
plot(gene.annotation.venn)
attr(gene.annotation.venn, "intersections")
gene.annotation.venn
#' Get all sets (both shared and unique to each tool)
uniq.gene.counts <- attr(gene.annotation.venn, "intersections")
#' Convert list into dataframe
uniq.gene.counts <- stack(uniq.gene.counts)%>%
  rename("locus_tag_cluster" = "values",
         "annotation_type" = "ind")%>%
  count(annotation_type, sort=TRUE)%>%
  mutate(annotation_type = gsub(":all", "", annotation_type),
         annotation_type = gsub("^all$", "None", annotation_type),
         annotation_type = gsub(":", "+", annotation_type))

head(uniq.gene.counts)%>%
  knitr::kable(format = "simple")%>%
  print()

#+ echo=FALSE
## 4. The barplot of all genes for each annotation (4 bars). ####
#'
#' ## The barplot of all genes for each annotation (4 bars).
uniq.gene.counts.bp<-uniq.gene.counts%>%
  mutate(annotation_type=factor(annotation_type,
                                levels=sort(uniq.gene.counts$annotation_type)))%>%
  ggplot(aes(x=annotation_type,y=n,fill=annotation_type))+
  geom_bar(stat="identity")+
  geom_text(aes(label=n),vjust = -0.3, nudge_y = -0.5,
            size=8)+ # add counts of genes
  scale_fill_viridis_d(option = "C")+
  labs(y="Number of unique genes",x="",
       fill="Annotation type")+
  theme_bw()+
  coord_cartesian(ylim = c(0,max(uniq.gene.counts$n)+200000),
                  expand = F)+
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=20), # size of y axis ticks
        axis.title = element_text(size = 20), # size of axis names
        plot.title = element_text(size = 25), # size of plot title
        plot.caption = element_text(size=23), # size of plot caption
        legend.text = element_text(size = 20), # size of legend text
        legend.title = element_text(size = 25), # size of legend title
        legend.position = "inside", # legend inside the plot
        legend.position.inside=c(0.25,0.8),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) # legend on the left top
#+ fig.height=7, fig.width=8
print(uniq.gene.counts.bp)
for(image.format in image.formats){
  ggsave(filename = paste("number-of-genes-in-total", image.format, sep = "."),
         path = mag.figures,
         plot=uniq.gene.counts.bp,
         width=8, height=6,units="in",
         dpi=300,device = image.format)
}

rm(gene.annotation.tool.list)
rm(gene.annotation.venn)
rm(uniq.gene.counts)
rm(uniq.gene.counts.bp)
gc()

#+ echo=FALSE
## 5. Most abundant genes in total with and without filtering by prevalence. ####
#'
#' ## Most abundant genes in total with and without filtering by prevalence.
top.abundant.genes <- gene.annotation.df %>%
  group_by(sample_id)%>%
  count(locus_tag_cluster,name = "n_locus_tags")%>% # we're counting locus_tag_cluster, which are locus_tags
  ungroup()%>%
  arrange(desc(n_locus_tags))
#' Group by KOs:
top.ko <- top.abundant.genes%>%
  left_join(kofamscan.df, by = "locus_tag_cluster")%>%
  mutate(ko_definition = ifelse(is.na(ko),"KO unavailable",ko_definition))%>%
  group_by(sample_id, ko)%>%
  summarise(n_locus_tags = sum(n_locus_tags))%>%
  ungroup() %>%
  left_join(kofamscan.df[,c("ko", "ko_definition")]%>%distinct(), by = "ko")%>%
  arrange(desc(n_locus_tags))

#' Group by CAZymes:
top.caz_subclasses <- top.abundant.genes %>%
  left_join(dbcan.df[,c("locus_tag_cluster","caz_with_aa_coord",
                        "caz_subclass")]%>%distinct(),
            relationship = "many-to-many")%>%
  group_by(sample_id,caz_subclass)%>%
  summarise(n_locus_tags = sum(n_locus_tags))%>%
  ungroup()%>%
  left_join(dbcan.df[,c("caz_subclass", "caz_class")]%>%distinct())%>%
  arrange(desc(n_locus_tags))

#' Count prevalence of locus_tag_clusters(the number of samples it was observed in):
top.abundant.genes.prevalence <- top.abundant.genes %>%
  distinct(sample_id, locus_tag_cluster)%>%
  count(locus_tag_cluster, name = "prevalence", sort = TRUE)

#' Same for KOs:
ko.prevalence <- top.abundant.genes %>%
  left_join(kofamscan.df, by = "locus_tag_cluster")%>%
  mutate(ko_definition = ifelse(is.na(ko),"KO unavailable",ko_definition))%>%
  distinct(sample_id, ko)%>%
  count(ko, name = "prevalence", sort = TRUE)

#' Same for CAZymes:
caz_subclass.prevalence <- top.abundant.genes %>%
  left_join(dbcan.df[,c("locus_tag_cluster","caz_with_aa_coord",
                        "caz_subclass")]%>%distinct(),
            relationship = "many-to-many")%>%
  distinct(sample_id, caz_subclass)%>%
  count(caz_subclass, name = "prevalence", sort = TRUE)

head(top.abundant.genes)%>%
  knitr::kable(format = "simple")%>%
  print()
head(top.ko)%>%
  knitr::kable(format = "simple")%>%
  print()
head(caz_subclass.prevalence)%>%
  knitr::kable(format = "simple")%>%
  print()
  

#' Initial plot shows that top genes are short and don't have KO annotation
#' (but prevalence is low too).
top15.ko.initial <- top.ko %>%
  filter(!is.na(ko))%>%
  arrange(desc(n_locus_tags))%>%
  slice_head(n=1000) %>% # heuristic because we expect 15 most abundant genes to be within 1000 rows
  group_by(ko)%>%
  mutate(ko_row_num = row_number()) %>%
  ungroup()%>%
  filter(ko_row_num == 1 )%>%
  select(-ko_row_num, -sample_id)%>%
  slice_head(n=15)
top15.ko.initial%>%
  knitr::kable(format = "simple")%>%
  print()


top.ko.initial.boxplot <- top.ko %>%
  filter(ko %in% top15.ko.initial$ko)%>%
  # mutate(locus_tag_cluster = factor(locus_tag_cluster,
  #                                   levels = unique(locus_tag_cluster)))%>%
  mutate(ko = factor(ko, levels = unique(ko)),
         ko_definition = factor(ko_definition, levels = unique(ko_definition)))%>%
  # ggplot(aes(x=locus_tag_cluster,y=tpm,fill=ko_definition))+
  ggplot(aes(x = ko, y = n_locus_tags, fill = ko_definition))+
  geom_boxplot()+
  geom_jitter()+
  # geom_text(subset(top.ko,
  #                  locus_tag_cluster%in%top15.ko.initial$locus_tag_cluster) %>%
  #             distinct(locus_tag_cluster,.keep_all = TRUE),
  #           inherit.aes = TRUE,
  #           mapping=aes(x=locus_tag_cluster,
  #                       # y=tpm,
  #                       y = n_locus_tags,
  #                       label=rep_gene_len),
  #           vjust = -0.3, nudge_y = -0.5,
  #           size=6)+ # add counts of genes
  theme_bw()+
  labs(x="Gene ID",
       # y="Median TPM",
       y="Number of identified sequences per sample",
       fill="KO annotation",
       title = "Initial plot shows low prevalence of abundant genes")+
  theme(plot.margin=unit(c(1,1,1,1.5), 'cm'),
        axis.text.x = element_text(size=10,angle = 45,hjust=1),
        axis.text.y = element_text(size=10), # size of y axis ticks
        axis.title = element_text(size = 10), # size of axis names
        plot.title = element_text(size = 15), # size of plot title
        plot.caption = element_text(size=13), # size of plot caption
        legend.text = element_text(size = 10), # size of legend text
        legend.title = element_text(size = 15), # size of legend title
        legend.position = "right" ,
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())# legend on the right

#+ fig.height=8, fig.width=11
print(top.ko.initial.boxplot)

for(image.format in image.formats){
  ggsave(filename = paste("top-ko-initial", image.format, sep = "."),
         path = mag.figures,
         plot=top.ko.initial.boxplot,
         width=10, height=5,units="in",
         units = "px",dpi=300,device = image.format)
}

#+ echo=FALSE
### 5.1 Keep genes with prevalence >2 samples out of 11. ####
#'
#' ### Keep genes with prevalence >2 samples out of 11.
top15.ko.filtered <- ko.prevalence %>%
  filter(prevalence > 2)%>%
  right_join(top.ko, by =join_by(ko))%>%
  filter(!is.na(ko), ko_definition != "KO unavailable")%>%
  arrange(desc(n_locus_tags))%>%
  distinct(ko, .keep_all = T)%>%
  slice_head(n=15)

top.ko.filtered.boxplot <- top.ko %>%
  filter(ko %in% top15.ko.filtered$ko)%>%
  mutate(ko_definition = ifelse(is.na(ko),"KO unavailable",ko_definition)
                       )%>%
  mutate(ko = factor(ko, levels = unique(ko)),
                    ko_definition = factor(ko_definition, levels = unique(ko_definition)))%>%
  # ggplot(aes(x=locus_tag_cluster,y=tpm,fill=ko_definition))+
  ggplot(aes(x = ko,y = n_locus_tags,fill=ko_definition))+
  geom_boxplot(outlier.size = NULL)+
  geom_jitter()+
  # geom_text(subset(top.ko,
  #                  locus_tag_cluster%in%top15.genes.filtered$locus_tag_cluster)%>%
  #             distinct(locus_tag_cluster,.keep_all = TRUE),
  #           inherit.aes = TRUE,
  #           mapping=aes(x=locus_tag_cluster,
  #                       # y=tpm,
  #                       y = n_locus_tags,
  #                       label=rep_gene_len),
  #           vjust = -0.3, nudge_y = -0.5,
  #           size=3)+ # add counts of genes
  theme_bw()+
  scale_fill_viridis_d(option = "C",
                       alpha = 0.5)+
  labs(x="Gene ID",
       # y="Median TPM",
       y="Number of identified sequences per sample",
       # title = "Only prevalent genes (>2/11 samples) are shown here",
       fill="KO annotation"
       )+
  theme(
    plot.margin=unit(c(0.2,0.2,0.2,0.5), 'cm'),
    axis.text.x = element_text(size=10,angle = 45,hjust=1),
    axis.text.y = element_text(size=10), # size of y axis ticks
    axis.title = element_text(size = 10), # size of axis names
    plot.title = element_text(size = 15), # size of plot title
    plot.caption = element_text(size=13), # size of plot caption
    legend.text = element_text(size = 10), # size of legend text
    legend.title = element_text(size = 15), # size of legend title
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.position = "right" )

#+ fig.height=5, fig.width=10
print(top.ko.filtered.boxplot)

for(image.format in image.formats){
  ggsave(filename = paste("top-genes-filtered-prevalence", image.format,sep = "."),
         path = mag.figures,
         plot=top.ko.filtered.boxplot,
         width=10, height=5,units="in",
         dpi=300,device = image.format)
}


#+ echo=FALSE
## 6. Most abundant CAZyme subclasses. ####
#'
#' ## Most abundant CAZyme subclasses.
#' Get the number of samples where a cazyme subclass was found
top15.caz_subclasses <- top.caz_subclasses%>%
  group_by(caz_subclass)%>%
  mutate(caz_subclass_row_num = row_number()) %>%
  ungroup()%>%
  filter(caz_subclass_row_num == 1 )%>%
  select(-caz_subclass_row_num)%>%
  filter(!is.na(caz_subclass))%>%
  slice_head(n=15)

# top.caz_subclasses<-gene.annotation.df%>%
#   select(sample,caz_subclass,caz_class,locus_tag_cluster,tpm)%>%
#   filter(!is.na(caz_subclass))%>%
#   group_by(caz_subclass)%>%
#   mutate(n_samples=n_distinct(sample))%>%
#   add_tally(name="prevalence")%>% # prevalence is the number of samples where a subclass was found
#   mutate(median_tpm=median(tpm))%>% # median TPM for each gene
#   ungroup%>%
#   distinct(sample,locus_tag_cluster,.keep_all = T)%>%
#   arrange(desc(median_tpm))
#   
# top.caz_subclasses<-top.caz_subclasses%>%
#   left_join(unique(cazymes.blast.df[,c("locus_tag_cluster",
#                                        "Kingdom","Species")])) # add BLAST annotation

# top15.caz_subclasses<-top.caz_subclasses%>%
#   select(caz_subclass,median_tpm,caz_class)%>%
#   distinct(caz_subclass,median_tpm,.keep_all = TRUE)%>%
#   arrange(desc(median_tpm))%>%
#   slice_head(n=15)

top.cazyme.plot <- top.caz_subclasses%>%
  filter(caz_subclass%in%top15.caz_subclasses$caz_subclass)%>%
  group_by(sample_id,caz_subclass)%>%
  summarise(n_locus_tags = sum(n_locus_tags))%>%
  ungroup()%>%
  left_join(top.caz_subclasses[,c("caz_subclass", "caz_class")]%>%distinct())%>%
  mutate(caz_subclass=factor(caz_subclass,levels=unique(top15.caz_subclasses$caz_subclass)),
                       caz_class=factor(caz_class,levels=unique(top15.caz_subclasses$caz_class)))%>%
  # ggplot(aes(x=caz_subclass,y=tpm,fill=caz_class))+
  ggplot(aes(x=caz_subclass,y=n_locus_tags,fill=caz_class))+
  geom_boxplot()+
  geom_jitter(cex = 0.8)+
  labs(x="CAZyme subclass",
       # y="Median TPM",
       y="Number of identified sequences per sample",
       fill="CAZyme class")+
  theme_bw()+
  scale_fill_viridis_d(option = "C",
                       alpha = 0.5)+
  # coord_cartesian(expand = FALSE)+
  theme(axis.text.x = element_text(size=12,angle=45,vjust = 1,hjust = 1),
                      axis.text.y = element_text(size=10), # size of y axis ticks
                      axis.title = element_text(size = 10), # size of axis names
                      plot.title = element_text(size = 15), # size of plot title
                      plot.caption = element_text(size=13), # size of plot caption
                      legend.text = element_text(size = 10), # size of legend text
                      legend.title = element_text(size = 15), # size of legend title
                      panel.grid.minor = element_blank(),
                      panel.grid.major = element_blank(),
                      legend.position = "right") # legend on the right

#+ fig.height=5, fig.width=10
print(top.cazyme.plot)

for(image.format in image.formats){
  ggsave(filename =paste( "top-cazymes", image.format, sep = "."),
         path = mag.figures,
         plot=top.cazyme.plot,
         width=10, height=5,units="in",
         dpi=300,device = image.format)
}

#+ echo=FALSE
## 7. Summary of the BLASTP results. ####
#'
#' ## Summary of the BLASTP results.
cazymes.blast.df%>%
  count(Kingdom,sort=TRUE)

#' Find species with the highest number of cazymes
top.cazyme_species<-cazymes.blast.df%>%
  add_count(Species,sort=TRUE)%>%
  distinct(Species,Kingdom,n)%>%
  rename(n_cazymes=n)
head(top.cazyme_species)
top.cazyme_species.fname <- file.path(mag.tables, "cazymes-blastp-taxonomy-top-species.tsv")
if(!file.exists(top.cazyme_species.fname)){
  write.table(top.cazyme_species,
              file=top.cazyme_species.fname,
              sep="\t",col.names = T,row.names = F)
  
}

#' Check how many species the top 15 CAZyme subclasses were found in BLASTP 
#' (showing the first ten):
cazymes.blast.df%>%
  filter(caz_db %in%top15.caz_subclasses$caz_subclass)%>%
  rename("caz_subclass" = "caz_db")%>%
  group_by(caz_subclass)%>%
  count(Species, name = "n_cazymes", sort = TRUE)%>%
  ungroup() %>%
  slice_head(n = 10)%>%
  knitr::kable(format = "simple")%>%
  print()

# top.caz_subclasses%>%
#   filter(caz_subclass%in%top15.caz_subclasses$caz_subclass)%>%
#   group_by(caz_subclass)%>%
#   add_count(Species)%>%
#   # count(Species,sort=TRUE)%>%
#   filter(caz_subclass%in%c("GH45","GH11","GH44"))%>%
#   distinct(caz_subclass,Species,Kingdom)

#' Only for GH subclass:
cazymes.blast.df%>%
  rename("caz_subclass" = "caz_db")%>%
  left_join(top.caz_subclasses[,c("caz_subclass", "caz_class")]%>%distinct(),
                          by = join_by("caz_subclass"))%>%
  filter(caz_class=="GH")%>%
  group_by(caz_subclass)%>%
  add_count(Species)%>%
  ungroup%>%
  # filter(caz_subclass%in%c("GH45","GH11","GH44"))%>%
  distinct(caz_subclass,Species,Kingdom)%>%
  count(Species,sort=TRUE)%>%
  head()%>%
  knitr::kable(format = "simple")%>%
  print()

# top.caz_subclasses%>%
#   filter(caz_class=="GH")%>%
#   group_by(caz_subclass)%>%
#   add_count(Species)%>%
#   ungroup%>%
#   # filter(caz_subclass%in%c("GH45","GH11","GH44"))%>%
#   distinct(caz_subclass,Species,Kingdom)%>%
#   count(Species,sort=TRUE)

#' Save data.
rm(top.ko.initial.boxplot)
rm(top.ko.filtered.boxplot)
rm(top.cazyme.plot)
rm(top15.caz_subclasses)
rm(top15.genes.initial)
rm(top15.ko.filtered)
rm(top15.ko.initial)
rm(cazymes.blast.df)
rm(representative.gene_lengths)
gc()
# save.image("./output/rdafiles/gene-annotation-workspace.RData")


#+ echo=FALSE
## 8. Find core CAZymes and KOs and compare with iMGMC. ####
#'
#' ## Find core CAZymes and KOs and compare with iMGMC.
# load("./output/rdafiles/gene-annotation-workspace.RData")
core.ko <- ko.prevalence %>%
  left_join(kofamscan.df)%>%
  rename("n_samples"="prevalence")%>%
  filter(!is.na(ko),
         n_samples>5)%>%
  filter(!is.na(ko))%>%
  distinct(ko,ko_definition)
head(core.ko)%>%
  knitr::kable(format = "simple")%>%
  print()
nrow(core.ko)

core.caz_subclasses <- caz_subclass.prevalence%>%
  left_join(dbcan.df, by = "caz_subclass")%>%
  rename("n_samples"="prevalence")%>%
  filter(n_samples>5)%>%
  select(caz_subclass,caz_class)%>%
  filter(!is.na(caz_subclass))%>%
  distinct(caz_subclass,.keep_all = TRUE)
head(core.caz_subclasses)%>%
  knitr::kable(format = "simple")%>%
  print()
nrow(core.caz_subclasses)
core.ko.fname <- file.path(mag.tables, "core-kos-by-n_samples.tsv")
core.caz_subclasses.fname <- file.path(mag.tables, "core-caz_subclasses-by-n_samples.tsv")
if(!file.exists(core.ko.fname)){
  write.table(core.ko,
              file=core.ko.fname,
              sep = "\t",
              row.names = F)
}
if(!file.exists(core.caz_subclasses.fname)){
  write.table(core.caz_subclasses,
              file=core.caz_subclasses.fname,
              sep = "\t",
              row.names = F)
}

core.ko.names<-core.ko%>%
  distinct(ko)%>%
  pull(ko)
core.caz_subclasses.names<-core.caz_subclasses%>%
  distinct(caz_subclass)%>%
  pull(caz_subclass)
#+ echo=FALSE
### 8.1 Compare core CAZymes. ####
#'
#' ### Compare core CAZymes.
#' Mouse core gut data are found here:
#' https://zenodo.org/records/3631711/files/iMGMC_map_functionality.tar.gz?download=1
#' https://pmc.ncbi.nlm.nih.gov/articles/PMC7059117/
#' Read the iMGMC core CAZyme dataset:
imgmc.caz<-read.table(imgmc.cazy.fname,
                      header = FALSE, sep = "\t")%>%
  as_tibble()
#' Tidy up the data:
imgmc.caz.vector<-imgmc.caz%>%
  rename("caz"="V2")%>%
  distinct(caz)%>%
  mutate(caz_split=sub("[A-Z0-9\\._]*\\|","",caz),
         caz_split=gsub("_[0-9]+","",caz_split),
         caz_split=gsub("\\+","|",caz_split))%>%
  pull(caz_split)
#' Get unique subclasses:
imgmc.caz.str<-paste(imgmc.caz.vector,collapse = "|")
imgmc.caz.all<-str_split(imgmc.caz.vector, "\\|")
imgmc.caz.all<-unlist(imgmc.caz.all,use.names = F)
imgmc.caz.all<-unique(imgmc.caz.all)
#' Find which core NMR subclasses are absent in core mouse data.
nmr.imgmc.core.caz.venn <- venn(list ("NMR" = core.caz_subclasses.names,
                                      "iMGMC" = imgmc.caz.all))
nmr.core.caz.vector <- attr(nmr.imgmc.core.caz.venn, "intersections")$NMR
length(nmr.core.caz.vector)
#' Save the core NMR subclasses.
nmr.core.caz.df <- nmr.core.caz.vector%>%
  as_tibble()%>%
  rename("caz_subclass"="value")%>%
  inner_join(core.caz_subclasses)
head(nmr.core.caz.df)
nmr.core.caz.df.fname <- file.path(mag.tables, "core-caz_subclasses-nmr_specific.tsv")
if(!file.exists(nmr.core.caz.df.fname)){
  write.table(nmr.core.caz.df,
              file = nmr.core.caz.df.fname,
              sep = "\t",
              row.names = F)
  
}

#+ echo=FALSE
### 8.2 Compare core KOs. ####
#'
#' ### Compare core KOs.
imgmc.ko<-read.table(imgmc.ko.fname,
                     header = FALSE, sep = "\t")%>%
  as_tibble()

imgmc.ko.vector<-imgmc.ko%>%
  rename("ko"="V2")%>%
  distinct(ko)%>%
  pull(ko)
#' Find the core NMR-specific KOs:
nmr.imgmc.core.ko.venn <- venn(list ("NMR" = core.ko.names,
                                     "iMGMC" = imgmc.ko.vector))
nmr.core.ko.vector<-attr(nmr.imgmc.core.ko.venn, "intersections")$NMR
length(nmr.core.ko.vector)

#' Save the dataset:
nmr.core.ko.df<-core.ko%>%
  distinct(ko,.keep_all = TRUE)%>%
  filter(ko%in%nmr.core.ko.vector)%>%
  select(ko,ko_definition)
head(nmr.core.ko.df)
nmr.core.ko.df.fname <- file.path(mag.tables, "core-kos-nmr_specific.tsv")
if(!file.exists(nmr.core.ko.df.fname)){
  write.table(nmr.core.ko.df,
              file = nmr.core.ko.df.fname,
              sep = "\t",
              row.names = F)
  
}

#+ echo=FALSE
## 9. Match gene annotation with taxonomic classification from QIIME2 and Kraken2. ####
#'
#' ## Match gene annotation with taxonomic classification from QIIME2 and Kraken2.
seqkit.all.df <- readRDS(file = seqkit.all.df.fname.rds)
gtdbtk.coverm.drep <- readRDS(file = gtdbtk.coverm.drep.fname.rds)
gtdbtk.coverm.metabat2 <- readRDS(file = gtdbtk.coverm.metabat2.fname.rds)
blastn.coverm.megahit <- readRDS(file = blastn.coverm.megahit.fname.rds)
gtdbtk.taxonomy <- readRDS(gtdbtk.taxonomy.clean.fname.rds)

#' Find p-251-o5 in BLASTN classification of contigs: 14 contigs
blastn.coverm.megahit%>%
  filter(grepl("p-251-o5",classification))%>%
  str
#' Prevotellaceae UCG-003: 32 contigs
blastn.coverm.megahit%>%
  filter(grepl("Prevotellaceae UCG-003",classification))%>%
  str

#' Treponema and Fibrobacter in representative MAGs (GTDB-tk): 8 and 1 MAGs
gtdbtk.taxonomy%>%
  filter(grepl("Treponema|Fibrobacter",Genus))%>%
  arrange(Genus)
#' Erysipelotrichaceae and Prevotellaceae: 5 MAGs
gtdbtk.taxonomy%>%
  filter(grepl("Erysipelotrichaceae|Prevotellaceae",Family))%>%
  arrange(Genus)

#+ echo=FALSE
### 9.1 Which individual MAGs have the highest number of unique CAZyme subclasses. #### 
#'
#' ### Which individual MAGs have the highest number of unique CAZyme subclasses.
mag.with.caz.df<-gtdbtk.taxonomy%>%
  left_join(seqkit.with_taxonomy)%>%
  filter(!is.na(secondary_cluster))%>%
  filter(is_drep_bin)%>%
  select(sample,contig_id,secondary_cluster,
         gtdbtk_result,
         Phylum,
         Family,
         Genus,Species,GC)%>%
  left_join(gene.annotation.df)%>%
  filter(!is.na(caz_subclass))%>%
  select(sample,caz,caz_subclass,caz_class,locus_tag_cluster,
         locus_tag,tpm,
         Family,Genus,secondary_cluster)
str(mag.with.caz.df)

# For each MAG:
mag.with.caz.df%>%
  group_by(secondary_cluster)%>%
  summarise(n_caz_in_mag=n_distinct(caz),
            n_caz_subclass_in_mag=n_distinct(caz_subclass),
            n_cas_class_in_mag=n_distinct(caz_class))%>%
  arrange(desc(n_caz_subclass_in_mag))%>%
  left_join(gtdbtk.taxonomy[,c("sample","bin_id","secondary_cluster","Genus","Species")])

#' For each family: focus on Treponema and Fibrobacter.
#' We get subclasses for the whole family, not each MAG.
mag.with.caz.df%>%
  group_by(Family)%>%
  summarise(n_caz_in_taxon=n_distinct(caz),
            n_caz_subclass_in_taxon=n_distinct(caz_subclass),
            n_cas_class_in_taxon=n_distinct(caz_class))%>%
  arrange(desc(n_caz_subclass_in_taxon))%>%
  left_join(gtdbtk.taxonomy[,c("sample","bin_id","secondary_cluster",
                               "Family","Genus","Species")])%>%
  filter(grepl("Treponema|Erysipelotrichaceae",Family))
#' For Prevotella and Fibrobacter:
mag.with.caz.df%>%
  group_by(Family)%>%
  summarise(n_caz_in_taxon=n_distinct(caz),
            n_caz_subclass_in_taxon=n_distinct(caz_subclass),
            n_cas_class_in_taxon=n_distinct(caz_class))%>%
  arrange(desc(n_caz_subclass_in_taxon))%>%
  left_join(gtdbtk.taxonomy[,c("sample","bin_id",
                               "secondary_cluster","Family","Genus","Species")])%>%
  filter(grepl("Prevotella|Fibrobacter",Genus))

#' Break down each MAG by CAZyme subclasses: get the highest number of CAZymes 
#' for each subclass in each MAG.
treponema.caz<-mag.with.caz.df%>%
  group_by(Family,caz_subclass)%>%
  summarise(n_caz_in_taxon=n_distinct(caz))%>%
  arrange(Family,desc(n_caz_in_taxon))%>%
  filter(grepl("Treponema",Family))
head(treponema.caz)

fibrobacter.caz<-mag.with.caz.df%>%
  group_by(Genus,caz_subclass)%>%
  summarise(n_caz_in_taxon=n_distinct(caz))%>%
  arrange(Genus,desc(n_caz_in_taxon))%>%
  filter(grepl("Fibrobacter",Genus))
head(fibrobacter.caz)

erysipelotrichaceae.caz<-mag.with.caz.df%>%
  group_by(Family,caz_subclass)%>%
  summarise(n_caz_in_taxon=n_distinct(caz))%>%
  arrange(Family,desc(n_caz_in_taxon))%>%
  filter(grepl("Erysipelotrichaceae",Family))
head(erysipelotrichaceae.caz)

prevotella.caz<-mag.with.caz.df%>%
  group_by(Genus,caz_subclass)%>%
  summarise(n_caz_in_taxon=n_distinct(caz))%>%
  arrange(Genus,desc(n_caz_in_taxon))%>%
  filter(grepl("Prevotella",Genus))
head(prevotella.caz)

#' How many CAZyme subclasses are found in Treponema but not Fibrobacter: 48
length(setdiff(treponema.caz$caz_subclass,fibrobacter.caz$caz_subclass))
#' Reverse: 23
length(setdiff(fibrobacter.caz$caz_subclass,treponema.caz$caz_subclass))

#' Among all CAZ subclasses in MAGs, four are not found in Treponema:
caz.not.in.treponema<-mag.with.caz.df%>%
  filter(!grepl("Treponema",Family))%>%
  distinct(caz_subclass)%>%
  pull
setdiff(treponema.caz$caz_subclass,caz.not.in.treponema)

#' Four are not found in Fibrobacter:
caz.not.in.fibrobacter<-mag.with.caz.df%>%
  filter(!grepl("Fibrobacter",Genus))%>%
  distinct(caz_subclass)%>%
  pull
setdiff(fibrobacter.caz$caz_subclass,caz.not.in.fibrobacter)

#' Erysipelotrichaceae has all:
caz.not.in.erysipelotrichaceae<-mag.with.caz.df%>%
  filter(!grepl("Erysipelotrichaceae",Family))%>%
  distinct(caz_subclass)%>%
  pull
setdiff(erysipelotrichaceae.caz$caz_subclass,caz.not.in.erysipelotrichaceae)

#' Two are not found in Prevotella
caz.not.in.prevotella<-mag.with.caz.df%>%
  filter(!grepl("Prevotella",Genus))%>%
  distinct(caz_subclass)%>%
  pull
setdiff(prevotella.caz$caz_subclass,caz.not.in.prevotella)


#+ echo=FALSE
## 10. Summary statistics of all high-quality MAGs including gene annotation. ####
#'
#' ## Summary statistics of all high-quality MAGs including gene annotation.
#' Calculate genome size, contig number, GC%, CDS number, 
#' CheckM2 completeness and contamination
metabat2.date_time<-"20250612_13_37_47"
checkm2.high_quality_mags<-read.table("./output/mag_assembly/checkm2_output/20250619_05_47_09_high_quality_mags.tsv",
                                      sep="\t",header = T)%>%
  as_tibble()%>%
  rename("genome"="Name")%>%
  mutate(genome=gsub(paste0(metabat2.date_time,"_"),"",genome))%>%
  separate_wider_delim(genome,delim = "_bin.",names = c("sample","bin_id"))%>%
  mutate(bin_id=as.integer(bin_id))

mag.stats.with_annotation.df<-seqkit.with_taxonomy%>%
  filter(is_drep_bin)%>%
  select(sample,contig_id,bin_id,contig_length,secondary_cluster,GC)%>%
  filter(!is.na(secondary_cluster))%>%
  group_by(secondary_cluster)%>%
  mutate(genome_size_mbp=sum(contig_length)/10^6,
         mag_gc=mean(GC),
         n_contig=n_distinct(contig_id))%>%
  ungroup%>%
  full_join(gene.annotation.df[,c("locus_tag_cluster","locus_tag","sample","contig_id",
                                  "gene_length","gene_type")])%>%
  filter(!is.na(secondary_cluster))%>%
  group_by(secondary_cluster)%>%
  mutate(n_cds = n_distinct(locus_tag,na.rm = T)) %>% # add gene annotations
  ungroup%>%
  distinct(secondary_cluster,.keep_all = T)%>%
  select(-contig_length,-contig_id,-GC)%>%
  left_join(checkm2.high_quality_mags)%>%
  select(-locus_tag_cluster, -locus_tag,-gene_length,-gene_type,
         -Completeness_Model_Used, -Translation_Table_Used,
         -Genome_Size,-GC_Content,-Additional_Notes)%>%
  left_join(gtdbtk.taxonomy)%>%
  arrange(secondary_cluster)
head(mag.stats.with_annotation.df)
nrow(mag.stats.with_annotation.df)
# 
# write.table(mag.stats.with_annotation.df,
#             file="./output/rtables/mag-stats-with-n_cds.tsv",
#             sep = "\t", row.names = F)

sessionInfo()
rm(list = ls(all=TRUE))
gc()