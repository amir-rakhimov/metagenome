library(tidyverse)
image.formats<-c("tiff","png")
barplot.directory<-"./images/barplots/"
boxplot.directory<-"./images/boxplots/"

gene.annotation.df<-readRDS(file="./output/rdafiles/gene-annotation-df.rds")
ko.defs<-readRDS(file="./output/rdafiles/ko-defs.rds")
representative.gene_lengths<-readRDS(file="./output/rdafiles/representative-gene_lengths.rds")
cazymes.blast.df<-read.table(file="./output/mag_assembly/blastp_output/20250802_12_29_03/20250802_12_29_03_blastp_cazy_db_cazymes_top_hits_with_taxa.tsv",
                             sep = "\t",header = F,fill=T)%>%
  as_tibble()
colnames(cazymes.blast.df)<-c("qacc","sacc", "qseqid","sseqid", 
                              "pident", "length", "mismatch", "gapopen", 
                              "qstart", "qend", "sstart", "send",
                              "evalue", "bitscore", "sframe","uniparc",
                              "caz_db","Kingdom","Species")
cazymes.blast.df<-cazymes.blast.df%>%
  mutate(caz_db=gsub(".*fasta\\|","",caz_db))%>%
  select(-qseqid,-sseqid)%>%
  rename("locus_tag_cluster"="qacc")%>%
  mutate(Kingdom=ifelse(Kingdom=="","unmapped",Kingdom),
         Species=ifelse(Kingdom=="unmapped","unmapped",Species))


# 3764752 CDS
gene.annotation.df%>%
  filter(gene_type=="CDS")%>%
  distinct(locus_tag, .keep_all =T )%>%
  nrow

# Min, max, mean gene lengths
gene.annotation.df%>%
  filter(gene_type=="CDS")%>%
  distinct(locus_tag_cluster, .keep_all =T )%>%
  select(gene_length)%>%
  summary()


# Number of unique genes by annotation in total (kofamscan, dbcan, all) (**My Figure 2**) ####
# https://www.nature.com/articles/s42003-021-02827-2
uniq.gene.counts.all<-gene.annotation.df%>%
  distinct(locus_tag_cluster,.keep_all = T)%>%
  count(annotation_type,sort=TRUE)

ko.vec<-gene.annotation.df%>%
  distinct(locus_tag_cluster,.keep_all = T)%>%
  filter(annotation_type%in%c("KofamScan","dbCAN+KofamScan"))%>%
  pull(locus_tag_cluster)
dbcan.vec<-gene.annotation.df%>%
  distinct(locus_tag_cluster,.keep_all = T)%>%
  filter(annotation_type%in%c("dbCAN","dbCAN+KofamScan"))%>%
  pull(locus_tag_cluster)
no_annot.vec<-gene.annotation.df%>%
  distinct(locus_tag_cluster,.keep_all = T)%>%
  filter(annotation_type=="None")%>%
  pull(locus_tag_cluster)

# gene.list<-list(KofamScan=ko.vec,
#                 dbCAN=dbcan.vec,
#                 None=no_annot.vec)
# ggVennDiagram::ggVennDiagram(gene.list,
#                              force_upset = TRUE, order.set.by = "name", 
#                              order.intersect.by = "none")
# 
# 
# vennD<-ggVennDiagram::ggVennDiagram(gene.list,
#                                     # label_alpha = 0,
#                                     # label_geom = "label",
#                                     # label_size = 4,
#                                     )+
#   scale_fill_distiller(palette = "Reds", direction = 1)+
#   scale_x_continuous(expand = expansion(mult = .2))



### The barplot of all genes for each annotation (4 bars) ####
uniq.gene.counts.all.bp<-uniq.gene.counts.all%>%
  mutate(annotation_type=factor(annotation_type,
                                levels=sort(uniq.gene.counts.all$annotation_type)))%>%
  ggplot(aes(x=annotation_type,y=n,fill=annotation_type))+
  geom_bar(stat="identity")+
  geom_text(aes(label=n),vjust = -0.3, nudge_y = -0.5,
            size=8)+ # add counts of genes
  scale_fill_viridis_d(option = "C")+
  labs(y="Number of unique genes",x="",
       fill="Annotation type")+
  theme_bw()+
  coord_cartesian(ylim = c(0,max(uniq.gene.counts.all$n)+200000),
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


## Plot Number of unique genes by annotation in each sample (kofamscan, dbcan, all) (**My Figure 2**) ####
nmr.ages<-readRDS("../amplicon_nmr/output/rdafiles/custom.md.ages.rds")
uniq.gene.counts.per_sample<-gene.annotation.df%>%
  distinct(locus_tag_cluster,.keep_all = T)%>%
  count(annotation_type,sample,sort=TRUE)%>%
  left_join(nmr.ages,by=join_by("sample"=="Sample"))
### The barplot ####
uniq.gene.counts.per_sample.bp<-uniq.gene.counts.per_sample%>%
  mutate(annotation_type=factor(annotation_type,
                                levels=sort(unique(uniq.gene.counts.per_sample$annotation_type))),
         age=paste0("(",age, " years)"),
         sample=paste(sample,age))%>%
  ggplot(aes(x=sample,y=n,fill=annotation_type))+
  geom_bar(stat="identity")+
  labs(y="Number of unique genes",x="Sample",
       fill="Annotation type")+
  theme_bw()+
  coord_cartesian(ylim = c(0,max(uniq.gene.counts.per_sample$n)+50000),
                  expand = F)+
  theme(axis.text.x = element_text(size=20,angle = 45,hjust=1), 
        axis.text.y = element_text(size=20), # size of y axis ticks
        axis.title = element_text(size = 20), # size of axis names
        plot.title = element_text(size = 25), # size of plot title
        plot.caption = element_text(size=23), # size of plot caption
        legend.text = element_text(size = 20), # size of legend text
        legend.title = element_text(size = 25), # size of legend title
        legend.position = "right" # legend on the right
  ) 
## Save the plots ####
for(image.format in image.formats){
  ggsave(paste0(barplot.directory,
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "number-of-genes-in-total",
                      sep = "-"),".",image.format),
         plot=uniq.gene.counts.all.bp,
         width = 2500,height = 1800,
         units = "px",dpi=300,device = image.format)
  ggsave(paste0(barplot.directory,
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "number-of-genes-per-sample",
                      sep = "-"),".",image.format),
         plot=uniq.gene.counts.per_sample.bp,
         width = 4000,height = 2200,
         units = "px",dpi=300,device = image.format)
}

rm(uniq.gene.counts.all)
rm(uniq.gene.counts.all.bp)
rm(uniq.gene.counts.per_sample)
rm(uniq.gene.counts.per_sample.bp)
gc()

# TPM plot: each sample must sum 10^6
# gene.annotation.df%>%
#   select(sample,locus_tag_cluster,tpm,annotation_type)%>%
#   distinct(sample,locus_tag_cluster,.keep_all = T)%>%
#   ggplot(aes(x=sample,y=tpm,fill=annotation_type))+
#   geom_bar(stat="identity")

# Most abundant genes in total with and without filtering by prevalence (**My Figure 3**) ####
# Use median TPM
top.abundant.genes<-gene.annotation.df%>%
  select(sample,locus_tag_cluster,tpm,ko,gene_length,locus_tag)%>%
  distinct(sample,locus_tag_cluster,.keep_all = TRUE)%>%
  group_by(locus_tag_cluster)%>%
  add_tally(name="prevalence")%>% # prevalence is the number of samples where a gene was found
  mutate(median_tpm=median(tpm))%>% # median TPM for each gene
  ungroup%>%
  arrange(desc(median_tpm))

## Initial plot shows that top genes are short and don't have KO annotation (but prevalence is low too) ####
top15.genes.initial<-top.abundant.genes%>%
  distinct(locus_tag_cluster,.keep_all = T)%>%
  slice_head(n=15)

# Get all the information for abundant genes, such as TPM in each sample, KO,
# gene length
top.abundant.genes<-top.abundant.genes%>%
  distinct(sample,locus_tag_cluster,.keep_all = T)%>%
  left_join(ko.defs)%>%
  mutate(locus_tag_cluster=factor(locus_tag_cluster,
                                  levels=unique(top.abundant.genes$locus_tag_cluster)),
         ko_definition=ifelse(is.na(ko),"KO unavailable",ko_definition),
         gene_length=paste(gene_length, "bp"))

# Get the highest TPM for representative genes, add columns with KO definitions
representative.gene_lengths.max_tpm<-top.abundant.genes%>%
  group_by(locus_tag_cluster)%>%
  filter(tpm==max(tpm))%>%
  ungroup%>%
  distinct(locus_tag_cluster,.keep_all = T)%>%
  select(locus_tag_cluster,tpm)%>%
  left_join(representative.gene_lengths)%>%
  left_join(top.abundant.genes[,c("locus_tag_cluster","ko_definition")])%>%
  mutate(rep_gene_len=paste(rep_gene_len, "bp"))%>%
  distinct(locus_tag_cluster,.keep_all = T)

top.abundant.genes.initial.boxplot<-top.abundant.genes%>%
  filter(locus_tag_cluster%in%top15.genes.initial$locus_tag_cluster)%>%
  ggplot(aes(x=locus_tag_cluster,y=tpm,fill=ko_definition))+
  geom_boxplot()+
  geom_jitter()+
  geom_text(subset(representative.gene_lengths.max_tpm,
                   locus_tag_cluster%in%top15.genes.initial$locus_tag_cluster),
            inherit.aes = TRUE,
            mapping=aes(x=locus_tag_cluster,
                        y=tpm,
                        label=rep_gene_len),
            vjust = -0.3, nudge_y = -0.5,
            size=6)+ # add counts of genes
  theme_bw()+
  labs(x="Gene ID",
       y="Median TPM",
       fill="KO annotation",
       title = "Initial plot shows low prevalence of abundant genes")+
  theme(plot.margin=unit(c(1,1,1,1.5), 'cm'),
        axis.text.x = element_text(size=20,angle = 45,hjust=1), 
        axis.text.y = element_text(size=20), # size of y axis ticks
        axis.title = element_text(size = 20), # size of axis names
        plot.title = element_text(size = 25), # size of plot title
        plot.caption = element_text(size=23), # size of plot caption
        legend.text = element_text(size = 20), # size of legend text
        legend.title = element_text(size = 25), # size of legend title
        legend.position = "right" )# legend on the right

for(image.format in image.formats){
  ggsave(paste0(boxplot.directory,
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "top-genes-initial",
                      sep = "-"),".",image.format),
         plot=top.abundant.genes.initial.boxplot,
         width = 5000,height = 3200,
         units = "px",dpi=300,device = image.format)
}


## Keep genes with prevalence >2 samples out of 11 ####
top15.genes.filtered<-top.abundant.genes%>%
  filter(prevalence>2)%>%
  # filter(!is.na(ko),!grepl("transposase|ribosom|replicative|polymerase|topoisomerase|gyrase",ko_definition))%>%
  distinct(locus_tag_cluster,.keep_all = T)%>%
  slice_head(n=15)

top.abundant.genes.filtered.boxplot<-top.abundant.genes%>%
  filter(locus_tag_cluster%in%top15.genes.filtered$locus_tag_cluster)%>%
  mutate(ko_definition=factor(ko_definition,
                              levels=unique(top15.genes.filtered$ko_definition)))%>%
  ggplot(aes(x=locus_tag_cluster,y=tpm,fill=ko_definition))+
  geom_boxplot(outlier.size = NULL)+
  geom_jitter()+
  geom_text(subset(representative.gene_lengths.max_tpm,
                   locus_tag_cluster%in%top15.genes.filtered$locus_tag_cluster),
            inherit.aes = TRUE,
            mapping=aes(x=locus_tag_cluster,
                        y=tpm,
                        label=rep_gene_len),
            vjust = -0.3, nudge_y = -0.5,
            size=3)+ # add counts of genes
  theme_bw()+
  scale_fill_viridis_d(option = "C",
                       alpha = 0.5)+
  labs(x="Gene ID",
       y="Median TPM",
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
        legend.position = "right" )# legend on the right

for(image.format in image.formats){
  ggsave(paste0(boxplot.directory,
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "top-genes-filtered-prevalence",
                      sep = "-"),".",image.format),
         plot=top.abundant.genes.filtered.boxplot,
         width=10, height=5,units="in",
         # width = 7000,height = 3200, units = "px",
         dpi=300,device = image.format)
}

## Filter by KO to remove house-keeping genes ####
top15.genes.filtered.by_ko<-top.abundant.genes%>%
  filter(prevalence>2)%>%
  filter(!is.na(ko),
         !grepl("transposase|ribosom|replicative|polymerase",ko_definition), # used this filter first time
         !grepl("topoisomerase|gyrase",ko_definition), # used this filter first time
         !grepl("chaperonin|endonuclease|resolvase|single-strand",ko_definition))%>% # added this filter recently
  distinct(locus_tag_cluster,.keep_all = T)%>%
  slice_head(n=15)

top.abundant.genes.filtered.by_ko.boxplot<-top.abundant.genes%>%
  filter(locus_tag_cluster%in%top15.genes.filtered.by_ko$locus_tag_cluster)%>%
  mutate(ko_definition=factor(ko_definition,
                              levels=unique(top15.genes.filtered.by_ko$ko_definition)))%>%
  ggplot(aes(x=locus_tag_cluster,y=tpm,fill=ko_definition))+
  geom_boxplot()+
  geom_jitter()+
  geom_text(subset(representative.gene_lengths.max_tpm,
                   locus_tag_cluster%in%top15.genes.filtered.by_ko$locus_tag_cluster),
            inherit.aes = TRUE,
            mapping=aes(x=locus_tag_cluster,
                        y=tpm,
                        label=rep_gene_len),
            vjust = -0.3, nudge_y = -0.5,
            size=6)+ # add counts of genes
  theme_bw()+
  labs(x="Gene ID",
       y="Median TPM",
       fill="KO annotation",
       subtitle = paste0("Genes related to core functions (transposases, ",
                         "ribosomal activity, replication, and DNA repair) are removed"),
       title = "Only prevalent genes (>2/11 samples) are shown here")+
  theme(plot.margin=unit(c(1,1,1,1.5), 'cm'),
        plot.subtitle = element_text(size = 20), # 
        axis.text.x = element_text(size=20,angle = 45,hjust=1), 
        axis.text.y = element_text(size=20), # size of y axis ticks
        axis.title = element_text(size = 20), # size of axis names
        plot.title = element_text(size = 25), # size of plot title
        plot.caption = element_text(size=23), # size of plot caption
        legend.text = element_text(size = 20), # size of legend text
        legend.title = element_text(size = 25), # size of legend title
        legend.position = "right" )# legend on the right

for(image.format in image.formats){
  ggsave(paste0(boxplot.directory,
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "top-genes-filtered-ko",
                      sep = "-"),".",image.format),
         plot=top.abundant.genes.filtered.by_ko.boxplot,
         width = 7000,height = 3200,
         units = "px",dpi=300,device = image.format)
}


# Most abundant CAZyme subclasses (**My Figure 4**) ####
# Get the number of samples where a cazyme subclass was found
# Figure 3a 
# https://www.nature.com/articles/s42003-021-02827-2
top.caz_subclasses<-gene.annotation.df%>%
  select(sample,caz_subclass,caz_class,locus_tag_cluster,tpm)%>%
  filter(!is.na(caz_subclass))%>%
  group_by(caz_subclass)%>%
  mutate(n_samples=n_distinct(sample))%>%
  add_tally(name="prevalence")%>% # prevalence is the number of samples where a subclass was found
  # Should we use distinct(locus_tag_cluster) to remove locus_tag_cluster from multiple samples?
  mutate(median_tpm=median(tpm))%>% # median TPM for each gene
  ungroup%>%
  distinct(sample,locus_tag_cluster,.keep_all = T)%>%
  arrange(desc(median_tpm))
  
top.caz_subclasses<-top.caz_subclasses%>%
  left_join(unique(cazymes.blast.df[,c("locus_tag_cluster",
                                       "Kingdom","Species")])) # add BLAST annotation

top15.caz_subclasses<-top.caz_subclasses%>%
  select(caz_subclass,median_tpm,caz_class)%>%
  distinct(caz_subclass,median_tpm,.keep_all = TRUE)%>%
  arrange(desc(median_tpm))%>%
  slice_head(n=15)

top.cazyme.plot<-top.caz_subclasses%>%
  filter(caz_subclass%in%top15.caz_subclasses$caz_subclass)%>%
  mutate(caz_subclass=factor(caz_subclass,levels=top15.caz_subclasses$caz_subclass),
         caz_class=factor(caz_class,levels=unique(top15.caz_subclasses$caz_class)))%>%
  ggplot(aes(x=caz_subclass,y=tpm,fill=caz_class))+
  geom_boxplot()+
  geom_jitter(cex = 0.8)+
  labs(y="Median TPM",
       x="CAZyme subclass",
       fill="CAZyme class")+
  theme_bw()+
  scale_fill_viridis_d(option = "C",
                       alpha = 0.5)+
  # coord_cartesian(expand = F)+
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


for(image.format in image.formats){
  ggsave(paste0(boxplot.directory,
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "top-cazymes",
                      sep = "-"),".",image.format),
         plot=top.cazyme.plot,
         width=10, height=5,units="in",
         # width = 4000,height = 2200, units = "px",
         dpi=300,device = image.format)
}

# Summary of the BLASTP results ####
cazymes.blast.df%>%
  count(Kingdom,sort=TRUE)

# Find species with the highest number of cazymes
top.cazyme_species<-cazymes.blast.df%>%
  add_count(Species,sort=TRUE)%>%
  distinct(Species,Kingdom,n)%>%
  rename(n_cazymes=n)

write.table(top.cazyme_species,
            file="./output/rtables/cazymes-blastp-taxonomy-top-species.tsv",
            sep="\t",col.names = T,row.names = F)


# Check how many species the top 15 CAZyme subclasses were found in ####
top.caz_subclasses%>%
  filter(caz_subclass%in%top15.caz_subclasses$caz_subclass)%>%
  group_by(caz_subclass)%>%
  add_count(Species)%>%
  # count(Species,sort=TRUE)%>%
  filter(caz_subclass%in%c("GH45","GH11","GH44"))%>%
  distinct(caz_subclass,Species,Kingdom)

top.caz_subclasses%>%
  filter(caz_class=="GH")%>%
  group_by(caz_subclass)%>%
  add_count(Species)%>%
  ungroup%>%
  # filter(caz_subclass%in%c("GH45","GH11","GH44"))%>%
  distinct(caz_subclass,Species,Kingdom)%>%
  count(Species,sort=TRUE)


# 
top.caz_subclasses%>%
  filter(Species%in%top.cazyme_species[1:10,]$Species)%>%
  add_count(caz_subclass,sort=TRUE)%>%
  ggplot(aes(x=Species,y=n,fill=caz_subclass))+
  geom_bar(stat="identity")



# Save data ####
rm(top.abundant.genes.filtered.boxplot)
rm(top.abundant.genes.filtered.by_ko.boxplot)
rm(top.abundant.genes.initial.boxplot)
rm(top.cazyme.plot)
rm(top15.caz_subclasses)
rm(top15.genes.filtered)
rm(top15.genes.filtered.by_ko)
rm(top15.genes.initial)
rm(cazymes.blast.df)
rm(gene.list)
rm(representative.gene_lengths)
rm(representative.gene_lengths.max_tpm)
rm(dbcan.vec)
rm(no_annot.vec)
rm(ko.vec)
rm(vennD)
gc()
save.image("./output/rdafiles/gene-annotation-workspace.RData")

# Find core CAZymes and KOs and compare with iMGMC ####
load("./output/rdafiles/gene-annotation-workspace.RData")
core.ko<-top.abundant.genes%>%
  rename("n_samples"="prevalence")%>%
  filter(!is.na(ko),
         n_samples>5)%>%
  select(ko,ko_definition)%>%
  distinct(ko,.keep_all = T)

core.caz_subclasses<-top.caz_subclasses%>%
  filter(n_samples>5)%>%
  select(caz_subclass,caz_class)%>%
  distinct(caz_subclass,.keep_all = T)

write.table(core.ko,
            file="./output/rtables/core-kos-by-n_samples.tsv",
            sep = "\t",
            row.names = F)
write.table(core.caz_subclasses,
            file="./output/rtables/core-caz_subclasses-by-n_samples.tsv",
            sep = "\t",
            row.names = F)

core.ko.names<-core.ko%>%
  distinct(ko)%>%
  pull(ko)

core.caz_subclasses.names<-core.caz_subclasses%>%
  distinct(caz_subclass)%>%
  pull(caz_subclass)

# Mouse core gut KOs are found here
# https://zenodo.org/records/3631711/files/iMGMC_map_functionality.tar.gz?download=1
# https://pmc.ncbi.nlm.nih.gov/articles/PMC7059117/
imgmc.caz<-read.table("./data/iMGMC_map_functionality/iMGMC-map-GeneID-CAZy.tab",
                                  header = F, sep = "\t")%>%
  as_tibble()


# imgmc.caz<-imgmc.caz%>%
imgmc.caz.vector<-imgmc.caz%>%
  rename("caz"="V2")%>%
  distinct(caz)%>%
  mutate(caz_split=sub("[A-Z0-9\\._]*\\|","",caz),
         caz_split=gsub("_[0-9]+","",caz_split),
         caz_split=gsub("\\+","|",caz_split))%>%
  pull(caz_split)

imgmc.caz.str<-paste(imgmc.caz.vector,collapse = "|")
imgmc.caz.all<-str_split(imgmc.caz.vector, "\\|")
imgmc.caz.all<-unlist(imgmc.caz.all,use.names = F)


imgmc.caz.all<-unique(imgmc.caz.all)
nmr.core.caz.vector<-setdiff(core.caz_subclasses.names,imgmc.caz.all)
length(nmr.core.caz.vector)

nmr.core.caz.df<-nmr.core.caz.vector%>%
  as_tibble()%>%
  rename("caz_subclass"="value")%>%
  inner_join(core.caz_subclasses)

write.table(nmr.core.caz.df,
            file="./output/rtables/core-caz_subclasses-nmr_specific.tsv",
            sep = "\t",
            row.names = F)

### KOs ####
imgmc.ko<-read.table("./data/iMGMC_map_functionality/iMGMC-map-GeneID-KeggKO.tab",
                     header = F, sep = "\t")%>%
  as_tibble()

imgmc.ko.vector<-imgmc.ko%>%
  rename("ko"="V2")%>%
  distinct(ko)%>%
  pull(ko)

nmr.core.ko.vector<-setdiff(core.ko.names,imgmc.ko.vector)
length(nmr.core.ko.vector)

nmr.core.ko.df<-core.ko%>%
  distinct(ko,.keep_all = T)%>%
  filter(ko%in%nmr.core.ko.vector)%>%
  select(ko,ko_definition)
write.table(nmr.core.ko.df,
            file="./output/rtables/core-kos-nmr_specific.tsv",
            sep = "\t",
            row.names = F)


#####
# Find taxa from 16S and Kraken2 data ####
seqkit.with_taxonomy<-readRDS(file = "./output/rdafiles/seqkit-with_taxonomy.rds")
gtdbtk.coverm.drep<-readRDS(file = "./output/rdafiles/gtdbtk-coverm-drep.rds")
gtdbtk.coverm.metabat2<-readRDS(file = "./output/rdafiles/gtdbtk-coverm-metabat2.rds")
blastn.coverm.megahit<-readRDS(file = "./output/rdafiles/blastn-coverm-megahit.rds")
gtdbtk.taxonomy<-readRDS("./output/rdafiles/gtdbtk-taxonomy.rds")

# BLASTN
blastn.coverm.megahit%>%
  filter(grepl("p-251-o5",classification))

blastn.coverm.megahit%>%
  filter(grepl("Prevotellaceae UCG-003",classification))

# GTDBTK
gtdbtk.taxonomy%>%
  filter(grepl("Treponema|Fibrobacter",Genus))%>%
  arrange(Genus)

gtdbtk.taxonomy%>%
  filter(grepl("Erysipelotrichaceae|Prevotellaceae",Family))%>%
  arrange(Genus)

# Which individual MAGs have the highest number of unique CAZyme subclasses ####
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

# For each MAG
mag.with.caz.df%>%
  group_by(secondary_cluster)%>%
  summarise(n_caz_in_mag=n_distinct(caz),
            n_caz_subclass_in_mag=n_distinct(caz_subclass),
            n_cas_class_in_mag=n_distinct(caz_class))%>%
  arrange(desc(n_caz_subclass_in_mag))%>%
  left_join(gtdbtk.taxonomy[,c("sample","bin_id","secondary_cluster","Genus","Species")])

# For each family: Treponema and Fibrobacter
# We get subclasses for the whole family, not each MAG
mag.with.caz.df%>%
  group_by(Family)%>%
  summarise(n_caz_in_taxon=n_distinct(caz),
            n_caz_subclass_in_taxon=n_distinct(caz_subclass),
            n_cas_class_in_taxon=n_distinct(caz_class))%>%
  arrange(desc(n_caz_subclass_in_taxon))%>%
  left_join(gtdbtk.taxonomy[,c("sample","bin_id","secondary_cluster",
                               "Family","Genus","Species")])%>%
  filter(grepl("Treponema|Erysipelotrichaceae",Family))

mag.with.caz.df%>%
  group_by(Family)%>%
  summarise(n_caz_in_taxon=n_distinct(caz),
            n_caz_subclass_in_taxon=n_distinct(caz_subclass),
            n_cas_class_in_taxon=n_distinct(caz_class))%>%
  arrange(desc(n_caz_subclass_in_taxon))%>%
  left_join(gtdbtk.taxonomy[,c("sample","bin_id",
                               "secondary_cluster","Family","Genus","Species")])%>%
  filter(grepl("Prevotella|Fibrobacter",Genus))

 

# Break down each MAG by CAZyme subclasses: get the highest number of CAZymes for each
# subclass in each MAG 
treponema.caz<-mag.with.caz.df%>%
  group_by(Family,caz_subclass)%>%
  summarise(n_caz_in_taxon=n_distinct(caz))%>%
  arrange(Family,desc(n_caz_in_taxon))%>%
  filter(grepl("Treponema",Family))

fibrobacter.caz<-mag.with.caz.df%>%
  group_by(Genus,caz_subclass)%>%
  summarise(n_caz_in_taxon=n_distinct(caz))%>%
  arrange(Genus,desc(n_caz_in_taxon))%>%
  filter(grepl("Fibrobacter",Genus))

erysipelotrichaceae.caz<-mag.with.caz.df%>%
  group_by(Family,caz_subclass)%>%
  summarise(n_caz_in_taxon=n_distinct(caz))%>%
  arrange(Family,desc(n_caz_in_taxon))%>%
  filter(grepl("Erysipelotrichaceae",Family))

prevotella.caz<-mag.with.caz.df%>%
  group_by(Genus,caz_subclass)%>%
  summarise(n_caz_in_taxon=n_distinct(caz))%>%
  arrange(Genus,desc(n_caz_in_taxon))%>%
  filter(grepl("Prevotella",Genus))

length(setdiff(treponema.caz$caz_subclass,fibrobacter.caz$caz_subclass))
length(setdiff(fibrobacter.caz$caz_subclass,treponema.caz$caz_subclass))

caz.not.in.treponema<-mag.with.caz.df%>%
  filter(!grepl("Treponema",Family))%>%
  distinct(caz_subclass)%>%
  pull

setdiff(treponema.caz$caz_subclass,caz.not.in.treponema)


caz.not.in.fibrobacter<-mag.with.caz.df%>%
  filter(!grepl("Fibrobacter",Genus))%>%
  distinct(caz_subclass)%>%
  pull
setdiff(fibrobacter.caz$caz_subclass,caz.not.in.fibrobacter)


caz.not.in.erysipelotrichaceae<-mag.with.caz.df%>%
  filter(!grepl("Erysipelotrichaceae",Family))%>%
  distinct(caz_subclass)%>%
  pull
setdiff(erysipelotrichaceae.caz$caz_subclass,caz.not.in.erysipelotrichaceae)


caz.not.in.prevotella<-mag.with.caz.df%>%
  filter(!grepl("Prevotella",Genus))%>%
  distinct(caz_subclass)%>%
  pull
setdiff(prevotella.caz$caz_subclass,caz.not.in.prevotella)


# Genes in taxa with low GC% ####
# 
# gtdbtk.taxonomy%>%
#   filter(Phylum=="Desulfobacterota"|Phylum=="Methanobacteriota")%>%
#   left_join(seqkit.with_taxonomy)%>%
#   select(sample,contig_id,secondary_cluster,gtdbtk_result,Phylum,Genus,Species,GC)%>%
#   left_join(gene.annotation.df)%>%
#   filter(!is.na(ko))%>%
#   group_by(Phylum,secondary_cluster)%>%
#   distinct(Phylum,ko,tpm,secondary_cluster)%>%
#   left_join(ko.defs)%>%
#   add_count(ko_definition,sort = T)%>%
#   filter(ko_definition!="uncharacterized protein")%>%
#   slice_head(n=20)%>%
#   ggplot(aes(x=ko_definition,y=tpm,fill=Phylum))+
#   geom_boxplot()+coord_flip()



# For BLASTP tsv is not needed anymore. but I'll analyse later ####
# head(top15.caz_subclasses)
# rep.ids<-gene.annotation.df%>%
#   filter(caz_subclass%in%c("GH45","GH11"))%>%
#   distinct(locus_tag_cluster)
# write.table(rep.ids,file="output/rtables/cazymes-for-blast.tsv",
#             row.names = F,col.names = F,quote = F)


# KEGG Pathways (**My Figure 5**) ####
# BiocManager::install("pathview")

library(tidyverse)
library(KEGGREST)
library(pathview)

# TODO: in all data vs MAGs
# Sulfur metabolism [PATH:ko00920]
# Nitrogen metabolism [PATH:ko00910]
# Thiamine metabolism [PATH:ko00730]
# Riboflavin metabolism [PATH:ko00740]
# Vitamin B6 metabolism [PATH:ko00750]
# Nicotinate and nicotinamide metabolism [PATH:ko00760]
# Pantothenate and CoA biosynthesis [PATH:ko00770]
# Biotin metabolism [PATH:ko00780]
# Lipoic acid metabolism [PATH:ko00785]
# Folate biosynthesis [PATH:ko00790]
# One carbon pool by folate [PATH:ko00670]
# Retinol metabolism [PATH:ko00830]
# Porphyrin metabolism [PATH:ko00860]
# Ubiquinone and other terpenoid-quinone biosynthesis [PATH:ko00130]

# 4.4 Cellular community - prokaryotes
# 02024 Quorum sensing
# 05111 Biofilm formation - Vibrio cholerae
# 02025 Biofilm formation - Pseudomonas aeruginosa
# 02026 Biofilm formation - Escherichia coli


# kofamscan.df<-read.table(kofamscan.fname,header = T,sep = "\t")%>%
#   as_tibble()
colnames(kofamscan.df)[1]<-"locus_tag_cluster"
kofamscan.df<-gene.annotation.df%>%
  filter(!is.na(ko))%>%
  distinct(locus_tag_cluster,ko)

kofamscan.df<-seqkit.with_taxonomy%>%
  distinct(sample,contig_id,gtdbtk_result,secondary_cluster)%>%
  filter(grepl("Treponema",gtdbtk_result))%>%
  left_join(gene.annotation.df)%>%
  distinct(ko)%>%
  filter(!is.na(ko))%>%
  left_join(ko.defs)

kegg.pw.ids<-c("ko00920"="Sulfur metabolism",
               "ko00910"="Nitrogen metabolism",
               "ko00730"="Thiamine metabolism",
               "ko00740"="Riboflavin metabolism",
               "ko00750"="Vitamin B6 metabolism",
               "ko00760"="Nicotinate and nicotinamide metabolism",
               "ko00770"="Pantothenate and CoA biosynthesis",
               "ko00780"="Biotin metabolism",
               "ko00785"="Lipoic acid metabolism",
               "ko00790"="Folate biosynthesis",
               "ko00670"="One carbon pool by folate",
               "ko00830"="Retinol metabolism",
               "ko00860"="Porphyrin metabolism",
               "ko00130"="Ubiquinone and other terpenoid-quinone biosynthesis")
kegg.pw.ids<-c("ko00910"="Nitrogen metabolism")
# kegg.pw.ids<-c("ko02024"="Quorum sensing",
#                "ko05111"="Biofilm formation - Vibrio cholerae",
#                "ko02025"="Biofilm formation - Pseudomonas aeruginosa",
#                "ko02026"="Biofilm formation - Escherichia coli",
#                "ko03070"="Bacterial secretion system")
for (i in seq_along(kegg.pw.ids)){
  kegg.pathway.id<-names(kegg.pw.ids[i])
  kegg.pathway.name<-unname(kegg.pw.ids[i])
  
  # kegg.pw<-keggGet(kegg.pathway.id)[[1]]
  # kegg.pathway.orthologs<-names(kegg.pw$ORTHOLOGY)
  # df.kegg<-tibble(ko=ko.ids.in_pw,
  #                   ko_pathway=pathway.id)
  # 
  # length(intersect(ko.sulfur,kofamscan.df$ko))==length(ko.sulfur)
  # inner_join(df.sulfur,kofamscan.df,by="ko")
  
  kegg.pathway.id.num<-gsub("ko","",kegg.pathway.id)
  file.suffix<-tolower(gsub(" ","_",kegg.pathway.name))
  
  pv.out <- pathview(gene.data = kofamscan.df$ko, 
                     pathway.id = kegg.pathway.id.num,
                     species = "ko", 
                     out.suffix =file.suffix, 
                     kegg.native = T,
                     kegg.dir = "./images/ko_pathways")
  file.rename(from=file.path(".",paste(kegg.pathway.id,file.suffix,"png",sep=".")),
              to=file.path(".","images/ko_pathways",paste(kegg.pathway.id,file.suffix,"png",sep=".")))
}



#### Top 40 CAZyme classes in different colonies
top40.caz_subclasses<-gene.annotation.df%>%
  filter(!is.na(locus_tag_cluster))%>%
  group_by(caz_subclass)%>%
  summarise(median_tpm=median(tpm))%>%
  arrange(desc(median_tpm))%>%
  slice_head(n=40)%>%
  mutate(caz_subclass=factor(caz_subclass,levels=caz_subclass))

nmr.relations<-read.table("../amplicon_nmr/data/metadata/pooled-metadata/nmr-relations.tsv",
                          header = T)

caz.colony.tpm<-gene.annotation.df%>%
  left_join(nmr.relations,by=join_by(sample==Sample))%>%
  left_join(nmr.ages,by=join_by(sample==Sample))%>%
  filter(!is.na(caz_subclass))%>%
  group_by(sample,caz_subclass)%>%
  mutate(tpm=sum(tpm))%>%
  distinct(sample,caz_subclass,.keep_all = T)

foo<-caz.colony.tpm%>%
  ungroup%>%
  select(caz_subclass,tpm,sample)%>%
  pivot_wider(names_from = caz_subclass,
              values_from = tpm,
              values_fill = 0)%>%
  as.data.frame()%>%
  column_to_rownames("sample")
heatmap(as.matrix(foo),scale = "row")








# Read counts ####
## Raw read counts ####
raw.read.counts<-read.table("./output/mag_assembly/seqkit_output/20250728_14_27_54_trim_decontam_read_count.tsv",
                            header = T)%>%
  as_tibble()
raw.read.counts<-raw.read.counts%>%
  filter(grepl("R1",file))%>%
  mutate(Sample=gsub("data\\/bowtie2_decontam_fastq\\/","",file),
         Sample=gsub("_trim_decontam","",Sample),
         Sample=gsub("\\.fastq\\.gz","",Sample),
         Sample_for_plot=gsub("_R1","",Sample))%>%
  mutate(counts_in_mln=num_seqs/10^6)%>%
  relocate(Sample,Sample_for_plot,num_seqs,counts_in_mln)

raw.read.count.plot<-ggplot(raw.read.counts,aes(x=Sample_for_plot,y=counts_in_mln))+
  geom_bar(stat="identity")+
  scale_y_continuous(limits=c(0, round(max(raw.read.counts$counts_in_mln),digits = -1)),
                     breaks=seq(0, round(max(raw.read.counts$counts_in_mln),digits = -1),5))+
  labs(y="Number of read pairs (millions)",
       x="Sample",
       title="Number of read pairs (millions) in the raw data \n(used to quantify gene abundances)")+
  theme_bw()+
  coord_cartesian(expand = FALSE)

ggsave(paste0(barplot.directory,
              paste(paste(format(Sys.time(),format="%Y%m%d"),
                          format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                    "raw-read-pairs-count-barplot",
                    sep = "-"),".png"),
       plot=raw.read.count.plot,
       width = 2000,height = 1500,
       units = "px",dpi=300,device = "png")

## Mapped reads ####
# Re-import data because I removed it previously
# gene counts
htseq.gene_counts.df<-read.table(htseq.gene_counts.fname,
                                 header = F,sep="\t")%>%
  as_tibble()
# 3,847,253 rows, like in prokka.contig.gene.map.df
colnames(htseq.gene_counts.df)<-c("locus_tag","raw_count")
# gene mappings to samples
prokka.contig.gene.map.df<-read.table(prokka.contig.gene.map.fname,
                                      header = F,sep="\t")%>%
  as_tibble()
colnames(prokka.contig.gene.map.df)<-c("sample","contig_id",
                                       "locus_tag","gene_length","gene_type")


mapped.read.counts<-htseq.gene_counts.df%>%
  left_join(prokka.contig.gene.map.df)%>%
  group_by(sample)%>%
  summarise(num_mapped_reads=sum(raw_count))

mapped.read.counts<-mapped.read.counts%>%
  mutate(counts_in_mln=num_mapped_reads/10^6)


mapped.read.count.plot<-mapped.read.counts%>%
  ggplot(aes(x=sample,y=counts_in_mln))+
  geom_bar(stat="identity")+
  scale_y_continuous(limits=c(0, 15))+
  labs(y="Number of reads (millions)",
       x="Sample",
       title="Number of reads (millions) mapped to genes by HTSeq")+
  theme_bw()+
  coord_cartesian(expand = FALSE)
ggsave(paste0(barplot.directory,
              paste(paste(format(Sys.time(),format="%Y%m%d"),
                          format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                    "mapped-read-count-barplot",
                    sep = "-"),".png"),
       plot=mapped.read.count.plot,
       width = 2000,height = 1500,
       units = "px",dpi=300,device = "png")

# raw.read.counts%>%
#   select(Sample_for_plot,num_seqs)%>%
#   left_join(mapped.read.counts,by=join_by("Sample_for_plot"=="sample"))%>%
#   pivot_longer(num_seqs:num_mapped_reads,names_to = "data_type")%>%
#   mutate(counts_in_mln=value/10^6)%>%
#   ggplot(aes(x=Sample_for_plot,y=counts_in_mln,fill=data_type))+
#   geom_bar(stat="identity",position="dodge")+
#   scale_y_continuous(limits=c(0, round(max(raw.read.counts$counts_in_mln),digits = -1)),
#                      breaks=seq(0, round(max(raw.read.counts$counts_in_mln),digits = -1),5))+
#   labs(y="Number of read pairs (millions)",
#        x="Sample",
#        title="Number of read pairs (millions) in the raw data \n(used to quantify gene abundances)")+
#   theme_bw()+
#   coord_cartesian(expand = FALSE)




#####
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
  full_join(gene.annotation.df)%>%
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


# Sanity check
mag.stats.with_annotation.df%>%arrange(-genome_size_mbp)%>%
  group_by(Phylum)%>%
  mutate(median_genome_size_mbp=median(genome_size_mbp))%>%
  ungroup%>%
  arrange(-median_genome_size_mbp)%>%
  mutate(Phylum=factor(Phylum,levels=unique(Phylum)))%>%
  ggplot(aes(x=Phylum,y=genome_size_mbp))+
  geom_boxplot()+
  geom_jitter()


write.table(mag.stats.with_annotation.df,
            file="./output/rtables/mag-stats-with-n_cds.tsv",
            sep = "\t", row.names = F)
