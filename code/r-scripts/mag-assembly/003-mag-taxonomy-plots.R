library(tidyverse)
library(Polychrome)
library(ggtext)

seqkit.with_taxonomy<-readRDS(file = "./output/rdafiles/seqkit-with_taxonomy.rds")
gtdbtk.coverm.drep<-readRDS(file = "./output/rdafiles/gtdbtk-coverm-drep.rds")
gtdbtk.coverm.metabat2<-readRDS(file = "./output/rdafiles/gtdbtk-coverm-metabat2.rds")
blastn.coverm.megahit<-readRDS(file = "./output/rdafiles/blastn-coverm-megahit.rds")
gtdbtk.taxonomy<-readRDS("./output/rdafiles/gtdbtk-taxonomy.rds")
gene.annotation.df<-readRDS(file="./output/rdafiles/gene-annotation-df.rds")

barplot.directory<-"./images/barplots/" # set the path where barplots will
boxplot.directory<-"./images/boxplots/" # set the path where boxplots will
image.formats<-c("png","tiff")

metabat2.date_time<-"20250612_13_37_47"
checkm2.high_quality_mags<-read.table("./output/mag_assembly/checkm2_output/20250619_05_47_09_high_quality_mags.tsv",
                                      sep="\t",header = T)%>%
  as_tibble()%>%
  rename("genome"="Name")%>%
  mutate(genome=gsub(paste0(metabat2.date_time,"_"),"",genome))%>%
  separate_wider_delim(genome,delim = "_bin.",names = c("sample","bin_id"))%>%
  mutate(bin_id=as.integer(bin_id))

# Basic taxonomy barplot (GTDB-tk) (**Figure 5**) ####
avg.taxon.len <- round(mean(nchar(gtdbtk.coverm.drep$classification)))

gtdbtk.coverm.drep.for_plot<-gtdbtk.coverm.drep%>%
  group_by(classification)%>%
  mutate(mean_relab=mean(relative_abundance))%>%
  ungroup()%>%
  separate_wider_delim(cols = classification,
                       delim = ";",
                       names = c("Kingdom","Phylum","Class","Order",
                                 "Family","Genus","Species"),
                       cols_remove = FALSE)%>%
  mutate(taxon = case_when(mean_relab < 0.7 ~ "Remainder (mean relative abundance < 0.7%)",
                           Species == "Fimivicinus sp900544375"~ paste0(Genus,  " (", Family,")"),
                           Genus == "Phil12" ~ paste0(Genus,  " (", Order,")"),
                           Species =="unclassified" & Genus!="unclassified" ~ paste0(Genus, " (", Family,")"),
                           Species != "unclassified" & Genus!="unclassified" ~ Species,
                           Species =="unclassified"  & Genus=="unclassified" ~ Family
                           ),
         taxon = ifelse(! grepl("\\(", taxon) & !grepl("Remainder",taxon)& grepl(" ",taxon), 
                        paste0("<i>",taxon,"</i>"), taxon),
         taxon = case_when(Kingdom =="Archaea" & !grepl("Remainder",taxon) ~ paste0("<span style='color: violetred4'><b>",taxon,"</b></span>"),
                           Kingdom =="Bacteria" & !grepl("Remainder",taxon) ~ paste0("<span style='color: orange'><b>",taxon,"</b></span>"),
                           grepl("Remainder", taxon) ~ taxon
                           ))%>%
  arrange(Family)%>%
  mutate(taxon=factor(taxon,levels=unique(taxon)))%>%
  mutate(taxon=fct_relevel(taxon,
                                    "Remainder (mean relative abundance < 0.7%)",
                                    after = 0))
# Color palette is used for all plots
set.seed(1)
plot.cols<-createPalette(length(levels(unique(gtdbtk.coverm.drep.for_plot$taxon))),
                         seedcolors =rainbow(7))# input: number of rows

plot.cols<-c("#C1CDCD",plot.cols[1:length(plot.cols)-1])
col.vec<-setNames(plot.cols,levels(gtdbtk.coverm.drep.for_plot$taxon))

gtdbtk.coverm.drep.plot<-gtdbtk.coverm.drep.for_plot%>%
  ggplot(aes(x=sample,y=relative_abundance,fill=taxon))+
  geom_bar(stat="identity")+
  scale_fill_manual(values = col.vec)+ # custom fill that is based on our 
  guides(fill=guide_legend(ncol=1))+
  theme(legend.position = "right")+
  theme_bw()+
  coord_cartesian(expand=FALSE) +
  labs(x="",
       y="Relative Abundance (%)",
       # title="Relative abundances of representative MAGs",
       fill="Classification"
       )+
  theme(axis.text.x = element_text(size=13),
        axis.text.y = element_text(size=13), # size of y axis ticks
        axis.title = element_text(size = 15), # size of axis names
        plot.caption = element_text(size=13), # size of plot caption
        legend.title = element_text(size = 15) # size of legend title
  )+
  theme(plot.title = element_text(size=25), 
        panel.spacing = unit(0.8, "cm"), # increase distance between facets
        legend.text = element_markdown(size = 13), # size of legend text
        legend.position = "right",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

#

gtdbtk.coverm.metabat2.for_plot<-gtdbtk.coverm.metabat2%>%
  group_by(classification)%>%
  mutate(mean_relab=mean(relative_abundance))%>%
  ungroup()%>%
  # mutate(classification=ifelse(mean_relab<1,"Remainder",classification))%>%
  mutate(classification=ifelse(mean_relab<0.7,
                               "Remainder (mean relative abundance < 0.7%)",
                               classification))%>%
  mutate(classification=gsub(";","; ",classification),
         classification=str_wrap(classification,width=avg.taxon.len),
         classification=gsub("; ",";",classification))%>%
  arrange(classification)%>%
  mutate(classification=factor(classification,levels=unique(classification)))%>%
  mutate(classification=fct_relevel(classification,
                                    "Remainder (mean relative abundance < 0.7%)",
                                    after = 0))

set.seed(1)
plot.cols<-createPalette(length(levels(unique(gtdbtk.coverm.metabat2.for_plot$classification))),
                         seedcolors =rainbow(7))# input: number of rows

plot.cols<-c("#C1CDCD",plot.cols[1:length(plot.cols)-1])
col.vec<-setNames(plot.cols,levels(gtdbtk.coverm.metabat2.for_plot$classification))

gtdbtk.coverm.metabat2.plot<-gtdbtk.coverm.metabat2.for_plot%>%
  ggplot(aes(x=sample,y=relative_abundance,fill=classification))+
  geom_bar(stat="identity")+
  scale_fill_manual(values = col.vec)+ # custom fill that is based on our 
  guides(fill=guide_legend(ncol=1))+
  theme_bw()+
  theme(legend.position = "right")+
  # labs(title = Sys.time())+
  coord_cartesian(expand=FALSE) +
  labs(x="",
       y="Relative Abundance (%)",
       # title="Relative abundances of all high-quality MAGs (representative + redundant)",
       fill="Classification")+
  theme(#plot.margin=unit(c(1,1,1,1.5), 'cm'),
        axis.text.x = element_text(size=13),
        axis.text.y = element_text(size=13), # size of y axis ticks
        axis.title = element_text(size = 15), # size of axis names
        plot.caption = element_text(size=13), # size of plot caption
        legend.title = element_text(size = 15) # size of legend title
  )+
  theme(plot.title = element_text(size=25), 
        panel.spacing = unit(0.8, "cm"), # increase distance between facets
        legend.text = element_text(size = 8), # size of legend text
        legend.position = "right",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) 

## Taxonomy barplot BLASTN (**Figure 5b**) ####
blastn.coverm.megahit.for_plot<-blastn.coverm.megahit%>%
  filter(grepl("Eukaryota",classification))%>%
  mutate(classification = str_extract(classification, "[^;]+;[^;]+;[^;]+$"))%>%
  separate_wider_delim(classification,
                       delim = ";",
                       names = c("Family","Genus", "Species"),
                       cols_remove = F)%>%
  mutate(Species = gsub(" \\(sweet potato\\)","",Species))%>%
  mutate(taxon = ifelse(grepl("sp\\.", Species),
                        paste0(Genus, " (", Family,")"),
                        paste0("<i>",Species,"</i>", " (", Genus,")")),
         taxon = str_wrap(taxon,width=avg.taxon.len))

set.seed(1)
blastn.plot.cols<-createPalette(length(levels(unique(blastn.coverm.megahit.for_plot$classification))),
                         seedcolors =rainbow(7))# input: number of rows

# blastn.plot.cols<-c("#C1CDCD",blastn.plot.cols[1:length(plot.cols)-1])
blastn.col.vec<-setNames(blastn.plot.cols,levels(blastn.coverm.megahit.for_plot$classification))

blastn.coverm.megahit.plot<-blastn.coverm.megahit.for_plot%>%
  ggplot(aes(x=sample,y=tpm,fill=taxon))+
  geom_bar(stat="identity")+
  coord_cartesian(expand = F)+
  labs(x="Sample",
       y="TPM",
       # title="TPM of eukaryotic contigs",
       fill="Classification")+
  theme_bw()+
  scale_fill_manual(values = blastn.col.vec)+ # custom fill that is based on our 
  theme(axis.text.x = element_text(size=13),
        axis.text.y = element_text(size=13), # size of y axis ticks
        axis.title = element_text(size = 15), # size of axis names
        plot.title = element_text(size = 25), # size of plot title
        plot.caption = element_text(size=13), # size of plot caption
        legend.text = element_markdown(size = 13), # size of legend text
        legend.title = element_text(size = 15), # size of legend title
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "right") # legend on the right

for(image.format in image.formats){
  ggsave(paste0(barplot.directory,
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "gtdbtk-coverm-drep-barplot",
                      sep = "-"),".",image.format),
         plot=gtdbtk.coverm.drep.plot,
         width=11, height=5.5,units="in",
         # width = 5000,height = 2200, units = "px",
         dpi=300,device = image.format)
  
  ggsave(paste0(barplot.directory,
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "gtdbtk-coverm-metabat2-barplot",
                      sep = "-"),".",image.format),
         plot=gtdbtk.coverm.metabat2.plot,
         width=11, height=7,units="in",
         # width = 5000,height = 2200, units = "px",
         dpi=300,device = image.format)
  ggsave(paste0(barplot.directory,
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "blastn-coverm-megahit-barplot",
                      sep = "-"),".",image.format),
         plot=blastn.coverm.megahit.plot,
         width=11, height=5,units="in",
         # width = 5000,height = 2200,units = "px",
         dpi=300,device = image.format)
}


## Taxa with most unique genes ####
seqkit.with_taxonomy.with_genes%>%
  distinct(locus_tag_cluster,.keep_all = T)%>%
  add_count(sample,bin_id,name = "n_genes")%>% # count locus_tag_cluster
  ungroup%>%
  select(sample,bin_id,n_genes,gtdbtk_result)%>%
  distinct()%>%
  arrange(-n_genes)

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


write.table(mag.stats.with_annotation.df,
            file="./output/rtables/mag-stats-with-n_cds.tsv",
            sep = "\t", row.names = F)



# MAG genome sizes ####
mag.genome.sizes.plot<-mag.stats.with_annotation.df%>%
  group_by(Phylum)%>%
  mutate(median_genome_size_mbp=median(genome_size_mbp))%>%
  ungroup%>%
  arrange(-median_genome_size_mbp)%>%
  select(-median_genome_size_mbp)%>%
  mutate(Phylum=factor(Phylum,levels=unique(Phylum)))%>%
  ggplot(aes(x=Phylum,y=genome_size_mbp))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()+
  labs(y="Genome size (Mbp)"#,
       # title="MAG genome sizes"
       )+
  theme(axis.text.x = element_text(angle=45,size=14,hjust=1),# rotate 
        # the x-axis labels by 45 degrees and shift to the right
        axis.text.y = element_text(size=15), # size of y axis ticks
        axis.title = element_text(size = 14), # size of axis names
        plot.title = element_text(size=10), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.caption = element_text(size=13), # size of plot caption
        legend.title = element_text(size = 15) # size of legend title
  )
for(image.format in image.formats){
  ggsave(filename = paste0("./images/boxplots/gtdbtk-mag-genome-sizes.",image.format),
         dpi =300,
         width=10, height=7,units="in",
         # width = 4000, height = 2000,units="px",
         mag.genome.sizes.plot,
         device = image.format)
}

# GC% in MAGs (**My Figure 3**) ####
## Boxplot ####
mag.gc.plot<-mag.stats.with_annotation.df%>%
  group_by(Phylum)%>%
  mutate(median_gc=median(mag_gc))%>%
  ungroup%>%
  arrange(-median_gc)%>%
  select(-median_gc)%>%
  mutate(Phylum=factor(Phylum,levels=unique(Phylum)))%>%
  ggplot(aes(x=Phylum,y=mag_gc))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()+
  labs(y="GC%"#,
       # title="MAG GC%",
       )+
  theme(
    plot.margin = unit(c(0.2,0.2,0.2,0.5), "cm"),
        axis.text.x = element_text(angle=45,size=14,hjust=1),# rotate 
        # the x-axis labels by 45 degrees and shift to the right
        axis.text.y = element_text(size=13), # size of y axis ticks
        axis.title = element_text(size = 14), # size of axis names
        plot.title = element_text(size=10), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.caption = element_text(size=13), # size of plot caption
        legend.title = element_text(size = 15) # size of legend title
  )
for(image.format in image.formats){
  ggsave(filename = paste0("./images/boxplots/gtdbtk-mag-gc.",image.format),
         dpi =300,
         width=8, height=7,units="in",
         # width = 4000, height = 2000, units="px",
         mag.gc.plot,
         device = image.format)
}


## GC stats ####
seqkit.with_taxonomy%>%
  filter(is_drep_bin)%>%
  select(sample,contig_id,bin_id,contig_length,secondary_cluster,GC)%>%
  left_join(gtdbtk.taxonomy)%>%
  filter(!is.na(secondary_cluster))%>%
  group_by(secondary_cluster)%>%
  mutate(mag_gc=mean(GC))%>%
  ungroup%>%
  distinct(secondary_cluster,.keep_all = T)%>%
  group_by(Phylum)%>%
  mutate(min_gc=min(mag_gc),
         max_gc=max(mag_gc),
         mean_gc=mean(mag_gc),
         sd_gc=sd(mag_gc),
         n=n())%>%
  ungroup%>%
  distinct(Phylum,min_gc,max_gc,mean_gc,sd_gc,n)%>%
  arrange(-mean_gc)



# GC% in BLASTN ####
blastn.coverm.megahit.eukaryotes.gc.plot<-blastn.coverm.megahit%>%
  filter(grepl("Eukaryota",classification))%>%
  left_join(seqkit.with_taxonomy[,c("sample","contig_id","GC")])%>%
  # filter(Kingdom=="Eukaryota")%>%
  mutate(classification=sub(".*;[^;]*;([^;]*;[^;]*)$", "\\1", classification))%>%
  separate_wider_delim(delim = ";",cols=classification,
                       names=c("Genus","Species"))%>%
  # group_by(Genus)%>%
  select(sample,contig_id,GC,Genus)%>%
  arrange(-GC)%>%
  mutate(Genus=factor(Genus,levels=unique(Genus)))%>%
  ggplot(aes(x=Genus,y=GC))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()+
  labs(y="Contig GC %"#,
       # title="GC% of eukaryotic contigs identified by BLASTN"
       )+
  theme(axis.text.x = element_text(angle=45,size=12,hjust=1),
        axis.text.y = element_text(size=12), # size of y axis ticks
        axis.title = element_text(size = 14), # size of axis names
        plot.title = element_text(size = 10), # size of plot title
        plot.caption = element_text(size=10), # size of plot caption
        legend.text = element_text(size = 10), # size of legend text
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.title = element_text(size = 10), # size of legend title
        legend.position = "right") # legend on the right

for(image.format in image.formats){
  ggsave(paste0(boxplot.directory,
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "blastn-coverm-megahit-eukaryotes-gc-boxplot",
                      sep = "-"),".",image.format),
         plot=blastn.coverm.megahit.eukaryotes.gc.plot,
         width=8, height=6,units="in",
         # width = 4500,height = 2000, units = "px",
         dpi=300,device = image.format)
}


