#' ---
#' output: 
#'   bookdown::html_document2:
#'      toc: true
#' ---

#' ```{r, setup 003-mag-taxonomy-plots.R, include=FALSE}
#' knitr::opts_knit$set(root.dir = '/home/rakhimov/projects/metagenome')
#' ```
#+ echo=FALSE
# MAG taxonomy barplots. ####
#' # MAG taxonomy barplots.
#'
#+ echo=FALSE
## Introduction ####
#'
#' ## Introduction
#' Plot taxonomy barplots from GTDB-tk and BLAST.
#' 
#+ echo=FALSE
## 1. Load necessary libraries. ####
#'
#' ## Load necessary libraries.
# install.packages(c("tidyverse","Polychrome","ggtext"))
library(tidyverse)
library(Polychrome)
library(ggtext)
#+ echo=FALSE
## 2. Import datasets. ####
#'
#' ## Import datasets.
#' Seqkit statistics on every contig and bin with taxonomy classification and
#' relative abudnance from CoverM:
seqkit.with_taxonomy<-readRDS(file = "./output/rdafiles/seqkit-with_taxonomy.rds")
#' Taxonomy classification of each bin from dRep
gtdbtk.coverm.drep<-readRDS(file = "./output/rdafiles/gtdbtk-coverm-drep.rds")
#' Taxonomy classification of each bin from MetaBAT2
gtdbtk.coverm.metabat2<-readRDS(file = "./output/rdafiles/gtdbtk-coverm-metabat2.rds")
#' Taxonomy classification of each contig from BLASTN
blastn.coverm.megahit<-readRDS(file = "./output/rdafiles/blastn-coverm-megahit.rds")
#' Tidy taxonomy table
gtdbtk.taxonomy<-readRDS("./output/rdafiles/gtdbtk-taxonomy.rds")
gene.annotation.df<-readRDS(file="./output/rdafiles/gene-annotation-df.rds")
#' Output paths
barplot.directory<-"./images/barplots/" 
boxplot.directory<-"./images/boxplots/" 
image.formats<-c("png","tiff")
#' Add CheckM2 quality information for all MAGs:
metabat2.date_time<-"20250612_13_37_47"
checkm2.high_quality_mags<-read.table("./output/mag_assembly/checkm2_output/20250619_05_47_09_high_quality_mags.tsv",
                                      sep="\t",header = T)%>%
  as_tibble()%>%
  rename("genome"="Name")%>%
  mutate(genome=gsub(paste0(metabat2.date_time,"_"),"",genome))%>%
  separate_wider_delim(genome,delim = "_bin.",names = c("sample","bin_id"))%>%
  mutate(bin_id=as.integer(bin_id))
#+ echo=FALSE
## 3. Taxonomy barplots (GTDB-tk) (**Figure 5**). ####
#'
#' ## Taxonomy barplots (GTDB-tk) (**Figure 5**).
#' Prepare the dRep classification dataset for plotting: hide rare taxa,
#' italicise names, highlight non-bacterial species in the legend.
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
#' Generate a color palette.
set.seed(1)
plot.cols<-createPalette(length(levels(unique(gtdbtk.coverm.drep.for_plot$taxon))),
                         seedcolors =rainbow(7))# input: number of rows

plot.cols<-c("#C1CDCD",plot.cols[1:length(plot.cols)-1])
col.vec<-setNames(plot.cols,levels(gtdbtk.coverm.drep.for_plot$taxon))
#' The barplot:
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
#+ fig.height=5.5, fig.width=11
print(gtdbtk.coverm.drep.plot + 
        ggtitle("Relative abundances of representative MAGs"))

#' MetaBAT2 classification plot wasn't included in the paper, so we leave it
#' as it is
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
#+ fig.height=7, fig.width=11
print(gtdbtk.coverm.metabat2.plot +
        ggtitle("Relative abundances of all high-quality MAGs (representative + redundant)"))

#+ echo=FALSE
## 4. Taxonomy barplots (BLASTN) (**Figure 5**). ####
#'
#' ## Taxonomy barplots (BLASTN) (**Figure 5**).
#' Prepare the dataset in a similar manner to previous plots:
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
#+ fig.height=5, fig.width=11
print(blastn.coverm.megahit.plot + 
        ggtitle("Classification of eukaryotic contigs with BLASTN"))
# for(image.format in image.formats){
#   ggsave(paste0(barplot.directory,
#                 paste(paste(format(Sys.time(),format="%Y%m%d"),
#                             format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                       "gtdbtk-coverm-drep-barplot",
#                       sep = "-"),".",image.format),
#          plot=gtdbtk.coverm.drep.plot,
#          width=11, height=5.5,units="in",
#          # width = 5000,height = 2200, units = "px",
#          dpi=300,device = image.format)
#   
#   ggsave(paste0(barplot.directory,
#                 paste(paste(format(Sys.time(),format="%Y%m%d"),
#                             format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                       "gtdbtk-coverm-metabat2-barplot",
#                       sep = "-"),".",image.format),
#          plot=gtdbtk.coverm.metabat2.plot,
#          width=11, height=7,units="in",
#          # width = 5000,height = 2200, units = "px",
#          dpi=300,device = image.format)
#   ggsave(paste0(barplot.directory,
#                 paste(paste(format(Sys.time(),format="%Y%m%d"),
#                             format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                       "blastn-coverm-megahit-barplot",
#                       sep = "-"),".",image.format),
#          plot=blastn.coverm.megahit.plot,
#          width=11, height=5,units="in",
#          # width = 5000,height = 2200,units = "px",
#          dpi=300,device = image.format)
# }

#+ echo=FALSE
## 5. Plot MAG genome sizes. ####
#'
#' ## Plot MAG genome sizes. 
mag.genome.sizes.plot<-seqkit.with_taxonomy%>%
  filter(is_drep_bin)%>%
  select(sample,contig_id,bin_id,contig_length,secondary_cluster,GC)%>%
  filter(!is.na(secondary_cluster))%>%
  group_by(secondary_cluster)%>%
  mutate(genome_size_mbp=sum(contig_length)/10^6,
         mag_gc=mean(GC),
         n_contig=n_distinct(contig_id))%>%
  ungroup%>%
  distinct(secondary_cluster,.keep_all = T)%>%
  left_join(gtdbtk.taxonomy)%>%
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
#+ fig.height=7, fig.width=10
print(mag.genome.sizes.plot + 
        ggtitle ("MAG genome sizes"))
# for(image.format in image.formats){
#   ggsave(filename = paste0("./images/boxplots/gtdbtk-mag-genome-sizes.",image.format),
#          dpi =300,
#          width=10, height=7,units="in",
#          # width = 4000, height = 2000,units="px",
#          mag.genome.sizes.plot,
#          device = image.format)
# }

#+ echo=FALSE
## 6. Plot MAG GC%. ####
#'
#' ## Plot MAG GC%. 
mag.gc.plot<-seqkit.with_taxonomy%>%
  filter(is_drep_bin)%>%
  select(sample,contig_id,bin_id,contig_length,secondary_cluster,GC)%>%
  filter(!is.na(secondary_cluster))%>%
  group_by(secondary_cluster)%>%
  mutate(genome_size_mbp=sum(contig_length)/10^6,
         mag_gc=mean(GC),
         n_contig=n_distinct(contig_id))%>%
  ungroup%>%
  distinct(secondary_cluster,.keep_all = T)%>%
  left_join(gtdbtk.taxonomy)%>%
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
  labs(y="GC%")+
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
#+ fig.height=7, fig.width=8
print(mag.gc.plot + 
        ggtitle("MAG GC%"))
# for(image.format in image.formats){
#   ggsave(filename = paste0("./images/boxplots/gtdbtk-mag-gc.",image.format),
#          dpi =300,
#          width=8, height=7,units="in",
#          mag.gc.plot,
#          device = image.format)
# }


#+ echo=FALSE
### 6.1 GC% in detail. ####
#'
#' ### GC% in detail.
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

#' GC% in BLASTN.
blastn.coverm.megahit.eukaryotes.gc.plot<-blastn.coverm.megahit%>%
  filter(grepl("Eukaryota",classification))%>%
  left_join(seqkit.with_taxonomy[,c("sample","contig_id","GC")])%>%
  # filter(Kingdom=="Eukaryota")%>%
  mutate(classification=sub(".*;[^;]*;([^;]*;[^;]*)$", "\\1", classification))%>%
  separate_wider_delim(delim = ";",cols=classification,
                       names=c("Genus","Species"))%>%
  select(sample,contig_id,GC,Genus)%>%
  arrange(-GC)%>%
  mutate(Genus=factor(Genus,levels=unique(Genus)))%>%
  ggplot(aes(x=Genus,y=GC))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()+
  labs(y="Contig GC %")+
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
#+ fig.height=6, fig.width=8
print(blastn.coverm.megahit.eukaryotes.gc.plot + 
        ggtitle("GC% of eukaryotic contigs identified by BLASTN"))
# for(image.format in image.formats){
#   ggsave(paste0(boxplot.directory,
#                 paste(paste(format(Sys.time(),format="%Y%m%d"),
#                             format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                       "blastn-coverm-megahit-eukaryotes-gc-boxplot",
#                       sep = "-"),".",image.format),
#          plot=blastn.coverm.megahit.eukaryotes.gc.plot,
#          width=8, height=6,units="in",
#          # width = 4500,height = 2000, units = "px",
#          dpi=300,device = image.format)
# }
sessionInfo()
rm(list = ls(all=TRUE))
gc()