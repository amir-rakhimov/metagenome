# Make figures wms ####
library(tidyverse)
library(Polychrome)
library(ggtext)
library(cowplot)
## Specifying parameters and directory/file names #### 
agglom.rank<-"Species" # this is the taxonomic rank that was used for agglomeration
barplot.directory<-"./images/barplots/" # set the path where barplots will
# be saved
image.formats<-c("png","tiff")
pipeline.name<-"kraken2"
rdafiles.directory<-"./output/rdafiles"
rtables.directory<-"./output/rtables"
metadata.directory<-"../amplicon_nmr/output/rdafiles" # path for custom.md.ages.rds
# Import abundance table from 002-phyloseq-kraken2.R
if(pipeline.name=="kraken2"){
  ps.q.agg.date_time<-"20241003_13_52_43"
}
ps.q.agg<-readRDS(file=file.path(
  rdafiles.directory,
  paste(ps.q.agg.date_time,"phyloseq",pipeline.name,agglom.rank,
        "table.rds",sep = "-")))
# Import metadata
custom.md<-readRDS(file.path(metadata.directory,"custom.md.ages.rds"))%>%
  filter(sequencing_type == "Naked mole-rat whole metagenome sequencing")

pretty.level.names<-c("NMR" = "*Heterocephalus glaber*") # better labels for facets)

custom.levels<-intersect(names(pretty.level.names),custom.md$class)
# Filter the phyloseq object to retain animal hosts from custom.levels
ps.q.agg<-ps.q.agg%>%
  filter(class%in%custom.levels,Abundance!=0)

# Figure 4 ####
# Barplot ####
agglom.rank.col<-which(colnames(ps.q.agg) =="Species")
preceding.rank.col<-which(colnames(ps.q.agg) ==agglom.rank)-1
preceding.rank<-colnames(ps.q.agg)[preceding.rank.col]
custom_order <- c("Kingdom", "Phylum", "Class", "Order", "Family","Genus","Species")
# If we agglomerate by higher level (Order,Class, etc), need to adjust the rank
if(agglom.rank%in%custom_order){
  agglom.rank.index<-match(agglom.rank,custom_order)
  custom_order<-custom_order[1:agglom.rank.index-1]
}

avg.taxon.len <- ps.q.agg%>%
  select(all_of(agglom.rank))%>%
  pull%>%
  nchar%>%
  mean%>%
  round

new.df<-ps.q.agg%>%
  unite("taxon",Kingdom:all_of(agglom.rank),sep = ";",remove = FALSE)%>%
  mutate(is_unclassified = grepl
         ("Kingdom|Phylum|Class|Order|Family|Genus|Species",taxon))%>% # add column to show that a taxon was unclassified
  mutate(is_unclassified=ifelse(is_unclassified==TRUE &!grepl
                                ("Kingdom|Phylum|Class|Order|Family|Genus|Species",
                                  get(agglom.rank)),
                                FALSE, is_unclassified))%>%
  mutate(
    taxon=ifelse(is_unclassified,
                 get(agglom.rank),
                 paste0(get(agglom.rank)#, " (",get(preceding.rank),
                        #")"
                 )),
    taxon = ifelse(mapply(grepl,get(preceding.rank),get(agglom.rank)),
                   taxon,
                   paste0(get(agglom.rank), " (", get(preceding.rank),")")),
    taxon=ifelse(MeanRelativeAbundance<1,
                 "Remainder (Mean relative abundance < 1%)",
                 taxon))%>% # if a taxon is classified, we change it to show the preceding rank, e.g Family (Genus)
  mutate(taxon=gsub("_"," ",taxon),
         taxon=str_wrap(taxon,width=avg.taxon.len))%>%
  mutate(taxon=ifelse(Kingdom=="Bacteria",
                      paste0("<span style='color: orange'><b><i>",taxon,"</i></b></span>"),
                      taxon),
         taxon=ifelse(Kingdom=="Eukaryota",
                      paste0("<span style='color: palegreen4'><b><i>",taxon,"</i></b></span>"),
                      taxon),
         taxon=ifelse(Kingdom=="Archaea",
                      paste0("<span style='color: violetred4'><b><i>",taxon,"</i></b></span>"),
                      taxon),
         taxon=ifelse(grepl("Remainder",taxon),
                      "Remainder (Mean relative abundance < 1%)",taxon))%>%
  group_by_at(c("is_unclassified",preceding.rank))%>%
  arrange(desc(is_unclassified),
          get(preceding.rank),
          taxon,
          .by_group = TRUE)%>% # sort within unclassified and classified taxa (by group)
  mutate(taxon=factor(taxon,levels=unique(taxon)),
         taxon=taxon%>%
           fct_relevel(c("Remainder (Mean relative abundance < 1%)",
                         unique(unlist(lapply(custom_order,
                                              function(pat){levels(.)[str_detect(levels(.),
                                                                                 fixed(pat))]
                                              })))
           ),
           after=0))%>% # change taxon column to factor. Then, we set "Remainder" as the first level, unclassified taxa next, then classified taxa. we use the str_detect to find unclassified taxa because these have Kingdom or Phylum in the name
  ungroup()

# Create a color palette
# We need to choose colors for the taxa in our barplot. They should be 
# distinguishable, so we can't choose similar colors. Or at least we shouldn't 
# put them next to each other.
# 
# We will use `createPalette` function from the `Polychrome` package.
# The function takes the number of colors for the palette and seed colors that 
# will be used for palette generation. In this case, we use the number of
# rows in our legend (taxa) and rainbow colors (to avoid having similar colors 
#                                               next to each other). 
# We also need to set the random seed because the output
# is a bit random. The output is a vector of colors.

set.seed(1)
plot.cols<-createPalette(length(levels(unique(new.df$taxon))),
                         seedcolors =rainbow(7))# input: number of rows

plot.cols<-c("#C1CDCD",plot.cols[1:length(plot.cols)-1])
col.vec<-setNames(plot.cols,levels(new.df$taxon))

## Plot the barplot ####
mainplot<-new.df%>%
  left_join(custom.md[,c("Sample","age","host_birthday")])%>%
  arrange(age,host_birthday)%>%
  mutate(age = ifelse(age==1,paste0("(",age,")"),
                      paste0("(",age,")")),
         age = factor(age),
         Sample = paste(Sample,age),
         Sample = factor(Sample,levels=unique(Sample)))%>%
  
  ggplot(aes(x=Sample,y=RelativeAbundance,fill=taxon))+
  geom_bar(stat="identity")+
  scale_fill_manual(values = col.vec)+ # custom fill that is based on our 
  coord_cartesian(expand=FALSE) +
  # custom palette
  xlab("") +
  ylab("Relative Abundance (%)")+
  labs(fill="Taxon")+
  theme_bw()+
  guides(fill=guide_legend(ncol=1))+ # legend as one column
  # ggtitle(paste(agglom.rank,"level gut microbiota profiles of fecal samples from naked mole-rats"))+
  theme(plot.margin=unit(c(0.2,0,0.1,0.1), 'cm'),
        # axis.line = element_blank(), #TODO: what does it do?
        strip.text.x = ggtext::element_markdown(size = 10),# the name of
        # each facet will be recognised as a markdown object, so we can
        # add line breaks (cause host names are too long)
        panel.spacing = unit(0.8, "cm"), # increase distance between facets
        axis.text.x = element_text(angle=45,size=9,hjust=1),# rotate
        # the x-axis labels by 45 degrees and shift to the right
        axis.text.y = element_text(size=10), # size of y axis ticks
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12), # size of axis names
        plot.title = element_text(size = 12), # size of plot title
        plot.caption = element_text(size=12), # size of plot caption
        legend.box.spacing = unit(0.01, units="cm"),
        legend.key.size = unit(0.5, units= "cm"),
        legend.text = element_markdown(size = 9), # size of legend text
        legend.title = element_text(size = 9), # size of legend title
        legend.position = "right") # legend on the right

# Relative abundance per kingdom (including bacteria) ####
if(pipeline.name=="kraken2"){
  unclassified_reads<-read.table("./output/kraken2_pipeline/20240409_17_32_40_unclassified_reads.tsv",
                                 header = T)
  colnames(unclassified_reads)[which(colnames(unclassified_reads)=="Unclassified")]<-"Abundance"
  unclassified_reads$Kingdom<-"Unclassified"
}

kingdom.plot.with_bact<-ps.q.agg%>%
  group_by(Sample,Kingdom)%>%
  summarise(Abundance=sum(Abundance))%>%
  bind_rows(unclassified_reads)%>%
  group_by(Kingdom)%>%
  summarise(TotalAbundance=sum(Abundance))%>%
  mutate(TotalAbundancePerTax=TotalAbundance/sum(TotalAbundance)*100)%>%
  arrange(-TotalAbundancePerTax)%>%
  ggplot(aes(x=reorder(Kingdom,-TotalAbundancePerTax),
             y=TotalAbundancePerTax,
             fill=Kingdom))+
  geom_bar(stat = "identity")+
  scale_fill_viridis_d(option = "C")+
  ylab("Relative abundance (%)")+
  xlab("")+
  theme_bw()+
  # ggtitle("Whole metagenome sequencing profile of naked mole-rat gut microbiota
  #        (Kingdom level)")+
  coord_cartesian(expand = c("bottom" = F,
                             "left" = F,
                             "right" = F,
                             "top" = T),
                  ylim = c(0,80))+
  geom_text(aes(x=Kingdom, 
                y=TotalAbundancePerTax, 
                label=round(TotalAbundancePerTax,digits = 3)),
            vjust=-0.5,size=2)+ # add text with percentage above the bar
  theme(plot.margin=unit(c(0.2,0,0.1,0.1), 'cm'),
        axis.line = element_blank(), 
        axis.text.x = element_text(size=5),
        axis.text.y = element_text(size=5), # size of y axis ticks
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 7), # size of axis names
        plot.title = ggtext::element_markdown(size = 25), # the plot 
        # title will be recognised as a markdown object, so we can
        # add line breaks (cause host names are too long)
        plot.caption = element_text(size=23),# size of plot caption
        legend.text = element_text(size = 10),# size of legend text
        legend.title = element_text(size = 12), # size of legend title
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) # no legend 


# Save figure 4 ####
# figure4 <-plot_grid(kingdom.plot.with_bact,mainplot, 
#                     ncol=1)

figure4<-ggdraw(mainplot)+
  draw_plot(kingdom.plot.with_bact, 
            x = 0.02,
            y = 0.62,
            width = 0.85,
            height = 0.5,
            scale=0.4,
            hjust = 0.17)#+
  # draw_plot_label(
  #   # c("b", "a"),
  #   c(0.02, 0.12),
  #   c(1, 0.98),
  #   size = 20
  # )

figure4
# save the plots: kingdom.plot.with_bact
for(image.format in image.formats){
  ggsave(filename = paste(paste(format(Sys.time(),format="%Y%m%d"),
                                format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                          paste("figure4",image.format,sep = "."),sep = "-"),
         path = "images",
         plot = figure4,
         device = image.format,
         dpi=300,
         width = 170,
         height = 130, #190
         units = "mm",
         scale = 1)
}

# Figure 5 ####
seqkit.with_taxonomy<-readRDS(file = "./output/rdafiles/seqkit-with_taxonomy.rds")
gtdbtk.coverm.drep<-readRDS(file = "./output/rdafiles/gtdbtk-coverm-drep.rds")
gtdbtk.coverm.metabat2<-readRDS(file = "./output/rdafiles/gtdbtk-coverm-metabat2.rds")
blastn.coverm.megahit<-readRDS(file = "./output/rdafiles/blastn-coverm-megahit.rds")
gtdbtk.taxonomy<-readRDS("./output/rdafiles/gtdbtk-taxonomy.rds")
gene.annotation.df<-readRDS(file="./output/rdafiles/gene-annotation-df.rds")

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
                           Species == "Fimivicinus sp900544375"~ paste0("<i>",Genus,"</i>",  " (", "<i>",Family,"</i>",")"),
                           Genus == "Phil12" ~ paste0(Genus,  " (", "<i>",Order,"</i>",")"),
                           grepl("UBA",Genus)~ paste0(Genus,  " (", "<i>",Family,"</i>",")"),
                           grepl("MGBC",Genus)~ paste0(Genus,  " (", "<i>",Family,"</i>",")"),
                           Species =="unclassified" & Genus!="unclassified" ~ paste0("<i>",Genus,"</i>",  " (", "<i>",Family,"</i>",")"),
                           Species != "unclassified" & Genus!="unclassified" ~ Species,
                           Species =="unclassified"  & Genus=="unclassified" ~ paste0("<i>",Family,"</i>")
  # mutate(taxon = case_when(mean_relab < 0.7 ~ "Remainder (mean relative abundance < 0.7%)",
  #                          Species == "Fimivicinus sp900544375"~ paste0(Genus,  " (", Family,")"),
  #                          Genus == "Phil12" ~ paste0(Genus,  " (", Order,")"),
  #                          Species =="unclassified" & Genus!="unclassified" ~ paste0(Genus, " (", Family,")"),
  #                          Species != "unclassified" & Genus!="unclassified" ~ Species,
  #                          Species =="unclassified"  & Genus=="unclassified" ~ Family
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
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=10), # size of y axis ticks
        axis.title = element_text(size = 11), # size of axis names
        plot.caption = element_text(size=10), # size of plot caption
        legend.title = element_text(size = 12) # size of legend title
  )+
  theme(plot.title = element_text(size=10), 
        panel.spacing = unit(0.8, "cm"), # increase distance between facets
        legend.text = element_markdown(size = 10), # size of legend text
        legend.position = "right",
        legend.key.size = unit(0.4, units= "cm"),
        legend.box.spacing = unit(0.01, units="cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

#
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
                        paste0("<i>",Genus,"</i>",  " (", "<i>",Family,"</i>",")"),
                        paste0("<i>",Species,"</i>", " (", "<i>",Genus,"</i>",")")),
         taxon = str_wrap(taxon,width=avg.taxon.len))%>%
  mutate(is_unclassified = ifelse(grepl("sp\\.", Species),
                                  TRUE,
                                  FALSE))%>%
  arrange(is_unclassified,Species)%>%
  mutate(taxon = factor(taxon, levels = unique(taxon)))%>%
  select(-is_unclassified)
  # mutate(taxon = ifelse(grepl("sp\\.", Species),
  #                       paste0(Genus, " (", Family,")"),
  #                       paste0("<i>",Species,"</i>", " (", Genus,")")),
  #        taxon = str_wrap(taxon,width=avg.taxon.len))

set.seed(1)
blastn.plot.cols<-createPalette(length(levels(blastn.coverm.megahit.for_plot$taxon)),
                                seedcolors =rainbow(7))# input: number of rows

# blastn.plot.cols<-c("#C1CDCD",blastn.plot.cols[1:length(plot.cols)-1])
blastn.col.vec<-setNames(blastn.plot.cols,levels(blastn.coverm.megahit.for_plot$taxon))

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
  theme(axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size=10), # size of y axis ticks
        axis.title = element_text(size = 11), # size of axis names
        plot.title = element_text(size = 12), # size of plot title
        plot.caption = element_text(size=10), # size of plot caption
        legend.text = element_markdown(size = 10), # size of legend text
        legend.title = element_text(size = 12), # size of legend title
        legend.box.spacing = unit(0.01, units="cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "right") # legend on the right

# 
# ko.defs<-readRDS(file="./output/rdafiles/ko-defs.rds")
# representative.gene_lengths<-readRDS(file="./output/rdafiles/representative-gene_lengths.rds")
# cazymes.blast.df<-read.table(file="./output/mag_assembly/blastp_output/20250802_12_29_03/20250802_12_29_03_blastp_cazy_db_cazymes_top_hits_with_taxa.tsv",
#                              sep = "\t",header = F,fill=T)%>%
#   as_tibble()
# colnames(cazymes.blast.df)<-c("qacc","sacc", "qseqid","sseqid", 
#                               "pident", "length", "mismatch", "gapopen", 
#                               "qstart", "qend", "sstart", "send",
#                               "evalue", "bitscore", "sframe","uniparc",
#                               "caz_db","Kingdom","Species")
# cazymes.blast.df<-cazymes.blast.df%>%
#   mutate(caz_db=gsub(".*fasta\\|","",caz_db))%>%
#   select(-qseqid,-sseqid)%>%
#   rename("locus_tag_cluster"="qacc")%>%
#   mutate(Kingdom=ifelse(Kingdom=="","unmapped",Kingdom),
#          Species=ifelse(Kingdom=="unmapped","unmapped",Species))
# # Number of unique genes by annotation in total (kofamscan, dbcan, all) (**My Figure 2**) ####
# # https://www.nature.com/articles/s42003-021-02827-2
# uniq.gene.counts.all<-gene.annotation.df%>%
#   distinct(locus_tag_cluster,.keep_all = T)%>%
#   count(annotation_type,sort=TRUE)
# 
# ### The barplot of all genes for each annotation (4 bars) ####
# uniq.gene.counts.all.bp<-uniq.gene.counts.all%>%
#   mutate(annotation_type=factor(annotation_type,
#                                 levels=sort(uniq.gene.counts.all$annotation_type)))%>%
#   ggplot(aes(x=annotation_type,y=n,fill=annotation_type))+
#   geom_bar(stat="identity")+
#   geom_text(aes(label=n),vjust = -0.3, nudge_y = -0.5,
#             size=8)+ # add counts of genes
#   scale_fill_viridis_d(option = "C")+
#   labs(y="Number of unique genes",x="",
#        fill="Annotation type")+
#   theme_bw()+
#   coord_cartesian(ylim = c(0,max(uniq.gene.counts.all$n)+200000),
#                   expand = F)+
#   theme(axis.text.x = element_blank(),
#         axis.text.y = element_text(size=20), # size of y axis ticks
#         axis.title = element_text(size = 20), # size of axis names
#         plot.title = element_text(size = 25), # size of plot title
#         plot.caption = element_text(size=23), # size of plot caption
#         legend.text = element_text(size = 20), # size of legend text
#         legend.title = element_text(size = 25), # size of legend title
#         legend.position = "inside", # legend inside the plot
#         legend.position.inside=c(0.15,0.8),
#         panel.grid.minor = element_blank(),
#         panel.grid.major = element_blank()) # legend on the left top


# Save figures ####
figure5<-plot_grid(gtdbtk.coverm.drep.plot,
                   blastn.coverm.megahit.plot, 
                   # uniq.gene.counts.all.bp,
                   labels = "auto",
                   label_size = 20,
                   ncol=1)
# save the plots: kingdom.plot.with_bact
for(image.format in image.formats){
  ggsave(filename = paste(paste(format(Sys.time(),format="%Y%m%d"),
                                format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                          paste("figure5",image.format,sep = "."),sep = "-"),
         path = "images",
         plot = figure5,
         device = image.format,
         dpi=300,
         width = 170,
         height = 180, # 190,
         units = "mm",
         scale = 1# 2.2
         )
}
