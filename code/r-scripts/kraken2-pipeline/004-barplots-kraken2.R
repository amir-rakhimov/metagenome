#' ---
#' output: 
#'   bookdown::html_document2:
#'      toc: true
#' ---

#' ```{r, setup 004-barplots-kraken2.R, include=FALSE}
#' knitr::opts_knit$set(root.dir = '/home/rakhimov/projects/metagenome')
#' ```
#+ echo=FALSE
# Barplots ####
#' # Barplots
#'
#+ echo=FALSE
## Introduction ####
#'
#' ## Introduction
#' After importing the kraken2 output files into R and processing them 
#' with phyloseq, it's time to explore the taxonomic composition of our data.
#' We will use the Polychrome package to create a custom palette for the 
#' barplots.

#' ```
#+ echo=FALSE
## 1. Load necessary libraries. ####
#'
#' ## Load necessary libraries.
# install.packages(c("tidyverse","ggtext","Polychrome"))
library(tidyverse)
library(Polychrome)
library(ggtext)
#+ echo=FALSE
## 2. Specifying parameters and directory/file names. #### 
#'
#' ## Specifying parameters and directory/file names. 
#' The taxonomic rank that was used for agglomeration:
agglom.rank<-"Species" 
#' Specify paths and image formats:
barplot.directory<-"./images/barplots/"
rdafiles.directory<-"./output/rdafiles"
rtables.directory<-"./output/rtables"
image.formats<-c("png","tiff")
#' Path for custom metadata:
metadata.directory<-"../amplicon_nmr/output/rdafiles" # path for custom.md.ages.rds
#' Import abundance table as an rds file (NOT rarefied): 
ps.q.agg.date_time<-"20241003_13_52_43"
ps.q.agg<-readRDS(file=file.path(
  rdafiles.directory,
  paste(ps.q.agg.date_time,"phyloseq","kraken2",agglom.rank,
        "table.rds",sep = "-")))
#' Import metadata:
custom.md<-readRDS(file.path(metadata.directory,"custom.md.ages.rds"))%>%
  filter(sequencing_type == "Naked mole-rat whole metagenome sequencing")

#' Check the total number of taxa (not filtered by mean relative abundance)
ps.q.agg%>%
  ungroup()%>%
  select(matches(paste0("^",agglom.rank,"$")))%>%
  pull(.)%>%
  unique()%>%
  sort()%>%
  length()
#' select(matches(paste0("^",agglom.rank,"$"))) will select variables that 
#' match a pattern (regex).
#' So, if we are agglomerating by Phylum, the command will select the Phylum column

#+ echo=FALSE
## 3. Set up the barplot. #### 
#'
#' ## Set up the barplot. 
#' Find the column number that precedes Species (Genus column).
if(agglom.rank=="Species"){
  agglom.rank.col<-which(colnames(ps.q.agg) =="Species")
  preceding.rank.col<-which(colnames(ps.q.agg) ==agglom.rank)-1
  preceding.rank<-colnames(ps.q.agg)[preceding.rank.col]
}else{
  agglom.rank.col<-which(colnames(ps.q.agg) ==agglom.rank)
  preceding.rank.col<-which(colnames(ps.q.agg) ==agglom.rank)-1
  preceding.rank<-colnames(ps.q.agg)[preceding.rank.col]
}
custom_order <- c("Kingdom", "Phylum", "Class", "Order", "Family","Genus","Species")
#' If we agglomerate by higher level (Order,Class, etc), need to adjust the rank
if(agglom.rank%in%custom_order){
  agglom.rank.index<-match(agglom.rank,custom_order)
  custom_order<-custom_order[1:agglom.rank.index-1]
}
#' Find the average number of characters in the taxa to wrap long names:
avg.taxon.len <- ps.q.agg%>%
  select(all_of(agglom.rank))%>%
  pull%>%
  nchar%>%
  mean%>%
  round
#' Hide rare taxa, wrap long names, highlight non-bacterial taxa.
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

#' Create a color palette.  
#' We need to choose colors for the taxa in our barplot. They should be
#' distinguishable, so we can't choose similar colors. Or at least we shouldn't
#' put them next to each other.
#' 
#' We will use `createPalette` function from the `Polychrome` package.
#' The function takes the number of colors for the palette and seed colors that
#' will be used for palette generation. In this case, we use the number of
#' rows in our legend (taxa) and rainbow colors (to avoid having similar colors
#' next to each other).
#' We also need to set the random seed because the output
#' is a bit random. The output is a vector of colors.

set.seed(1)
plot.cols<-createPalette(length(levels(unique(new.df$taxon))),
                         seedcolors =rainbow(7))# input: number of rows
#' Grey color is added ahead of others because it corresponds to Remainder 
#' taxa. 
plot.cols<-c("#C1CDCD",plot.cols[1:length(plot.cols)-1])
col.vec<-setNames(plot.cols,levels(new.df$taxon))

#+ echo=FALSE
## 3. Plot the barplot. ####
#'
#' ## Plot the barplot. 
mainplot<-new.df%>%
  left_join(custom.md[,c("Sample","age","host_birthday")])%>%
  arrange(age,host_birthday)%>%
  mutate(age = ifelse(age==1,paste0("(",age," year old)"),
                      paste0("(",age," years old)")),
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
  theme(plot.margin=unit(c(0.4,0.4,0.4,1), 'cm'),
        # axis.line = element_blank(), #TODO: what does it do?
        strip.text.x = ggtext::element_markdown(size = 10),# the name of
        # each facet will be recognised as a markdown object, so we can
        # add line breaks (cause host names are too long)
        panel.spacing = unit(0.8, "cm"), # increase distance between facets
        axis.text.x = element_text(angle=45,size=13,hjust=1),# rotate
        # the x-axis labels by 45 degrees and shift to the right
        axis.text.y = element_text(size=13), # size of y axis ticks
        axis.title = element_text(size = 15), # size of axis names
        plot.title = element_text(size = 25), # size of plot title
        plot.caption = element_text(size=23), # size of plot caption
        legend.text = element_markdown(size = 13), # size of legend text
        legend.title = element_text(size = 17), # size of legend title
        legend.position = "right") # legend on the right
#+ fig.height=8, fig.width=11
print(mainplot+
        ggtitle(paste(agglom.rank,"level gut microbiota profiles of fecal samples from naked mole-rats"))+
        theme(plot.title = element_text(size=15))
)

# for(image.format in image.formats){
#   ggsave(paste0(barplot.directory,
#                 paste(paste(format(Sys.time(),format="%Y%m%d"),
#                             format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                       "kraken2","barplot","NMR",
#                       agglom.rank,
#                       sep = "-"),".",image.format),
#          plot=mainplot,
#          width=11, height=7,units="in",
#          dpi=300,device = image.format)
# }

#+ echo=FALSE
## 4. Relative abundance per kingdom (including bacteria). ####
#'
#' ## Relative abundance per kingdom (including bacteria).
#' Read table with unclassified read numbers.
unclassified_reads<-read.table("./output/kraken2_pipeline/20240409_17_32_40_unclassified_reads.tsv",
                               header = T)
colnames(unclassified_reads)[which(colnames(unclassified_reads)=="Unclassified")]<-"Abundance"
unclassified_reads$Kingdom<-"Unclassified"

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
  ylab("Relative abundance from all samples (%)")+
  xlab("")+
  theme_bw()+
  coord_cartesian(expand = c("bottom" = F,
                             "left" = F,
                             "right" = F))+
  geom_text(aes(x=Kingdom, 
                y=TotalAbundancePerTax, 
                label=round(TotalAbundancePerTax,digits = 3)),
            vjust=-0.5,size=5)+ # add text with percentage above the bar
  theme(axis.line = element_blank(), 
        axis.text.x = element_text(size=13),
        axis.text.y = element_text(size=15), # size of y axis ticks
        axis.title = element_text(size = 15), # size of axis names
        plot.title = ggtext::element_markdown(size = 25), # the plot 
        # title will be recognised as a markdown object, so we can
        # add line breaks (cause host names are too long)
        plot.caption = element_text(size=23),# size of plot caption
        legend.text = element_text(size = 13),# size of legend text
        legend.title = element_text(size = 17), # size of legend title
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) # no legend 
#+ fig.height=8, fig.width=11
print(kingdom.plot.with_bact+
        ggtitle("Whole metagenome sequencing profile of naked mole-rat gut microbiota
               (Kingdom level)")+
        theme(plot.title = ggtext::element_markdown(size=15))
)
# for(image.format in image.formats){
#   ggsave(paste0(barplot.directory,
#                 paste(paste(format(Sys.time(),format="%Y%m%d"),
#                             format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                       "NMR","kingdom.plot.with_unclas",
#                       sep = "-"),".",image.format),
#          plot=kingdom.plot.with_bact,
#          width=11, height=6,units="in",
#          dpi=300,device = image.format)
# }

#' Check library size per sample (not included)
ps.q.agg%>%
  group_by(Sample,Kingdom)%>%
  summarise(Abundance=sum(Abundance))%>%
  bind_rows(unclassified_reads)%>%
  group_by(Sample)%>%
  summarise(TotalSample=sum(Abundance))%>%
  ggplot(aes(x=Sample,y=TotalSample,fill=TotalSample))+
  geom_bar(stat = "identity")

sessionInfo()
rm(list = ls(all=TRUE))
gc()