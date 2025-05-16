# Creating barplots ####

# After importing the kraken2 files into R and processing them with phyloseq,
# it's time to explore the taxonomic composition of our data.
# We will use the Polychrome package to create a custom palette for the 
# barplots.
library(phyloseq)
library(tidyverse)
library(Polychrome)
library(ggtext)
library(ggVennDiagram)
venn.directory<-"images/venn/"
## Specifying parameters and directory/file names for shotgun data #### 
agglom.rank.wms<-"Phylum"# this is the taxonomic rank that was used for agglomeration
barplot.directory<-"./images/barplots/" # set the path where barplots will
# be saved
## Specifying parameters and directory/file names for amplicon data #### 
authorname<-"pooled" # name of the folder with QIIME2 output
agglom.rank<-agglom.rank.wms # this is the taxonomic rank that was used for agglomeration
truncationlvl<-"234" #  truncation level that we chose in QIIME2
read.end.type<-"single" # single reads or paired reads: decided in QIIME2

# Load the Workspace from phyloseq (output of 001-phyloseq-qiime2.R)
load(paste0("../amplicon_nmr/output/rdafiles/",paste(authorname,read.end.type,"qiime2",
                                       truncationlvl,agglom.rank,
                                       "phyloseq-workspace.RData",sep = "-")))
ps.q.agg.amp<-ps.q.agg
agglom.rank.amp<-agglom.rank
objects.to.keep<-c("agglom.rank.amp","ps.q.agg.amp",
                         "agglom.rank.wms","barplot.directory",
                   "venn.directory")
objects.to.keep<-which(ls()%in%objects.to.keep)
rm(list = ls()[-objects.to.keep])

load(paste0("./output/rdafiles/",paste("kraken2",
                                  agglom.rank.wms,
                                  "phyloseq-workspace.RData",sep = "-")))
ps.q.agg.wms<-ps.q.agg
rm(ps.q.agg)
rm(agglom.rank)

# pretty.facet.labels<-c("NMR" = "*Heterocephalus glaber*")  
# custom.levels<-intersect(names(pretty.facet.labels),custom.md$class)
pretty.facet.labels<-c("NMR_amplicon" = "*Heterocephalus glaber* 16S data",
                       "NMR_wms"= "*Heterocephalus glaber* WMS data")
custom.levels<-names(pretty.facet.labels)
ps.q.agg.amp<-ps.q.agg.amp%>%
  filter(Sample%in%ps.q.agg.wms$Sample)%>% # amplicon data has more samples
  # filter(class%in%custom.levels,Abundance!=0)%>%
  filter(Abundance!=0)%>%
  mutate(class="NMR_amplicon")
ps.q.agg.wms<-ps.q.agg.wms%>%
  # filter(class%in%custom.levels,Abundance!=0)%>%
  filter(Abundance!=0)%>%
  mutate(class="NMR_wms")

# bind two dataframes
ps.q.agg<-rbind(ps.q.agg.amp,ps.q.agg.wms)
ps.q.agg<-ps.q.agg%>%
  filter(Kingdom=="Bacteria")%>%
  group_by(class,Sample)%>%
  mutate(TotalSample=sum(Abundance))%>%
  group_by_at(c("class","Sample",agglom.rank.wms))%>%
  mutate(RelativeAbundance=Abundance/TotalSample*100)%>%
  group_by(class)%>%
  mutate(TotalClass=sum(Abundance))%>%
  mutate(MeanRelativeAbundance = Abundance/TotalClass*100)%>%
  select(-TotalSample,-TotalClass)

# Find unique and shared taxa
amp.set<-ps.q.agg%>%
  group_by(class)%>%
  filter(class=="NMR_amplicon")%>%
  # filter(MeanRelativeAbundance>=0.1)%>%
  select(all_of(agglom.rank.amp))%>%
  unique()%>%
  pull()

wms.set<-ps.q.agg%>%
  group_by(class)%>%
  filter(class=="NMR_wms")%>%
  # filter(MeanRelativeAbundance>=0.1)%>%
  select(all_of(agglom.rank.wms))%>%
  unique()%>%
  pull()
wms.uniq<-setdiff(wms.set,amp.set)
amp.uniq<-setdiff(amp.set,wms.set)

taxa.sep<-list(amp=amp.set,wms=wms.set)
# Quick venn
venn.taxa.all<-
  ggVennDiagram(taxa.sep,
                # category.names = c("Phyla from \namplicon sequencing",
                #                  "Phyla from \nwhole metagenome sequencing"),
                # category.names = c(),
                set_size = 10, # size of set labels
                label = "count",# what kind of numbers are shown
                label_size = 10 # size of numbers
                )+  
  scale_fill_gradient(low="grey90",high = "red")+
  # expand axis to show long set labels
  scale_x_continuous(expand = expansion(mult = .3))+
  theme(axis.line = element_blank(), 
        plot.title = ggtext::element_markdown(size = 25), # the plot 
        # title will be recognised as a markdown object, so we can
        # add line breaks (cause host names are too long)
        plot.caption = element_text(size=23),# size of plot caption
        legend.text = element_text(size = 20),# size of legend text
        legend.title = element_text(size = 25), # size of legend title
        legend.position = "none") # legend on the right
ggsave(paste0(venn.directory,paste(Sys.Date(),"venn.all",
                                   agglom.rank.wms,sep = "-"),
              ".png"),
       plot = venn.taxa.all,
       width = 4000,height = 6000,
       units = "px",dpi=300,device = "png")

# >=0.1 %
amp.set.0.1pc<-ps.q.agg%>%
  group_by_at(c("class",agglom.rank.amp))%>%
  filter(class=="NMR_amplicon")%>%
  filter(MeanRelativeAbundance>=0.1)%>%
  ungroup()%>%
  select(all_of(agglom.rank.amp))%>%
  unique()%>%
  pull()

wms.set.0.1pc<-ps.q.agg%>%
  group_by_at(c("class",agglom.rank.wms))%>%
  filter(class=="NMR_wms")%>%
  filter(MeanRelativeAbundance>=0.1)%>%
  select(all_of(agglom.rank.wms))%>%
  unique()%>%
  pull()
taxa.sep.0.1pc<-list(amp.set.0.1pc=amp.set.0.1pc,
               wms.set.0.1pc=wms.set.0.1pc)

ggVennDiagram(taxa.sep.0.1pc,
              # category.names = c("Phyla from \namplicon sequencing",
              #                  "Phyla from \nwhole metagenome sequencing"),
              # category.names = c(),
              set_size = 8, # size of set labels
              label = "count",# what kind of numbers are shown
              label_size = 8 # size of numbers
              )+  
  scale_fill_gradient(low="grey90",high = "red")+
  # expand axis to show long set labels
  scale_x_continuous(expand = expansion(mult = .3))+
  theme(axis.line = element_blank(), 
        plot.title = ggtext::element_markdown(size = 25), # the plot 
        # title will be recognised as a markdown object, so we can
        # add line breaks (cause host names are too long)
        plot.caption = element_text(size=23),# size of plot caption
        legend.text = element_text(size = 20),# size of legend text
        legend.title = element_text(size = 25), # size of legend title
        legend.position = "none") # legend on the rightamp.set.0.1pc
amp.set.0.1pc
wms.set.0.1pc



amp<-subset(ps.q.agg[,c("class","Sample","Phylum","RelativeAbundance")],
            class=="NMR_amplicon")
wms<-subset(ps.q.agg[,c("class","Sample","Phylum","RelativeAbundance")],
            class=="NMR_wms")

inner_join(amp,wms,by=c("Phylum","Sample"))%>%
  mutate(diff=abs(RelativeAbundance.x-RelativeAbundance.y))%>%
  group_by(Phylum)%>%
  summarise(means=mean(diff))

# Barplot
ps.q.agg%>%
  # group_by(class,Sample)%>%
  ungroup()%>%
  mutate(Taxon=ifelse(MeanRelativeAbundance>=0.1,Phylum,"Remainder"))%>%
  ggplot(aes(x=Sample,y=RelativeAbundance,fill=Taxon))+
  geom_bar(stat = "identity")+
  facet_grid(~class)
