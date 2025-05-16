## Import libraries ####
# install.packages(c("tidyverse","Polychrome","ggtext"))
library(tidyverse)
agglom.rank<-"Species"
pipeline.name<-"kraken2"
rdafiles.directory<-"./output/rdafiles"
rtables.directory<-"./output/rtables"
source("../amplicon_nmr/code/r-scripts/get_dominant_taxa_in_host.R")
source("../amplicon_nmr/code/r-scripts/create_summary_stats_table.R")
source("../amplicon_nmr/code/r-scripts/get_n_uniq_taxa_per_host.R")
source("../amplicon_nmr/code/r-scripts/add_relab_to_tax_df.R")
source("../amplicon_nmr/code/r-scripts/add_agegroup_to_tax_df.R")
source("../amplicon_nmr/code/r-scripts/get_unclassified_summary_stats.R")
source("../amplicon_nmr/code/r-scripts/ggplot_species.R")
source("../amplicon_nmr/code/r-scripts/add_zero_rows.R")

# Import datasets as rds files
if(pipeline.name=="kraken2"){
  ps.q.agg.date_time<-"20241003_13_52_43"
  ps.q.agg.genus.date_time<-"20241003_13_56_28"
  ps.q.agg.family.date_time<-"20241003_13_57_30"
  ps.q.agg.phylum.date_time<-"20241003_13_57_59"
}else if (pipeline.name=="singlem"){
  ps.q.agg.date_time<-"20240929_23_33_37"
  ps.q.agg.genus.date_time<-"20240929_23_33_59"
  ps.q.agg.family.date_time<-"20240929_23_34_25"
  ps.q.agg.phylum.date_time<-"20240929_23_34_50"
  ps.q.agg.kraken.date_time<-"20241003_13_52_43"
}

ps.q.agg<-readRDS(file=file.path(
  rdafiles.directory,
  paste(ps.q.agg.date_time,"phyloseq",pipeline.name,"Species",
        "table.rds",sep = "-")))
ps.q.agg.genus<-readRDS(file=file.path(
  rdafiles.directory,
  paste(ps.q.agg.genus.date_time,"phyloseq",pipeline.name,"Genus",
        "table.rds",sep = "-")))

ps.q.agg.family<-readRDS(file=file.path(
  rdafiles.directory,
  paste(ps.q.agg.family.date_time,"phyloseq",pipeline.name,"Family",
        "table.rds",sep = "-")))
ps.q.agg.phylum<-readRDS(file=file.path(
  rdafiles.directory,
  paste(ps.q.agg.phylum.date_time,"phyloseq",pipeline.name,"Phylum",
        "table.rds",sep = "-")))
# Import metadata
custom.md<-readRDS("../amplicon_nmr/output/rdafiles/use/custom.md.ages.rds")
custom.md<-custom.md%>%
  filter(Sample%in%ps.q.agg$Species)

sample.levels<-custom.md%>%
  ungroup()%>%
  select(Sample,age)%>%
  arrange(age)%>%
  distinct()%>%
  mutate(NewSample=paste0(Sample," (",age,")"))%>%
  mutate(Sample=factor(Sample,levels=Sample),
         NewSample=factor(NewSample,levels=NewSample))


## Total ASV/phyla/families/genera per class ####
n.species.per.host<-get_n_uniq_taxa_per_host(ps.q.agg,"Species")
n.phylum.per.host<-get_n_uniq_taxa_per_host(ps.q.agg.phylum,"Phylum")
n.family.per.host<-get_n_uniq_taxa_per_host(ps.q.agg.family,"Family")
n.genus.per.host<-get_n_uniq_taxa_per_host(ps.q.agg.genus,"Genus")

## Summary table ####
# Number of samples per host, total reads per host, mean library size,
# sd library size, asv per host, phyla per host, families per host, 
# genera per host
summary.table<-create_summary_stats_table(ps.q.agg,
                                          n.asv.table = n.species.per.host,
                                          n.phylum.table = n.phylum.per.host,
                                          n.family.table = n.family.per.host,
                                          n.genus.table = n.genus.per.host)
# write.table(summary.table,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  pipeline.name,"summary-table.tsv",sep="-")),
#             row.names = F,sep = "\t")

## Add relative abundance, mean relative abundance, total abundance, etc ####
ps.q.agg.phylum.relab<-add_relab_to_tax_df(ps.q.agg.phylum,"Phylum")
ps.q.agg.family.relab<-add_relab_to_tax_df(ps.q.agg.family,"Family")
ps.q.agg.genus.relab<-add_relab_to_tax_df(ps.q.agg.genus,"Genus")
ps.q.agg.relab<-add_relab_to_tax_df(ps.q.agg,"Species")

## Add agegroup (Must run for plotting) ####
# we create these new levels for plots
ps.q.agg.relab<-add_agegroup_to_tax_df(ps.q.agg.relab,"Species",custom.md)
ps.q.agg.genus.relab<-add_agegroup_to_tax_df(ps.q.agg.genus.relab,"Genus",custom.md)
ps.q.agg.family.relab<-add_agegroup_to_tax_df(ps.q.agg.family.relab,"Family",custom.md)

if(pipeline.name=="singlem"){
  ps.q.agg.kraken<-readRDS(file.path(rdafiles.directory,paste(ps.q.agg.kraken.date_time,"phyloseq-kraken2-Species-table.rds",sep="-")))
  ps.q.agg.kraken.relab<-add_relab_to_tax_df(ps.q.agg.kraken,"Species")
}
### Check how many taxa are unclassified in each NMR sample ####
all.ranks<-c("Kingdom", "Phylum", "Class", "Order", "Family","Genus","Species")
agglom.rank<-"Species"
unclassified.species.summary.stats.table<-get_unclassified_summary_stats(ps.q.agg.relab,
                                                                         "Species")

# write.table(unclassified.species.summary.stats.table,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  pipeline.name,"unclassified-species-summary-table.tsv",sep="-")),
#             row.names = F,sep = "\t")

if(pipeline.name=="singlem"){
  unclassified.species.summary.stats.table.kraken<-get_unclassified_summary_stats(ps.q.agg.kraken.relab,"Species")
  # Barplot that shows how many unique species were found in Kraken2 and SingleM, 
  # and how many species were unclassified in SingleM.
  unclassified.species.summary.stats.table%>%
    rename("ClassifiedSpecies"=NumClassifiedSpecies,
           "UnclassifiedSpecies"=NumUnclassifiedSpecies)%>%
    left_join(unclassified.species.summary.stats.table.kraken)%>%
    rename("Kraken2ClassifiedSpecies"=NumClassifiedSpecies)%>%
    pivot_longer(c("ClassifiedSpecies","UnclassifiedSpecies","Kraken2ClassifiedSpecies"),
                 names_to = "SpeciesType",
                 values_to = "NumberSpecies")%>%
    mutate(Pipeline=ifelse(SpeciesType=="Kraken2ClassifiedSpecies","Kraken2","SingleM"))%>%
    relocate(SpeciesType,NumberSpecies,Pipeline)%>%
    left_join(sample.levels)%>%
    ggplot(aes(x=NewSample,y=NumberSpecies,fill=SpeciesType))+
    geom_bar(stat = "identity")+
    facet_grid(~Pipeline)+
    coord_cartesian(expand=F)+
    labs(x="",
         y="Number of species")+
    theme_bw()+
    theme(plot.margin=unit(c(1,1,1,1.5), 'cm'),
          axis.line = element_blank(), #TODO: what does it do?
          strip.text.x = ggtext::element_markdown(size = 20),# the name of
          # each facet will be recognised as a markdown object, so we can
          # add line breaks (cause host names are too long)
          panel.spacing = unit(0.8, "cm"), # increase distance between facets
          axis.text.x = element_text(angle=45,size=20,hjust=1),# rotate
          # the x-axis labels by 45 degrees and shift to the right
          axis.text.y = element_text(size=20), # size of y axis ticks
          axis.title = element_text(size = 20), # size of axis names
          plot.title = element_text(size = 25), # size of plot title
          plot.caption = element_text(size=23), # size of plot caption
          legend.text = element_markdown(size = 20), # size of legend text
          legend.title = element_text(size = 25), # size of legend title
          legend.position = "right") # legend under the plot
  
  ggsave(filename=paste0("./images/barplots/",
                         paste(format(Sys.time(),format="%Y%m%d"),
                               format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                         "-singlem-vs-kraken2.png"),
         plot = last_plot(),device = "png",width = 4800,height = 3000,
         units = "px",
         dpi=300)
  
}




# Distinct species in kraken2 vs singlem
ps.q.agg%>%
  group_by(class)%>%
  filter(Abundance!=0)%>%
  distinct(Species)%>%
  summarise(Taxon_count=n())%>%
  arrange(-Taxon_count)
if(pipeline.name=="singlem"){
  ps.q.agg.kraken%>%
  group_by(class)%>%
  filter(Abundance!=0)%>%
  distinct(Species)%>%
  summarise(Taxon_count=n())%>%
  arrange(-Taxon_count)
}
### Check the number of distinct OTU/taxa ####
ps.q.agg%>%
  group_by(class)%>%
  filter(Abundance!=0)%>%
  distinct(get(agglom.rank))%>%
  summarise(Taxon_count=n())%>%
  arrange(-Taxon_count)
ps.q.agg.genus%>%
  group_by(class)%>%
  filter(Abundance!=0)%>%
  distinct(Genus)%>%
  summarise(Taxon_count=n())%>%
  arrange(-Taxon_count)

### Most abundant phyla ####
ps.q.agg.dominant.phyla<-get_dominant_taxa_in_host(tax.df=ps.q.agg.phylum,
                                           tax.rank="Phylum",
                                           host="NMR",
                                           nonbacterial.table = F)
# write.table(ps.q.agg.dominant.phyla,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  pipeline.name,"dominant-phyla.tsv",sep="-")),
#             row.names = F,sep = "\t")

### Most abundant non-bacterial phyla  ####
ps.q.agg.dominant.phyla.nonbact<-get_dominant_taxa_in_host(tax.df=ps.q.agg.phylum,
                                                   tax.rank="Phylum",
                                                   host="NMR",
                                                   nonbacterial.table = T)
# write.table(ps.q.agg.dominant.phyla.nonbact,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  pipeline.name,"dominant-phyla-nonbact.tsv",sep="-")),
#             row.names = F,sep = "\t")

### Most abundant families ####
ps.q.agg.dominant.families<-get_dominant_taxa_in_host(tax.df=ps.q.agg.family,
                                                      tax.rank ="Family",
                                                      host="NMR",
                                                      nonbacterial.table = F)
# write.table(ps.q.agg.dominant.families,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  pipeline.name,"dominant-families.tsv",sep="-")),
#             row.names = F,sep = "\t")

### Most abundant non-bacterial families ####
ps.q.agg.dominant.families.nonbact<-get_dominant_taxa_in_host(tax.df=ps.q.agg.family,
                                                              tax.rank="Family",
                                                              host="NMR",
                                                              nonbacterial.table = T)
# write.table(ps.q.agg.dominant.families.nonbact,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  pipeline.name,"dominant-families-nonbact.tsv",sep="-")),
#             row.names = F,sep = "\t")
### Most abundant genera ####
ps.q.agg.dominant.genera<-get_dominant_taxa_in_host(tax.df=ps.q.agg.genus,
                                                    tax.rank="Genus",
                                                    host="NMR",
                                                    nonbacterial.table = F)
# write.table(ps.q.agg.dominant.genera,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  pipeline.name,"dominant-genera.tsv",sep="-")),
#             row.names = F,sep = "\t")

### Most abundant non-bacterial genera ####
ps.q.agg.dominant.genera.nonbact<-get_dominant_taxa_in_host(tax.df=ps.q.agg.genus,
                                                            tax.rank="Genus",
                                                            host="NMR",
                                                            nonbacterial.table = T)
# write.table(ps.q.agg.dominant.genera.nonbact,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  pipeline.name,"dominant-genera-nonbact.tsv",sep="-")),
#             row.names = F,sep = "\t")

### Most abundant species overall ####
ps.q.agg.dominant.species<-get_dominant_taxa_in_host(tax.df=ps.q.agg,
                                                     tax.rank="Species",
                                                     host="NMR",
                                             nonbacterial.table = F)
# write.table(ps.q.agg.dominant.species,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  pipeline.name,"dominant-species-all.tsv",sep="-")),
#             row.names = F,sep = "\t")
### Most abundant non-bacterial species ####
ps.q.agg.dominant.species.nonbact<-get_dominant_taxa_in_host(tax.df=ps.q.agg,
                                                             tax.rank="Species",
                                                             host="NMR",
                                                             nonbacterial.table = T)
# write.table(ps.q.agg.dominant.species.nonbact,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  pipeline.name,"dominant-species-nonbact.tsv",sep="-")),
#             row.names = F,sep = "\t")




# Compare with Debebe et al. ####
### How much Bacteroidaceae are in NMR  ####
bacteroidaceae.nmr<-ps.q.agg.family.relab%>%
  filter(Family=="Bacteroidaceae",class=="NMR")%>%
  mutate(min=min(RelativeAbundance),
         max=max(RelativeAbundance),
         mean=TotalAgglomRank/TotalClass*100,
         sd=sd(RelativeAbundance),
         n=n())%>%
  select(Family,min,max,mean,sd,n)%>%
  distinct()
# write.table(bacteroidaceae.nmr,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  pipeline.name,"bacteroidaceae-nmr-table.tsv",sep="-")),
#             row.names = F,sep = "\t")

### What are the most dominant Bacteroidota families in NMR ####
bacteroidota.nmr<-ps.q.agg.family.relab%>%
  filter(Phylum=="Bacteroidota",class=="NMR")%>%
  group_by(Family)%>%
  distinct(Family,.keep_all = T)%>%
  arrange(-MeanRelativeAbundance)%>%
  select(Family,MeanRelativeAbundance)
# write.table(bacteroidota.nmr,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  pipeline.name,"bacteroidota-nmr-table.tsv",sep="-")),
#             row.names = F,sep = "\t")

### Spirochaetaceae and Treponema ####
spirochaetaceae.nmr<-ps.q.agg.family.relab%>%
  filter(Family=="Spirochaetaceae",class=="NMR")%>%
  mutate(min=min(RelativeAbundance),
            max=max(RelativeAbundance),
            mean=TotalAgglomRank/TotalClass*100,
            sd=sd(RelativeAbundance),
            n=n())%>%
  select(Family,min,max,mean,sd,n)%>%
  distinct()
spirochaetota.nmr<-ps.q.agg.phylum.relab%>%
  filter(Phylum=="Spirochaetota",class=="NMR")%>%
  mutate(min=min(RelativeAbundance),
         max=max(RelativeAbundance),
         mean=TotalAgglomRank/TotalClass*100,
         sd=sd(RelativeAbundance),
         n=n())%>%
  select(Phylum,min,max,mean,sd,n)%>%
  distinct()
# write.table(spirochaetaceae.nmr,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  pipeline.name, "spirochaetaceae-nmr-table.tsv",sep="-")),
#             row.names = F,sep = "\t")
# write.table(spirochaetota.nmr,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  pipeline.name, "spirochaetota-nmr-table.tsv",sep="-")),
#             row.names = F,sep = "\t")

treponema.nmr<-ps.q.agg.genus.relab%>%
  filter(Genus=="Treponema",class=="NMR")%>%
  mutate(min=min(RelativeAbundance),
         max=max(RelativeAbundance),
         mean=TotalAgglomRank/TotalClass*100,
         sd=sd(RelativeAbundance),
         n=n())%>%
  select(Genus,min,max,mean,sd,n)%>%
  distinct()
treponema.nmr.names<-ps.q.agg.relab%>%
  filter(Genus=="Treponema",class=="NMR")%>%
  select(Species,MeanRelativeAbundance)%>%
  distinct()%>%
  arrange(-MeanRelativeAbundance)
# write.table(treponema.nmr,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  pipeline.name,"treponema-nmr-table.tsv",sep="-")),
#             row.names = F,sep = "\t")
# write.table(treponema.nmr.names,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  pipeline.name,"treponema-nmr-names.tsv",sep="-")),
#             row.names = F,sep = "\t")

#### How many ASVs in Treponema ####
ps.q.agg%>%
  ungroup()%>%
  filter(Genus=="Treponema",class=="NMR")%>%
  distinct(OTU)%>%
  tally

### Mogibacteriaceae is renamed to Anaerovoracaceae ####
mogibacteriaceae_anaerovoracaceae.all<-ps.q.agg.family.relab%>%
  filter(Family=="Anaerovoracaceae")%>%
  mutate(min=min(RelativeAbundance),
         max=max(RelativeAbundance),
         mean=TotalAgglomRank/TotalClass*100,
         sd=sd(RelativeAbundance),
         n=n())%>%
  select(Family,min,max,mean,sd,n)%>%
  distinct()%>%
  arrange(-mean)
# write.table(mogibacteriaceae_anaerovoracaceae.all,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                 pipeline.name, "mogibacteriaceae_anaerovoracaceae-all-table.tsv",sep="-")),
#             row.names = F,sep = "\t")



## Setup plots ####
library(Polychrome)
library(ggtext)

pretty.level.names<-names(table(custom.md$old_agegroup))
names(pretty.level.names)<-names(table(custom.md$agegroup))
custom.levels<-names(pretty.level.names)
gg.labs.name<-"Age group"
gg.title.groups<-"age groups"


set.seed(1)
custom.fill<-createPalette(length(custom.levels),
                           seedcolors = c("#EE2C2C","#5CACEE","#00CD66",
                                          "#FF8C00","#BF3EFF", "#00FFFF",
                                          "#FF6EB4","#00EE00","#EEC900",
                                          "#FFA07A"))
names(custom.fill)<-custom.levels
swatch(custom.fill)
all.ranks<-c("Kingdom", "Phylum", "Class", "Order", "Family","Genus","Species")

## Sulfur metabolising bacteria ####
# We don't know where the "sulf" pattern is (Phylum, Order, Genus, or anywhere else),
# so we use sapply with grepl
# https://stackoverflow.com/questions/47941680/grepl-across-multiple-specified-columns
sulfur.bact.nmr<-
  ps.q.agg.relab[!!rowSums(sapply(ps.q.agg.relab[,all.ranks], grepl, pattern = "sulf|thio") ),]%>%
  group_by(Species)%>%
  mutate(min=min(RelativeAbundance),
         max=max(RelativeAbundance),
         mean=TotalAgglomRank/TotalClass*100,
         sd=sd(RelativeAbundance),
         n=n())%>%
  select(Kingdom:Species,MeanRelativeAbundance,sd)%>%
  distinct()%>%
  arrange(-MeanRelativeAbundance)

# write.table(sulfur.bact.nmr,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  pipeline.name,"sulfur.bact.nmr-nmr-table.tsv",sep="-")),
#             row.names = F,sep = "\t")

### Plot sulfur metabolising bacteria ####
# Add zero rows
# ps.q.agg.relab<-ps.q.agg.relab%>%ungroup()
ps.q.agg.relab<-add_zero_rows(unique(sulfur.bact.nmr$Species),ps.q.agg.relab,"Species")
ps.q.agg.relab<-ps.q.agg.relab%>%
  group_by(Sample)%>%
  fill(agegroup,.direction = "down")
ggplot_species(taxa.to.plot=sulfur.bact.nmr$Species,
               tax.df=ps.q.agg.relab,
               tax.rank="Species",
               sample.order=sample.levels,
               group.names = pretty.level.names,
               grouping.variable = "agegroup",
               metadata.df = custom.md,
               ggplot.fill.name = gg.labs.name,
               ggplot.fill.vector = custom.fill)+
  ggtitle(paste0("Relative abundance of sulfur-utilizing bacteria in different naked mole-rat age groups"))


for (image_format in c("png","tiff")){
  ggsave(paste0("./images/barplots/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "sulfur-bacteria-nmr",
                      sep = "-"),".",image_format),
         plot=last_plot(),
         width = 7000,height = 9000,
         units = "px",dpi=300,device = image_format)
}


## Plot Treponema ####
treponema.sp.nmr<-ps.q.agg.relab%>%
  filter(Genus=="Treponema",class=="NMR")%>%
  group_by(Species)%>%
  select(Species,MeanRelativeAbundance)%>%
  arrange(-MeanRelativeAbundance)%>%
  ungroup()%>%
  distinct()%>%
  pull(Species)
ps.q.agg.relab<-ps.q.agg.relab%>%ungroup()
ps.q.agg.relab<-add_zero_rows(treponema.sp.nmr,ps.q.agg.relab,"Species")
ps.q.agg.relab<-ps.q.agg.relab%>%
  group_by(Sample)%>%
  fill(agegroup,.direction = "down")

ggplot_species(taxa.to.plot=treponema.sp.nmr,
               tax.df=ps.q.agg.relab,
               tax.rank="Species",
               sample.order=sample.levels,
               group.names = pretty.level.names,
               grouping.variable = "agegroup",
               metadata.df = custom.md,
               ggplot.fill.name = gg.labs.name,
               ggplot.fill.vector = custom.fill)+
  ggtitle(paste0("Relative abundance of Treponema genus members"))



for (image_format in c("png","tiff")){
  ggsave(paste0("./images/barplots/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "treponema.sp-nmr",
                      sep = "-"),".",image_format),
         plot=last_plot(),
         width = 7000,height = 9000,
         units = "px",dpi=300,device = image_format)
}


## Compare with other papers ####
### Debebe culturing paper ####
debebe.species<-c("Bacteroides_thetaiotaomicron",
                   "Bacteroides_ovatus",
                   "Mammaliicoccus_sciuri")
ps.q.agg.relab<-ps.q.agg.relab%>%ungroup()
ps.q.agg.relab<-add_zero_rows(debebe.species,ps.q.agg.relab,"Species")
ps.q.agg.relab<-ps.q.agg.relab%>%
  group_by(Sample)%>%
  fill(agegroup,.direction = "down")

ps.q.agg.relab%>%
  filter(Species%in%debebe.species&!is.na(MeanRelativeAbundance))%>%
  ungroup()%>%
  select(Species,MeanRelativeAbundance)%>%
  distinct()


ggplot_species(taxa.to.plot=debebe.species,
               tax.df=ps.q.agg.relab,
               tax.rank="Species",
               sample.order=sample.levels,
               group.names = pretty.level.names,
               grouping.variable = "agegroup",
               metadata.df = custom.md,
               ggplot.fill.name = gg.labs.name,
               ggplot.fill.vector = custom.fill)+
  ggtitle(paste0("Relative abundance of species from Debebe et al. (2016)"))


### Debebe et al. (16S) ####
debebe16s.families<-c("Anaerovoracaceae",
  "Desulfovibrionaceae",
  "Dethiosulfovibrionaceae", 
  "Desulfobulbaceae",
  "Desulfurococcaceae",
  "Haloarculaceae")

ps.q.agg.family.relab<-ps.q.agg.family.relab%>%ungroup()
ps.q.agg.family.relab<-add_zero_rows(debebe16s.families,ps.q.agg.family.relab,"Family")
ps.q.agg.relab<-ps.q.agg.relab%>%
  group_by(Sample)%>%
  fill(agegroup,.direction = "down")

ggplot_species(taxa.to.plot=debebe16s.families,
               tax.df=ps.q.agg.family.relab,
               tax.rank="Family",
               sample.order=sample.levels,
               group.names = pretty.level.names,
               grouping.variable = "agegroup",
               metadata.df = custom.md,
               ggplot.fill.name = gg.labs.name,
               ggplot.fill.vector = custom.fill)+
  ggtitle(paste0("Relative abundance of species from Debebe et al. (16S)"))

ps.q.agg.family.relab%>%
  filter(Family%in%debebe16s.families&!is.na(MeanRelativeAbundance))%>%
  ungroup()%>%
  select(Family,MeanRelativeAbundance)%>%
  distinct()

### From my data ####
species.to.plot<-c("Bacteroides_xylanisolvens",
  "Enterococcus_faecium",
  "Anaerostipes_caccae",
  "Phascolarctobacterium_faecium",
  "Bacteroides_thetaiotaomicron",
  "Parabacteroides_distasonis",
  "Clostridioides_difficile",
  "Ipomoea_triloba",
  "Cryptomeria_japonica",
  "Macadamia_integrifolia")

ps.q.agg.relab<-ps.q.agg.relab%>%ungroup()
ps.q.agg.relab<-add_zero_rows(species.to.plot,ps.q.agg.relab,"Species")
ps.q.agg.relab<-ps.q.agg.relab%>%
  group_by(Sample)%>%
  fill(agegroup,.direction = "down")

my.species.plot<-ggplot_species(taxa.to.plot=species.to.plot,
               tax.df=ps.q.agg.relab,
               tax.rank="Species",
               sample.order=sample.levels,
               group.names = pretty.level.names,
               grouping.variable = "agegroup",
               metadata.df = custom.md,
               ggplot.fill.name = gg.labs.name,
               ggplot.fill.vector = custom.fill)+
  ggtitle(paste0("Relative abundance of species from my data"))

for (image_format in c("png","tiff")){
  ggsave(paste0("./images/barplots/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "myspecies.plot-nmr",
                      sep = "-"),".",image_format),
         plot=my.species.plot,
         width = 7000,height = 5000,
         units = "px",dpi=300,device = image_format) 
}


### Groups 1, 2, and 3 ####
group1.genera<-c("Faecalibacterium",
                 "Roseburia",
                 "Coprococcus",
                 "Prevotella")
group1.species<-c("Eubacterium_rectale") # Not in our data
group2.increased.genera<-c("Eggerthella",
                           "Bilophila",
                           "Desulfovibrio",
                           "Fusobacterium",
                           "Anaerotruncus",
                           "Streptococcus",
                           "Escherichia")
group2.unhealthy.genera<-c("Eggerthella",
                           "Coprobacillus",
                           "Streptococcus",
                           "Bilophila",
                           "Actinomyces",
                           "Desulfovibrio",
                           "Campylobacter",
                           "Veillonella",
                           "Enterococcus")
group2.unhealthy.families<-c("Atopobiaceae",
                             "Enterobacteriaceae")
group2.unhealthy.species<-c("Bacteroides_fragilis",
                            "Clostridium_hathewayi",
                            "Enterocloster_clostridioformis", #Clostridium_clostridioforme
                            "Enterocloster_bolteae", # Clostridium_bolteae
                            "Clostridium_cindens",
                            "[Ruminococcus]_torques", #Ruminococcus_torques
                            "Mediterraneibacter_gnavus", # Ruminococcus_gnavus
                            "Clostridioides_difficile")

group3.genera<-c("Akkermansia",
                 "Odoribacter",
                 "Butyricimonas",
                 "Butyrivibrio",
                 "Barnesiella",
                 "Oscillospira")
group3.families<-c("Christensenellaceae")

group1.genera%in%ps.q.agg.genus.relab$Genus
group1.species%in%ps.q.agg.relab$Species
group2.increased.genera%in%ps.q.agg.genus.relab$Genus
group2.unhealthy.genera%in%ps.q.agg.genus.relab$Genus
group2.unhealthy.families%in%ps.q.agg.family.relab$Family
group2.unhealthy.species%in%ps.q.agg.relab$Species
group3.genera%in%ps.q.agg.genus.relab$Genus
group3.families%in%ps.q.agg.family.relab$Family

# Mean relative abundance of group 1 genera by age group
show.mean.abundances.and.plot<-function(taxa.to.plot,
                                        tax.df,
                                        tax.rank,
                                        sample.order,
                                        group.names,
                                        grouping.variable,
                                        metadata.df,
                                        ggplot.fill.name,
                                        ggplot.fill.vector){
  tax.df%>%
    filter(get(tax.rank)%in%taxa.to.plot&!is.na(MeanRelativeAbundance))%>%
    group_by_at(c(tax.rank,"agegroup"))%>%
    select(all_of(c(tax.rank,"MeanRelativeAbundance","MeanRelativeAbundanceAgegroup")))%>%
    ungroup()%>%
    distinct()%>%
    print()
  tax.df%>%
    filter(get(tax.rank)%in%taxa.to.plot&!is.na(MeanRelativeAbundance))%>%
    ungroup()%>%
    select(all_of(c(tax.rank,"MeanRelativeAbundance")))%>%
    distinct()%>%
    arrange(-MeanRelativeAbundance)%>%
    pull(tax.rank) -> taxa.to.plot
  tax.df<-add_zero_rows(taxa.to.plot,tax.df,tax.rank)
  tax.df<-tax.df%>%
    group_by(Sample)%>%
    fill(agegroup,.direction = "down")%>%
    ungroup()
  return(ggplot_species(taxa.to.plot=taxa.to.plot,
                        tax.df=tax.df,
                        tax.rank=tax.rank,
                        sample.order = sample.order,
                        group.names = group.names,
                        grouping.variable = grouping.variable,
                        metadata.df = metadata.df,
                        ggplot.fill.name = ggplot.fill.name,
                        ggplot.fill.vector = ggplot.fill.vector))
}

# ps.q.agg.genus.relab%>%
#   filter(Genus%in%group1.genera&!is.na(MeanRelativeAbundance))%>%
#   group_by(Genus,agegroup)%>%
#   select(Genus,MeanRelativeAbundance,MeanRelativeAbundanceAgegroup)%>%
#   ungroup()%>%
#   distinct()
# ps.q.agg.genus.relab%>%
#   filter(Genus%in%group1.genera&!is.na(MeanRelativeAbundance))%>%
#   ungroup()%>%
#   select(Genus,MeanRelativeAbundance)%>%
#   distinct()%>%
#   arrange(-MeanRelativeAbundance)%>%
#   pull(Genus) -> group1.genera
# 
# ps.q.agg.genus.relab<-add_zero_rows(group1.genera,ps.q.agg.genus.relab,"Genus")
# 
# ggplot_species(group1.genera,ps.q.agg.genus.relab,"Genus")+
#   ggtitle(paste0("Relative abundance of Group 1 members"))

group1.genera
group1.genera.plot<-show.mean.abundances.and.plot(taxa.to.plot =group1.genera,
                                        tax.df=ps.q.agg.genus.relab,
                                        tax.rank="Genus",
                                        sample.order=sample.levels,
                                        group.names = pretty.level.names,
                                        grouping.variable = "agegroup",
                                        metadata.df = custom.md,
                                        ggplot.fill.name = gg.labs.name,
                                        ggplot.fill.vector = custom.fill)
group1.genera.plot<-group1.genera.plot+
  ggtitle(paste0("Relative abundance of Group 1 members"))
group1.genera.plot # 4 plots


group1.species
group1.species.plot<-show.mean.abundances.and.plot(taxa.to.plot =group1.species,
                                                   tax.df=ps.q.agg.relab,
                                                   tax.rank="Species",
                                                   sample.order=sample.levels,
                                                   group.names = pretty.level.names,
                                                   grouping.variable = "agegroup",
                                                   metadata.df = custom.md,
                                                   ggplot.fill.name = gg.labs.name,
                                                   ggplot.fill.vector = custom.fill)
group1.species.plot<-group1.species.plot+
  ggtitle(paste0("Relative abundance of Group 1 members"))
group1.species.plot # not found

# Group 2 increased genera
group2.increased.genera
group2.increased.genera.plot<-show.mean.abundances.and.plot(taxa.to.plot =group2.increased.genera,
                                                   tax.df =ps.q.agg.genus.relab,
                                                   tax.rank="Genus",
                                                   sample.order=sample.levels,
                                                   group.names = pretty.level.names,
                                                   grouping.variable = "agegroup",
                                                   metadata.df = custom.md,
                                                   ggplot.fill.name = gg.labs.name,
                                                   ggplot.fill.vector = custom.fill)
group2.increased.genera.plot<-group2.increased.genera.plot+
  ggtitle(paste0("Relative abundance of Group 2 members increased with age"))
group2.increased.genera.plot # 4 plots

# Group 2 unhealthy genera
group2.unhealthy.genera
group2.unhealthy.genera.plot<-show.mean.abundances.and.plot(taxa.to.plot=group2.unhealthy.genera,
                                                            tax.df=ps.q.agg.genus.relab,
                                                            tax.rank="Genus",
                                                            sample.order=sample.levels,
                                                            group.names = pretty.level.names,
                                                            grouping.variable = "agegroup",
                                                            metadata.df = custom.md,
                                                            ggplot.fill.name = gg.labs.name,
                                                            ggplot.fill.vector = custom.fill)
group2.unhealthy.genera.plot<-group2.unhealthy.genera.plot+
  ggtitle(paste0("Relative abundance of Group 2 members associated with unhealthy aging"))
group2.unhealthy.genera.plot # 7 plots

# Group 2 unhealthy families
group2.unhealthy.families
group2.unhealthy.families.plot<-show.mean.abundances.and.plot(taxa.to.plot=group2.unhealthy.families,
                                                              tax.df=ps.q.agg.family.relab,
                                                            tax.rank="Family",
                                                            sample.order=sample.levels,
                                                            group.names = pretty.level.names,
                                                            grouping.variable = "agegroup",
                                                            metadata.df = custom.md,
                                                            ggplot.fill.name = gg.labs.name,
                                                            ggplot.fill.vector = custom.fill)
group2.unhealthy.families.plot<-group2.unhealthy.families.plot+
  ggtitle(paste0("Relative abundance of Group 2 members associated with unhealthy aging"))
group2.unhealthy.families.plot # 2 plots

# Group 2 unhealthy species
group2.unhealthy.species
group2.unhealthy.species.plot<-show.mean.abundances.and.plot(taxa.to.plot =group2.unhealthy.species,
                                                             tax.df=ps.q.agg.relab,
                                                              tax.rank="Species",
                                                             sample.order=sample.levels,
                                                             group.names = pretty.level.names,
                                                             grouping.variable = "agegroup",
                                                             metadata.df = custom.md,
                                                             ggplot.fill.name = gg.labs.name,
                                                             ggplot.fill.vector = custom.fill)

group2.unhealthy.species.plot<-group2.unhealthy.species.plot+
  ggtitle(paste0("Relative abundance of Group 2 members associated with unhealthy aging"))
group2.unhealthy.species.plot # 6 plots


# Group 3 genera
group3.genera
group3.genera.plot<-show.mean.abundances.and.plot(taxa.to.plot =group3.genera,
                                                  tax.df = ps.q.agg.genus.relab,
                                                  tax.rank = "Genus",
                                                  sample.order=sample.levels,
                                                  group.names = pretty.level.names,
                                                  grouping.variable = "agegroup",
                                                  metadata.df = custom.md,
                                                  ggplot.fill.name = gg.labs.name,
                                                  ggplot.fill.vector = custom.fill)
group3.genera.plot<-group3.genera.plot+
  ggtitle(paste0("Relative abundance of Group 3 members"))
group3.genera.plot # 5 plots


# Group 3 genera
group3.families.plot<-show.mean.abundances.and.plot(taxa.to.plot = group3.families,
                                                    tax.df = ps.q.agg.family.relab,
                                                  tax.rank = "Family",
                                                  sample.order=sample.levels,
                                                  group.names = pretty.level.names,
                                                  grouping.variable = "agegroup",
                                                  metadata.df = custom.md,
                                                  ggplot.fill.name = gg.labs.name,
                                                  ggplot.fill.vector = custom.fill)
group3.families.plot<-group3.families.plot+
  ggtitle(paste0("Relative abundance of Group 3 members"))
group3.families.plot # 1 plots

for (image_format in c("png","tiff")){
  ggsave(paste0("./images/barplots/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "group1.genera-nmr",
                      sep = "-"),".",image_format),
         plot=group1.genera.plot,
         width = 7000,height = 4000,
         units = "px",dpi=300,device = image_format) # 4 plots
  ggsave(paste0("./images/barplots/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "group2.increased.genera-nmr",
                      sep = "-"),".",image_format),
         plot=group2.increased.genera.plot,
         width = 7000,height = 4000,
         units = "px",dpi=300,device = image_format) # 4 plots
  ggsave(paste0("./images/barplots/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "group2.unhealthy.genera-nmr",
                      sep = "-"),".",image_format),
         plot=group2.unhealthy.genera.plot,
         width = 7000,height = 5000,
         units = "px",dpi=300,device = image_format) # 7 plots
  ggsave(paste0("./images/barplots/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "group2.unhealthy.families-nmr",
                      sep = "-"),".",image_format),
         plot=group2.unhealthy.families.plot,
         width = 7000,height = 2000,
         units = "px",dpi=300,device = image_format)    # 2 plots
  ggsave(paste0("./images/barplots/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "group2.unhealthy.species-nmr",
                      sep = "-"),".",image_format),
         plot=group2.unhealthy.species.plot,
         width = 7000,height = 5000,
         units = "px",dpi=300,device = image_format)    # 6 plots
  ggsave(paste0("./images/barplots/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "group3.genera-nmr",
                      sep = "-"),".",image_format),
         plot= group3.genera.plot ,
         width = 7000,height = 5000,
         units = "px",dpi=300,device = image_format)    # 5 plots
  ggsave(paste0("./images/barplots/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "group3.families-nmr",
                      sep = "-"),".",image_format),
         plot= group3.families.plot ,
         width = 4000,height = 2000,
         units = "px",dpi=300,device = image_format)    # 1 plot
  
}


# Bacterial species ####
# Save the table of bacterial abundances
ps.q.agg.bacteria<-ps.q.agg.relab%>%
  ungroup%>%
  filter(Kingdom=="Bacteria")%>%
  group_by(Phylum)%>%
  distinct(Species,.keep_all = T)%>%
  select(Kingdom:Species,MeanRelativeAbundance)
# write.table(ps.q.agg.bacteria,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  pipeline.name,"bacteria.tsv",sep="-")),
#             row.names = F,sep = "\t")


# Non-bacterial species ####
### 1. Eukaryotes ####
#### 1.1. How many eukaryotes do we have ####
ps.q.agg.relab%>%
  ungroup%>%
  filter(Kingdom=="Eukaryota")%>%
  distinct(Species)%>%
  tally

#### 1.2 Total abundance of eukaryotes in classified data ####
ps.q.agg.relab%>%
  ungroup%>%
  filter(Kingdom=="Eukaryota")%>%
  reframe(TotalKingdom=sum(Abundance)/TotalClass*100)%>%
  distinct

#### 1.3 Eukaryotic phyla ####
ps.q.agg.relab%>%
  ungroup%>%
  filter(Kingdom=="Eukaryota")%>%
  group_by(Phylum)%>%
  distinct(Species,.keep_all = T)%>%
  tally

#### 1.4 Eukaryotic species that aren't Streptophyta ####
ps.q.agg.relab%>%
  ungroup%>%
  filter(Kingdom=="Eukaryota",Phylum!="Streptophyta")%>%
  group_by(Phylum)%>%
  distinct(Species,.keep_all = T)%>%
  select(Kingdom,Phylum,Genus,Species,MeanRelativeAbundance)%>%
  group_by(Phylum)%>%
  mutate(NumSpecies=n())%>%
  arrange(-NumSpecies, Phylum)%>%
  select(-NumSpecies)

# Save the table of eukaryotic abundances
ps.q.agg.eukaryotes<-ps.q.agg.relab%>%
  ungroup%>%
  filter(Kingdom=="Eukaryota")%>%
  group_by(Phylum)%>%
  distinct(Species,.keep_all = T)%>%
  select(Kingdom:Species,MeanRelativeAbundance)
# write.table(ps.q.agg.eukaryotes,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  pipeline.name,"eukaryotes.tsv",sep="-")),
#             row.names = F,sep = "\t")


### 2. Archaea ####
#### 2.1. How many archaea do we have ####
ps.q.agg.relab%>%
  ungroup%>%
  filter(Kingdom=="Archaea")%>%
  distinct(Species)%>%
  tally

#### 2.2 Total abundance of archaea in classified data ####
ps.q.agg.relab%>%
  ungroup%>%
  filter(Kingdom=="Archaea")%>%
  reframe(TotalKingdom=sum(Abundance)/TotalClass*100)%>%
  distinct

#### 2.3 Archaea phyla ####
ps.q.agg.relab%>%
  ungroup%>%
  filter(Kingdom=="Archaea")%>%
  group_by(Phylum)%>%
  distinct(Species,.keep_all = T)%>%
  tally

#### 2.4 Archaea species ###
ps.q.agg.relab%>%
  ungroup%>%
  filter(Kingdom=="Archaea")%>%
  group_by(Phylum)%>%
  distinct(Species,.keep_all = T)%>%
  select(Kingdom,Phylum,Genus,Species,MeanRelativeAbundance)%>%
  group_by(Phylum)%>%
  mutate(NumSpecies=n())%>%
  arrange(-NumSpecies, Phylum)%>%
  select(-NumSpecies)

# Save the table of archaea abundances
ps.q.agg.archaea<-ps.q.agg.relab%>%
  ungroup%>%
  filter(Kingdom=="Archaea")%>%
  group_by(Phylum)%>%
  distinct(Species,.keep_all = T)%>%
  select(Kingdom:Species,MeanRelativeAbundance)
# write.table(ps.q.agg.archaea,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  pipeline.name,"archaea.tsv",sep="-")),
#             row.names = F,sep = "\t")


### 3. Viruses ####
#### 2.1. How many archaea do we have ####
ps.q.agg.relab%>%
  ungroup%>%
  filter(Kingdom=="Viruses")%>%
  distinct(Species)%>%
  tally

#### 2.2 Total abundance of viruses in classified data ####
ps.q.agg.relab%>%
  ungroup%>%
  filter(Kingdom=="Viruses")%>%
  reframe(TotalKingdom=sum(Abundance)/TotalClass*100)%>%
  distinct

#### 2.3 Viral phyla ####
ps.q.agg.relab%>%
  ungroup%>%
  filter(Kingdom=="Viruses")%>%
  group_by(Phylum)%>%
  distinct(Species,.keep_all = T)%>%
  tally

#### 2.4 Viral species ###
ps.q.agg.relab%>%
  ungroup%>%
  filter(Kingdom=="Viruses")%>%
  group_by(Phylum)%>%
  distinct(Species,.keep_all = T)%>%
  select(Kingdom,Phylum,Genus,Species,MeanRelativeAbundance)%>%
  group_by(Phylum)%>%
  mutate(NumSpecies=n())%>%
  arrange(-NumSpecies, Phylum)%>%
  select(-NumSpecies)

# Save the table of archaea abundances
ps.q.agg.viruses<-ps.q.agg.relab%>%
  ungroup%>%
  filter(Kingdom=="Viruses")%>%
  group_by(Phylum)%>%
  distinct(Species,.keep_all = T)%>%
  select(Kingdom:Species,MeanRelativeAbundance)
# write.table(ps.q.agg.viruses,
#             file=file.path(rtables.directory,
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  pipeline.name,"viruses.tsv",sep="-")),
#             row.names = F,sep = "\t")

library(ggrepel)
kingdom.boxplot<-ps.q.agg.relab%>%
  filter(Abundance!=0)%>%
  left_join(custom.md[,c("Sample","age","agegroup")])%>%
  group_by(Sample,Kingdom)%>%
  mutate(TotalRelativeAbundanceK=sum(RelativeAbundance))%>%
  select(Sample,Kingdom,age,agegroup,TotalRelativeAbundanceK)%>%
  mutate("Sample_age"=paste(Sample,age,sep = " ("))%>%
  mutate("Append"=" yrs old)")%>%
  unite("Sample_age",Sample_age:Append,sep = "")%>%
  distinct()%>%
  ggplot(aes(x=agegroup,y=TotalRelativeAbundanceK,fill=factor(agegroup)))+
  geom_boxplot()+
  facet_wrap(~Kingdom,scales = "free",nrow=1)+
  scale_color_manual(breaks = unname(pretty.level.names),
                     labels=unname(pretty.level.names))+
  scale_x_discrete(labels=pretty.level.names,
                   limits=custom.levels)+ # rename boxplot labels (x axis)
  scale_fill_manual(values = custom.fill)+
  labs(y="Relative Abundance (%)",
       x="",
       title="Relative abundance (%) of each kingdom")+
  theme(axis.title.y = element_text(size = 25),
        axis.title = element_text(size = 20),
        axis.text.y = ggtext::element_markdown(size=18),
        axis.text.x = element_text(size=20),
        strip.text.x = ggtext::element_markdown(size=20),
        plot.title = element_text(size = 27),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        legend.position = "none")
kingdom.boxplot.dots<-kingdom.boxplot+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.6)
kingdom.boxplot.labeled<-kingdom.boxplot+
  geom_label_repel(aes(label=Sample_age))


for (image_format in c("png","tiff")){
  ggsave(paste0("./images/boxplots/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "kingdom-boxplot-nmr",
                      sep = "-"),".",image_format),
         plot=kingdom.boxplot,
         width = 4000,height = 2000,
         units = "px",dpi=300,device = image_format)
  ggsave(paste0("./images/boxplots/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "kingdom-boxplot-nmr-dots",
                      sep = "-"),".",image_format),
         plot=kingdom.boxplot.dots,
         width = 4000,height = 2000,
         units = "px",dpi=300,device = image_format)
  ggsave(paste0("./images/boxplots/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "kingdom-boxplot-nmr-labeled",
                      sep = "-"),".",image_format),
         plot=kingdom.boxplot.labeled,
         width = 5000,height = 2000,
         units = "px",dpi=300,device = image_format)
}

kingdom.df<-ps.q.agg.relab%>%
  left_join(custom.md[,c("Sample","age")])%>%
  group_by(Sample,Kingdom)%>%
  mutate(TotalRelativeAbundanceK=sum(RelativeAbundance))%>%
  select(Sample,Kingdom,age,agegroup,TotalRelativeAbundanceK)%>%
  mutate("Sample_age"=paste(Sample,age,sep = " ("))%>%
  mutate("Append"=" yrs old)")%>%
  unite("Sample_age",Sample_age:Append,sep = "")%>%
  distinct()

# Not significant results
for(kingdom.var in c("Bacteria","Eukaryota","Archaea","Viruses")){
  test.df<-kingdom.df%>%
    filter(Kingdom==kingdom.var)
  print(kingdom.var)
  w.test<-pairwise.wilcox.test(test.df$TotalRelativeAbundanceK,
                       test.df$agegroup,
                       p.adjust.method = "BH",
                       exact=FALSE)
  print(w.test)
}
