## Import libraries ####
# install.packages(
#   "microViz",
#   repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
# )
# install.packages(c("tidyverse"))
# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install("phyloseq")
library(phyloseq)
library(tidyverse)
library(microViz)

agglom.rank<-"Species"
# If you already created a workspace, just load it
summary.stats.phyloseq.date_time<-"20240619_14_25_09"
load(file.path("./output/rdafiles",paste(
  summary.stats.phyloseq.date_time,"kraken2",agglom.rank,
  "phyloseq-summary-stats",
  "workspace.RData",sep = "-")))

# If not, then create it
# CREATING DATA STARTS VVV####
# agglom.rank<-"Species"
# metadatadir<-paste0("../amplicon_nmr/data/metadata/pooled-metadata/") # directory with metadata
# 
# combined.report.date_time="20240515_10_04_04"
# combined.report.filename<-file.path("output/rtables",
#                                     paste(combined.report.date_time,"combined_report.tsv",sep = "_"))
# combined.report<-read.table(combined.report.filename,
#                             header = T,
#                             sep = "\t",
#                             fill=TRUE,
#                             quote = "",
#                             comment.char = "@",
#                             na.strings = c("","NA"))
# colnames(combined.report)<-gsub("X2D10","2D10",colnames(combined.report))
# colnames(combined.report)<-gsub("X2D14","2D14",colnames(combined.report))
# 
# dim(combined.report)
# head(combined.report)
# 
# df.otus<-combined.report%>%
#   unite("OTU",Kingdom:Species,sep = "|")%>%
#   column_to_rownames("OTU")
# 
# df.otus[is.na(df.otus)]<-0
# 
# df.taxa<-combined.report%>%
#   select(Kingdom:Species)%>%
#   unite("OTU",Kingdom:Species,sep = "|",remove = FALSE)%>%
#   column_to_rownames("OTU")
# 
# df.taxa$Kingdom<-
#   gsub("k__","",df.taxa$Kingdom)
# df.taxa$Phylum<-
#   gsub("p__","",df.taxa$Phylum)
# df.taxa$Class<-
#   gsub("c__","",df.taxa$Class)
# df.taxa$Order<-
#   gsub("o__","",df.taxa$Order)
# df.taxa$Family<-
#   gsub("f__","",df.taxa$Family)
# df.taxa$Genus<-
#   gsub("g__","",df.taxa$Genus)
# df.taxa$Species<-
#   gsub("s__","",df.taxa$Species)
# # Convert into phyloseq taxonomyTable object
# df.taxa<-tax_table(as.matrix(df.taxa))
# 
# df.taxa.pretty<-tax_fix(df.taxa)
# df.taxa<-df.taxa.pretty
# rm(df.taxa.pretty)
# 
# # Create metadata
# metadata.filename<-paste0(metadatadir,
#                           paste("filenames-single-pooled-raw-supercomp.tsv", sep = "-"))
# custom.md<-read.table(metadata.filename, header = T)
# colnames(custom.md)[1]<-"Sample" # set the first column name as Sample
# custom.md<-custom.md%>%
#   filter(Sample%in%colnames(combined.report))%>%
#   column_to_rownames(var = "Sample")%>%
#   select(-absolute.filepath)
# custom.md$class<-as.factor(custom.md$class)
# custom.md$animal<-as.factor(custom.md$animal)
# custom.md$sex<-as.factor(custom.md$sex)
# custom.md$birthday<-as.Date(custom.md$birthday)
# custom.md$age<-year(as.period(interval(custom.md$birthday,as.Date("2023-11-16"))))
# 
# # Create a phyloseq object
# ps.q<-phyloseq(otu_table(df.otus,taxa_are_rows = TRUE),
#                tax_table(df.taxa),
#                sample_data(custom.md))
# 
# ## Create a dataframe of absolute abundances ####
# ### Extract absolute abundances ####
# ps.q.agg<-ps.q %>%
#   psmelt() %>% # transform the phyloseq object into an R dataframe
#   filter(Species!="Homo_sapiens")
# 
# ps.q.agg.phylum<-ps.q %>%
#   tax_glom("Phylum",NArm = FALSE) %>% # agglomerate by agglom.rank
#   psmelt()%>% 
#   filter(OTU!="k__Eukaryota|p__Chordata|c__Mammalia|o__Primates|f__Hominidae|g__Homo|s__Homo_sapiens")
# 
# ps.q.agg.family<-ps.q %>%
#   tax_glom("Family",NArm = FALSE) %>% # agglomerate by agglom.rank
#   psmelt()%>% 
#   filter(OTU!="k__Eukaryota|p__Chordata|c__Mammalia|o__Primates|f__Hominidae|g__Homo|s__Homo_sapiens")
# 
# ps.q.agg.genus<-ps.q %>%
#   tax_glom("Genus",NArm = FALSE) %>% # agglomerate by agglom.rank
#   psmelt()%>% 
#   filter(OTU!="k__Eukaryota|p__Chordata|c__Mammalia|o__Primates|f__Hominidae|g__Homo|s__Homo_sapiens")
# 
# # Remove entries with zero Abundance
# ps.q.agg<-ps.q.agg%>%
#   filter(Abundance!=0)
# ps.q.agg.phylum<-ps.q.agg.phylum%>%
#   filter(Abundance!=0)
# ps.q.agg.family<-ps.q.agg.family%>%
#   filter(Abundance!=0)
# ps.q.agg.genus<-ps.q.agg.genus%>%
#   filter(Abundance!=0)
# 
# save.image(file.path("./output/rdafiles",paste(
#   paste(format(Sys.time(),format="%Y%m%d"),
#         format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#   "kraken2",agglom.rank,
#   "phyloseq-summary-stats",
#   "workspace.RData",sep = "-")))
# #CREATING DATA FINISHED^^^, GETTING SUMMARY STATS STARTS VVV####

## Total ASV/phyla/families/genera per class ####
n.species.per.host<-ps.q.agg%>%
  group_by(class,OTU)%>% # group by class (animal host),
  distinct(OTU)%>%
  group_by(class)%>%
  tally()
n.phylum.per.host<-ps.q.agg.phylum%>%
  group_by(class,Phylum)%>% # group by class (animal host),
  distinct(Phylum)%>%
  group_by(class)%>%
  tally()
n.family.per.host<-ps.q.agg.family%>%
  group_by(class,Family)%>% # group by class (animal host),
  distinct(Family)%>%
  group_by(class)%>%
  tally()
n.genus.per.host<-ps.q.agg.genus%>%
  group_by(class,Genus)%>% # group by class (animal host),
  distinct(Genus)%>%
  group_by(class)%>%
  tally()

## Summary table ####
# Number of samples per host, total reads per host, mean library size,
# sd library size, asv per host, phyla per host, families per host, 
# genera per host
summary.table<-ps.q.agg%>%
  group_by(class)%>%
  mutate(TotalSamplesPerHost=n_distinct(Sample))%>%
  mutate(TotalReadsPerHost=sum(Abundance))%>%
  group_by(Sample)%>%
  mutate(LibrarySize=sum(Abundance))%>%
  distinct(Sample,.keep_all = T)%>%
  group_by(class)%>%
  mutate(MeanLibrarySize =round(mean(LibrarySize)),
         SDLibrarySize=round(sd(LibrarySize)))%>%
  select(class,
         TotalSamplesPerHost,
         TotalReadsPerHost,
         MeanLibrarySize,
         SDLibrarySize)%>%
  distinct(class,.keep_all = T)%>%
  arrange(class)%>%
  left_join(n.species.per.host)%>%
  rename(SpeciesPerHost=n)%>%
  left_join(n.phylum.per.host)%>%
  rename(PhylaPerHost=n)%>%
  left_join(n.family.per.host)%>%
  rename(FamiliesPerHost=n)%>%
  left_join(n.genus.per.host)%>%
  rename(GeneraPerHost=n)

# write.table(summary.table,
#             file=file.path("./output/rtables",
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "kraken2-summary-table.tsv",sep="-")),
#             row.names = F,sep = "\t")



### Check how many taxa are unclassified in each NMR sample ####
all.ranks<-c("Kingdom", "Phylum", "Class", "Order", "Family","Genus","Species")
agglom.rank<-"Species"
ps.q.agg<-ps.q.agg%>%
  group_by(class,Sample)%>%
  mutate(TotalSample=sum(Abundance))%>%
  group_by_at(c("class","Sample",agglom.rank))%>%
  mutate(RelativeAbundance=Abundance/TotalSample*100)
ps.q.agg.genus<-ps.q.agg.genus%>%
  group_by(class,Sample)%>%
  mutate(TotalSample=sum(Abundance))%>%
  group_by(class,Sample,Genus)%>%
  mutate(RelativeAbundance=Abundance/TotalSample*100)


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
ps.q.agg.dominant.phyla<-ps.q.agg.phylum%>%
  filter(class=="NMR")%>%
  mutate(TotalClass=sum(Abundance))%>%
  group_by(Phylum)%>%
  mutate(TotalAgglomRank=sum(Abundance))%>%
  mutate(MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)%>%
  select(Phylum,MeanRelativeAbundance)%>%
  distinct()%>%
  arrange(-MeanRelativeAbundance)
# write.table(ps.q.agg.dominant.phyla,
#             file=file.path("./output/rtables",
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "kraken2-dominant-phyla.tsv",sep="-")),
#             row.names = F,sep = "\t")

### Most abundant non-bacterial phyla  ####
ps.q.agg.dominant.phyla.nonbact<-ps.q.agg.phylum%>%
  filter(class=="NMR")%>%
  mutate(TotalClass=sum(Abundance))%>%
  group_by(Phylum)%>%
  mutate(TotalAgglomRank=sum(Abundance))%>%
  mutate(MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)%>%
  filter(Kingdom!="Bacteria")%>%
  select(Phylum,MeanRelativeAbundance)%>%
  distinct()%>%
  arrange(-MeanRelativeAbundance)
# write.table(ps.q.agg.dominant.phyla.nonbact,
#             file=file.path("./output/rtables",
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "kraken2-dominant-phyla-nonbact.tsv",sep="-")),
#             row.names = F,sep = "\t")

### Most abundant families ####
ps.q.agg.dominant.families<-ps.q.agg.family%>%
  filter(class=="NMR")%>%
  mutate(TotalClass=sum(Abundance))%>%
  group_by(Family)%>%
  mutate(TotalAgglomRank=sum(Abundance))%>%
  mutate(MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)%>%
  select(Phylum,Family,MeanRelativeAbundance)%>%
  distinct()%>%
  arrange(-MeanRelativeAbundance)
# write.table(ps.q.agg.dominant.families,
#             file=file.path("./output/rtables",
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "kraken2-dominant-families.tsv",sep="-")),
#             row.names = F,sep = "\t")

### Most abundant non-bacterial families ####
ps.q.agg.dominant.families.nonbact<-ps.q.agg.family%>%
  filter(class=="NMR")%>%
  mutate(TotalClass=sum(Abundance))%>%
  group_by(Family)%>%
  mutate(TotalAgglomRank=sum(Abundance))%>%
  mutate(MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)%>%
  filter(Kingdom!="Bacteria")%>%
  select(Phylum,Family,MeanRelativeAbundance)%>%
  distinct()%>%
  arrange(-MeanRelativeAbundance)
# write.table(ps.q.agg.dominant.families.nonbact,
#             file=file.path("./output/rtables",
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "kraken2-dominant-families-nonbact.tsv",sep="-")),
#             row.names = F,sep = "\t")
### Most abundant genera ####
ps.q.agg.dominant.genera<-ps.q.agg.genus%>%
  ungroup()%>%
  filter(class=="NMR")%>%
  mutate(TotalClass=sum(Abundance))%>%
  group_by(Genus)%>%
  mutate(TotalAgglomRank=sum(Abundance))%>%
  mutate(MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)%>%
  select(Phylum,Genus,MeanRelativeAbundance)%>%
  distinct()%>%
  arrange(-MeanRelativeAbundance)
# write.table(ps.q.agg.dominant.genera,
#             file=file.path("./output/rtables",
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "kraken2-dominant-genera.tsv",sep="-")),
#             row.names = F,sep = "\t")

### Most abundant non-bacterial genera ####
ps.q.agg.dominant.genera.nonbact<-ps.q.agg.genus%>%
  ungroup()%>%
  filter(class=="NMR")%>%
  mutate(TotalClass=sum(Abundance))%>%
  group_by(Genus)%>%
  mutate(TotalAgglomRank=sum(Abundance))%>%
  mutate(MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)%>%
  filter(Kingdom!="Bacteria")%>%
  select(Phylum,Genus,MeanRelativeAbundance)%>%
  distinct()%>%
  arrange(-MeanRelativeAbundance)
# write.table(ps.q.agg.dominant.genera.nonbact,
#             file=file.path("./output/rtables",
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "kraken2-dominant-genera-nonbact.tsv",sep="-")),
#             row.names = F,sep = "\t")

### Most abundant species overall ####
ps.q.agg.dominant.species<-ps.q.agg%>%
  ungroup()%>%
  filter(class=="NMR")%>%
  mutate(TotalClass=sum(Abundance))%>%
  group_by(Species)%>%
  mutate(TotalAgglomRank=sum(Abundance))%>%
  mutate(MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)%>%
  select(Species,MeanRelativeAbundance)%>%
  distinct()%>%
  arrange(-MeanRelativeAbundance)
# write.table(ps.q.agg.dominant.species,
#             file=file.path("./output/rtables",
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "kraken2-dominant-species-all.tsv",sep="-")),
#             row.names = F,sep = "\t")
### Most abundant non-bacterial species ####
ps.q.agg.dominant.species.nonbact<-ps.q.agg%>%
  ungroup()%>%
  filter(class=="NMR")%>%
  mutate(TotalClass=sum(Abundance))%>%
  group_by(Species)%>%
  mutate(TotalAgglomRank=sum(Abundance))%>%
  mutate(MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)%>%
  filter(Kingdom!="Bacteria")%>%
  select(Species,MeanRelativeAbundance)%>%
  distinct()%>%
  arrange(-MeanRelativeAbundance)
# write.table(ps.q.agg.dominant.species.nonbact,
#             file=file.path("./output/rtables",
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "kraken2-dominant-species-nonbact.tsv",sep="-")),
#             row.names = F,sep = "\t")


## Add relative abundance, mean relative abundance, total abundance, etc ####
ps.q.agg.family.relab<-ps.q.agg.family%>%
  group_by(class,Sample)%>%
  mutate(TotalSample=sum(Abundance))%>%
  group_by_at(c("class","Sample","Family"))%>%
  mutate(RelativeAbundance=Abundance/TotalSample*100)%>%
  group_by(class)%>% # group by class (animal host),
  mutate(TotalClass=sum(Abundance))%>%
  group_by_at(c("class","Family"))%>%
  mutate(TotalAgglomRank=sum(Abundance))%>%
  mutate(MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)

ps.q.agg.genus.relab<-ps.q.agg.genus%>%
  group_by(class,Sample)%>%
  mutate(TotalSample=sum(Abundance))%>%
  group_by_at(c("class","Sample","Genus"))%>%
  mutate(RelativeAbundance=Abundance/TotalSample*100)%>%
  group_by(class)%>% # group by class (animal host),
  mutate(TotalClass=sum(Abundance))%>%
  group_by_at(c("class","Genus"))%>%
  mutate(TotalAgglomRank=sum(Abundance))%>%
  mutate(MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)
  
ps.q.agg.relab<-ps.q.agg%>%
  group_by(class,Sample)%>%
  mutate(TotalSample=sum(Abundance))%>%
  group_by_at(c("class","Sample","Species"))%>%
  mutate(RelativeAbundance=Abundance/TotalSample*100)%>%
  group_by(class)%>% # group by class (animal host),
  mutate(TotalClass=sum(Abundance))%>%
  group_by_at(c("class","Species"))%>%
  mutate(TotalAgglomRank=sum(Abundance))%>%
  mutate(MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)
  
## Add agegroup (Must run for plotting) ####
custom.md<-readRDS("./output/rdafiles/custom.md.ages.rds")
# we create these new levels for plots
ps.q.agg.relab <- ps.q.agg.relab %>%
  left_join(custom.md)%>%
  group_by(agegroup)%>% # group by class (animal host),
  mutate(TotalAgegroup=sum(Abundance))%>%
  group_by(agegroup,Species)%>%
  mutate(TotalAgglomRankAge=sum(Abundance))%>%
  mutate(MeanRelativeAbundanceAgegroup=TotalAgglomRankAge/TotalAgegroup*100)
ps.q.agg.genus.relab<-ps.q.agg.genus.relab%>%
  left_join(custom.md)%>%
  group_by(agegroup)%>% # group by class (animal host),
  mutate(TotalAgegroup=sum(Abundance))%>%
  group_by(agegroup,Genus)%>%
  mutate(TotalAgglomRankAge=sum(Abundance))%>%
  mutate(MeanRelativeAbundanceAgegroup=TotalAgglomRankAge/TotalAgegroup*100)
ps.q.agg.family.relab<-ps.q.agg.family.relab%>%
  left_join(custom.md)%>%
  group_by(agegroup)%>% # group by class (animal host),
  mutate(TotalAgegroup=sum(Abundance))%>%
  group_by(agegroup,Family)%>%
  mutate(TotalAgglomRankAge=sum(Abundance))%>%
  mutate(MeanRelativeAbundanceAgegroup=TotalAgglomRankAge/TotalAgegroup*100)


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
#             file=file.path("./output/rtables",
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "kraken2","bacteroidaceae-nmr-table.tsv",sep="-")),
#             row.names = F,sep = "\t")

### What are the most dominant Bacteroidota families in NMR ####
bacteroidota.nmr<-ps.q.agg.family.relab%>%
  filter(Phylum=="Bacteroidota",class=="NMR")%>%
  group_by(Family)%>%
  distinct(Family,.keep_all = T)%>%
  arrange(-MeanRelativeAbundance)%>%
  select(Family,MeanRelativeAbundance)
# write.table(bacteroidota.nmr,
#             file=file.path("./output/rtables",
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "kraken2","bacteroidota-nmr-table.tsv",sep="-")),
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
# write.table(spirochaetaceae.nmr,
#             file=file.path("./output/rtables",
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "kraken2", "spirochaetaceae-nmr-table.tsv",sep="-")),
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
# write.table(treponema.nmr,
#             file=file.path("./output/rtables",
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "kraken2","treponema-nmr-table.tsv",sep="-")),
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
#             file=file.path("./output/rtables",
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                 "kraken2", "mogibacteriaceae_anaerovoracaceae-all-table.tsv",sep="-")),
#             row.names = F,sep = "\t")



## Setup plots ####
library(Polychrome)
library(ggtext)

pretty.level.names<-names(table(ps.q.agg.relab$old_agegroup))
names(pretty.level.names)<-names(table(ps.q.agg.relab$agegroup))
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

sample.levels<-ps.q.agg.relab%>%
  ungroup()%>%
  select(Sample,agegroup)%>%
  arrange(agegroup)%>%
  distinct()
sample.levels$Sample<-factor(sample.levels$Sample,
                             levels=unique(sample.levels$Sample))

all.ranks<-c("Kingdom", "Phylum", "Class", "Order", "Family","Genus","Species")

## Function to add zero rows if a species doesn't exist in a sample ####
# 1. Go through each sample
# 2. Find which organisms in the taxa.vector are missing in the sample
# 3. Create a tibble (new row)
# In order to have column names that are stored as string we make use 
# of bang bang operator !! which forces the evaluation of it succeeding name
# We also need to use walrus := instead of = which are equivalent and 
# prompts you to supply name (as is the case with our variable name) on it LHS (left hand side)
# tibble(Sample = "sample", 
#        Species =taxa,
#        Abundance = 0,
#        RelativeAbundance=0,
#        !!agglom.rank:="Treponema")
# 4. Add the row at the end of the dataframe
# 5. Fill in NA values in the empty columns based on non-empty rows in the Sample column 

add_zero_rows<-function(taxa.vector,tax.table,tax.rank){
  tax.table<-tax.table%>%ungroup()
  for (sample.name in unique(tax.table$Sample)) {
    missing_taxa <- setdiff(taxa.vector, tax.table[[tax.rank]][tax.table$Sample == sample.name])
    if (length(missing_taxa) > 0) {
      for (taxa in missing_taxa) {
        if(tax.rank=="Species"){
          new_row <- tibble(Sample = sample.name, 
                            Species= taxa, 
                            Abundance = 0,
                            RelativeAbundance=0,
                            Genus=unique(pull(ps.q.agg.relab[which(ps.q.agg.relab$Species==taxa),"Genus"])))
        }else{
          new_row <- tibble(Sample = sample.name, 
                            !!tax.rank:= taxa, 
                            Abundance = 0,
                            RelativeAbundance=0)
        }
        
        tax.table <- tax.table %>% add_row(.before = nrow(df), !!!new_row)
      }
    }
  }
  # To fill the NA values in the empty columns based on non-empty rows in the Sample column 
  tax.table<- tax.table %>%
    group_by(Sample) %>%
    fill(age, .direction = "down")%>%
    fill(agegroup, .direction = "down")%>%
    fill(old_agegroup, .direction = "down")
  
}


# If we are plotting species from a vector of names, use filter(Species%in%species.to.plot)
# If we are plotting all species from a bigger group like entire genus, use
# filter(get(agglom.rank)==species.to.plot)
ggplot.species<-function(taxa.to.plot,
                         tax.table,
                         tax.rank){
  # if(tax.rank=="Species"){
  #   ggplot.object<-tax.table%>%
  #     filter(Species%in%taxa.to.plot)
  # }else{
  #   ggplot.object<-tax.table%>%
  #     filter(get(tax.rank)==taxa.to.plot)
  # }
  ggplot.object<-tax.table%>%
    filter(get(tax.rank)%in%taxa.to.plot)

  # order the plot by species vector
  taxa.to.plot<-gsub("_"," ",taxa.to.plot)
  taxa.to.plot<-paste0("<i>",taxa.to.plot,"</i>")
  
  ggplot.object<-ggplot.object%>%
    mutate(!!tax.rank:=gsub("_"," ",get(tax.rank)),
           !!tax.rank:=paste0("<i>",get(tax.rank),"</i>"),
           !!tax.rank:=factor(get(tax.rank),levels=taxa.to.plot))%>%
    mutate(Sample=factor(Sample,levels=sample.levels$Sample))%>%
    group_by_at(c("class",tax.rank))%>%
    ggplot(aes(x=Sample,
               y=RelativeAbundance,
               fill=factor(agegroup)))+
    geom_bar(stat="identity")+
    facet_wrap(~get(tax.rank),
               scales = "free",
               ncol = 2)+
    theme_bw()+
    labs(x="",
         y="Relative abundance (%)",
         fill=gg.labs.name)+
    scale_color_manual(breaks = unname(pretty.level.names),
                       labels=unname(pretty.level.names))+
    scale_x_discrete(labels=pretty.level.names,
                     limits=sample.levels$Sample)+
    scale_fill_manual(values = custom.fill,
                      labels=pretty.level.names)+
    theme(axis.title.y = element_text(size = 25),
          axis.title = element_text(size = 20),
          axis.text.y = ggtext::element_markdown(size=18),
          axis.text.x = element_text(size=20),
          strip.text.x = ggtext::element_markdown(size=20),
          plot.title = element_text(size = 27),
          legend.text = element_text(size = 20),
          legend.title = element_text(size = 25),
          legend.position = "right")
  return(ggplot.object)
}



## Sulfur metabolising bacteria ####
# We don't know where the "sulf" pattern is (Phylum, Order, Genus, or anywhere else),
# so we use sapply with grepl
# https://stackoverflow.com/questions/47941680/grepl-across-multiple-specified-columns
sulfur.bact.nmr<-
  ps.q.agg.relab[!!rowSums(sapply(ps.q.agg.relab[,all.ranks], grepl, pattern = "sulf") ),]%>%
  group_by(Species)%>%
  mutate(min=min(RelativeAbundance),
         max=max(RelativeAbundance),
         mean=TotalAgglomRank/TotalClass*100,
         sd=sd(RelativeAbundance),
         n=n())%>%
  select(Species,MeanRelativeAbundance,sd)%>%
  distinct()%>%
  arrange(-MeanRelativeAbundance)

# write.table(sulfur.bact.nmr,
#             file=file.path("./output/rtables",
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "kraken2","sulfur.bact.nmr-nmr-table.tsv",sep="-")),
#             row.names = F,sep = "\t")

### Plot sulfur metabolising bacteria ####
# Add zero rows
# ps.q.agg.relab<-ps.q.agg.relab%>%ungroup()
ps.q.agg.relab<-add_zero_rows(unique(sulfur.bact.nmr$Species),ps.q.agg.relab,"Species")

ggplot.species(sulfur.bact.nmr$Species,ps.q.agg.relab,"Species")+
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

ggplot.species(treponema.sp.nmr,ps.q.agg.relab,"Species")+
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

ps.q.agg.relab%>%
  filter(Species%in%debebe.species&!is.na(MeanRelativeAbundance))%>%
  ungroup()%>%
  select(Species,MeanRelativeAbundance)%>%
  distinct()


ggplot.species(debebe.species,ps.q.agg.relab,"Species")+
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

ggplot.species(debebe16s.families,ps.q.agg.family.relab,"Family")+
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

my.species.plot<-ggplot.species(species.to.plot,ps.q.agg.relab,"Species")+
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
                                        tax.table,
                                        tax.rank){
  tax.table%>%
    filter(get(tax.rank)%in%taxa.to.plot&!is.na(MeanRelativeAbundance))%>%
    group_by_at(c(tax.rank,"agegroup"))%>%
    select(all_of(c(tax.rank,"MeanRelativeAbundance","MeanRelativeAbundanceAgegroup")))%>%
    ungroup()%>%
    distinct()%>%
    print()
  tax.table%>%
    filter(get(tax.rank)%in%taxa.to.plot&!is.na(MeanRelativeAbundance))%>%
    ungroup()%>%
    select(all_of(c(tax.rank,"MeanRelativeAbundance")))%>%
    distinct()%>%
    arrange(-MeanRelativeAbundance)%>%
    pull(tax.rank) -> taxa.to.plot
  tax.table<-add_zero_rows(taxa.to.plot,tax.table,tax.rank)
  
  return(ggplot.species(taxa.to.plot,tax.table,tax.rank))
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
# ggplot.species(group1.genera,ps.q.agg.genus.relab,"Genus")+
#   ggtitle(paste0("Relative abundance of Group 1 members"))

group1.genera
group1.genera.plot<-show.mean.abundances.and.plot(group1.genera,
                                        ps.q.agg.genus.relab,
                                        "Genus")
group1.genera.plot<-group1.genera.plot+
  ggtitle(paste0("Relative abundance of Group 1 members"))
group1.genera.plot # 4 plots


group1.species
group1.species.plot<-show.mean.abundances.and.plot(group1.species,
                                                  ps.q.agg.relab,
                                                  "Species")
group1.species.plot<-group1.species.plot+
  ggtitle(paste0("Relative abundance of Group 1 members"))
group1.species.plot # not found

# Group 2 increased genera
group2.increased.genera
group2.increased.genera.plot<-show.mean.abundances.and.plot(group2.increased.genera,
                                                   ps.q.agg.genus.relab,
                                                   "Genus")
group2.increased.genera.plot<-group2.increased.genera.plot+
  ggtitle(paste0("Relative abundance of Group 2 members increased with age"))
group2.increased.genera.plot # 4 plots

# Group 2 unhealthy genera
group2.unhealthy.genera
group2.unhealthy.genera.plot<-show.mean.abundances.and.plot(group2.unhealthy.genera,
                                                            ps.q.agg.genus.relab,
                                                            "Genus")
group2.unhealthy.genera.plot<-group2.unhealthy.genera.plot+
  ggtitle(paste0("Relative abundance of Group 2 members associated with unhealthy aging"))
group2.unhealthy.genera.plot # 7 plots

# Group 2 unhealthy families
group2.unhealthy.families
group2.unhealthy.families.plot<-show.mean.abundances.and.plot(group2.unhealthy.families,
                                                              ps.q.agg.family.relab,
                                                            "Family")
group2.unhealthy.families.plot<-group2.unhealthy.families.plot+
  ggtitle(paste0("Relative abundance of Group 2 members associated with unhealthy aging"))
group2.unhealthy.families.plot # 2 plots

# Group 2 unhealthy species
group2.unhealthy.species
group2.unhealthy.species.plot<-show.mean.abundances.and.plot(group2.unhealthy.species,
                                                             ps.q.agg.relab,
                                                              "Species")

group2.unhealthy.species.plot<-group2.unhealthy.species.plot+
  ggtitle(paste0("Relative abundance of Group 2 members associated with unhealthy aging"))
group2.unhealthy.species.plot # 6 plots


# Group 3 genera
group3.genera
group3.genera.plot<-show.mean.abundances.and.plot(group3.genera,
                                                  ps.q.agg.genus.relab,
                                                             "Genus")
group3.genera.plot<-group3.genera.plot+
  ggtitle(paste0("Relative abundance of Group 3 members"))
group3.genera.plot # 5 plots


# Group 3 genera
group3.families.plot<-show.mean.abundances.and.plot(group3.families,
                                                    ps.q.agg.family.relab,
                                                  "Family")
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
                      " group3.genera-nmr",
                      sep = "-"),".",image_format),
         plot= group3.genera.plot ,
         width = 7000,height = 5000,
         units = "px",dpi=300,device = image_format)    # 5 plots
  ggsave(paste0("./images/barplots/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      " group3.families-nmr",
                      sep = "-"),".",image_format),
         plot= group3.families.plot ,
         width = 4000,height = 2000,
         units = "px",dpi=300,device = image_format)    # 1 plot
  
}



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
#             file=file.path("./output/rtables",
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "kraken2-eukaryotes.tsv",sep="-")),
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
#             file=file.path("./output/rtables",
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "kraken2-archaea.tsv",sep="-")),
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
#             file=file.path("./output/rtables",
#                            paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                                  "kraken2-viruses.tsv",sep="-")),
#             row.names = F,sep = "\t")

library(ggrepel)
kingdom.boxplot<-ps.q.agg.relab%>%
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
