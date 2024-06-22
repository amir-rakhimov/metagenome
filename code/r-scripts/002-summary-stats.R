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

# Total ASV/phyla/families/genera per class ####
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

# Summary table ####
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



### 8. Check how many taxa are unclassified in each NMR sample ####
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


# Check the number of distinct OTU/taxa
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

# Most abundant phyla ####
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

# Most abundant non-bacterial phyla  ####
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

# Most abundant families ####
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

# Most abundant non-bacterial families ####
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
# Most abundant genera ####
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

# Most abundant non-bacterial genera ####
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

# Most abundant species overall ####
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
# Most abundant non-bacterial species ####
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


# Compare with Debebe et al. ####
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
  
# How much Bacteroidaceae are in NMR  ####
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

# What are the most dominant Bacteroidota families in NMR ####
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

# Spirochaetaceae and Treponema ####
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

# How many ASVs in Treponema ####
ps.q.agg%>%
  ungroup()%>%
  filter(Genus=="Treponema",class=="NMR")%>%
  distinct(OTU)%>%
  tally

# Mogibacteriaceae is renamed to Anaerovoracaceae ####
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

# sulfur metabolising bacteria ####
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


# Plot sulfur metabolising bacteria ####
library(Polychrome)
library(ggtext)
ps.q.agg.relab<-ps.q.agg.relab%>%
  mutate(agegroup=cut(age, breaks =c(0,10,16),
                      right = FALSE))
# we create these new levels because maaslin is itsy bitsy
unique_levels <- ps.q.agg.relab %>%
  ungroup()%>%
  distinct(agegroup)%>%
  arrange(agegroup) %>%
  mutate(new_agegroup = paste0("agegroup", agegroup))%>%
  mutate(new_agegroup = gsub("\\(|\\)|\\[|\\]","",new_agegroup))%>%
  mutate(new_agegroup = gsub("\\,","_",new_agegroup))
ps.q.agg.relab <- ps.q.agg.relab %>%
  left_join(unique_levels, by = "agegroup")
colnames(ps.q.agg.relab)[which(colnames(ps.q.agg.relab)=="agegroup")]<-"old_agegroup"
colnames(ps.q.agg.relab)[which(colnames(ps.q.agg.relab)=="new_agegroup")]<-"agegroup"
# add age group to metadata
custom.md$Sample<-rownames(custom.md)
custom.md<-custom.md%>% 
  group_by(Sample)%>%
  mutate(birthday=as.Date(birthday))%>%
  mutate(age=year(as.period(interval(birthday,as.Date("2023-11-16")))))%>%
  left_join(unique(ps.q.agg.relab[,c("Sample","agegroup")]),by="Sample")%>%
  as.data.frame()
rownames(custom.md)<-custom.md$Sample



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

# Add zero rows
ps.q.agg.relab<-ps.q.agg.relab%>%ungroup()
for (sample in unique(ps.q.agg.relab$Sample)) {
  missing_species <- setdiff(sulfur.bact.nmr$Species, ps.q.agg.relab$Species[ps.q.agg.relab$Sample == sample])
  if (length(missing_species) > 0) {
    for (species in missing_species) {
      new_row <- tibble(Sample = sample, 
                        Species = species, 
                        Abundance = 0,
                        RelativeAbundance=0)
      ps.q.agg.relab <- ps.q.agg.relab %>% add_row(.before = nrow(df), !!!new_row)
    }
  }
}
# To fill the NA values in the empty columns based on non-empty rows in the Sample column 
ps.q.agg.relab<- ps.q.agg.relab %>%
  group_by(Sample) %>%
  fill(age, .direction = "down")%>%
  fill(agegroup, .direction = "down")%>%
  fill(old_agegroup, .direction = "down")


sample.levels<-ps.q.agg.relab%>%
  select(Sample,agegroup)%>%
  arrange(agegroup)%>%
  distinct()
sample.levels$Sample<-factor(sample.levels$Sample,
                             levels=sample.levels$Sample)
                 
ps.q.agg.relab%>%
  filter(Species%in%sulfur.bact.nmr$Species)%>%
  mutate(#`Age group`=old_agegroup,
         Species=gsub("_"," ",Species),
         Species=paste0("<i>",Species,"</i>"))%>%
  mutate(Sample=factor(Sample,levels=sample.levels$Sample))%>%
  group_by(class,Species)%>%
  ggplot(aes(x=Sample,
             y=RelativeAbundance,
             fill=factor(agegroup)))+
  geom_bar(stat="identity")+
  facet_wrap(~Species,
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
  theme(plot.margin=unit(c(1,1,1,2), 'cm'),
        axis.title = element_text(size = 20),
        axis.title.y = element_text(size = 25),
        axis.text.y = ggtext::element_markdown(size=18),
        axis.text.x = element_text(size=20),
        strip.text.x = ggtext::element_markdown(size=20),
        plot.title = element_text(size = 27),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        legend.position = "right")+
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


# Plot Treponema ####
treponema.sp.nmr<-ps.q.agg.relab%>%
  filter(Genus=="Treponema",class=="NMR")%>%
  select(Species)
ps.q.agg.relab<-ps.q.agg.relab%>%ungroup()
for (sample in unique(ps.q.agg.relab$Sample)) {
  missing_species <- setdiff(treponema.sp.nmr$Species, ps.q.agg.relab$Species[ps.q.agg.relab$Sample == sample])
  if (length(missing_species) > 0) {
    for (species in missing_species) {
      new_row <- tibble(Sample = sample, 
                        Species = species, 
                        Abundance = 0,
                        RelativeAbundance=0,
                        Genus="Treponema")
      ps.q.agg.relab <- ps.q.agg.relab %>% add_row(.before = nrow(df), !!!new_row)
    }
  }
}
# To fill the NA values in the empty columns based on non-empty rows in the Sample column 
ps.q.agg.relab<- ps.q.agg.relab %>%
  group_by(Sample) %>%
  fill(age, .direction = "down")%>%
  fill(agegroup, .direction = "down")%>%
  fill(old_agegroup, .direction = "down")

ps.q.agg.relab%>%
  filter(Genus=="Treponema")%>%
  mutate(#`Age group`=old_agegroup,
    Species=gsub("_"," ",Species),
    Species=paste0("<i>",Species,"</i>"))%>%
  mutate(Sample=factor(Sample,levels=sample.levels$Sample))%>%
  group_by(class,Species)%>%
  ggplot(aes(x=Sample,
             y=RelativeAbundance,
             fill=factor(agegroup)))+
  geom_bar(stat="identity")+
  facet_wrap(~Species,
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
  theme(plot.margin=unit(c(1,1,1,2), 'cm'),
        axis.title.y = element_text(size = 25),
        axis.title = element_text(size = 20),
        axis.text.y = ggtext::element_markdown(size=18),
        axis.text.x = element_text(size=20),
        strip.text.x = ggtext::element_markdown(size=20),
        plot.title = element_text(size = 27),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        legend.position = "right")+
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
