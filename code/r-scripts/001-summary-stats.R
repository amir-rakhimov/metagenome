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
phyloseq.date_time<-"20240517_14_03_12"
load(file.path("./output/rdafiles",paste(
  phyloseq.date_time,"kraken2",agglom.rank,
  "phyloseq-summary-stats",
  "workspace.RData",sep = "-")))

# If not, then create it
# CREATING DATA STARTS VVV####
# agglom.rank<-"Species"
# metadatadir<-paste0("../amplicon_nmr/data/metadata/pooled-metadata/") # directory with metadata
# 
# date_time="20240518_13_40_08"
# combined.report.filename<-file.path("output/rtables",
#                                     paste(date_time,"combined_report.tsv",sep = "_"))
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
# custom.md$age<-year(as.period(interval(custom.md$birthday,now())))
# 
# # Create a phyloseq object
# ps.q<-phyloseq(otu_table(df.otus,taxa_are_rows = TRUE),
#                tax_table(df.taxa),
#                sample_data(custom.md))
# 
# ## Create a dataframe of absolute abundances ####
# ### Extract absolute abundances ####
# ps.q.agg<-ps.q %>%
#   psmelt()  # transform the phyloseq object into an R dataframe
# 
# ps.q.agg.phylum<-ps.q %>%
#   tax_glom("Phylum",NArm = FALSE) %>% # agglomerate by agglom.rank
#   psmelt()
# ps.q.agg.family<-ps.q %>%
#   tax_glom("Family",NArm = FALSE) %>% # agglomerate by agglom.rank
#   psmelt()
# ps.q.agg.genus<-ps.q %>%
#   tax_glom("Genus",NArm = FALSE) %>% # agglomerate by agglom.rank
#   psmelt()
# 
# # Remove entries with zero Abundance
# ps.q.agg<-ps.q.agg%>%
#   filter(Abundance!=0)
# ps.q.agg.phylum<-ps.q.agg.phylum%>%
#   filter(Abundance!=0)
# ps.q.agg.family<-ps.q.agg.family%>%
#   filter(Abundance!=0)
# 
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
write.table(ps.q.agg.dominant.species,
            file=file.path("./output/rtables",
                           paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "kraken2-dominant-species-all.tsv",sep="-")),
            row.names = F,sep = "\t")
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
write.table(ps.q.agg.dominant.species.nonbact,
            file=file.path("./output/rtables",
                           paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "kraken2-dominant-species-nonbact.tsv",sep="-")),
            row.names = F,sep = "\t")

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
write.table(ps.q.agg.dominant.phyla,
            file=file.path("./output/rtables",
                           paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "kraken2-dominant-phyla.tsv",sep="-")),
            row.names = F,sep = "\t")

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
write.table(ps.q.agg.dominant.phyla.nonbact,
            file=file.path("./output/rtables",
                           paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "kraken2-dominant-phyla-nonbact.tsv",sep="-")),
            row.names = F,sep = "\t")

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
write.table(ps.q.agg.dominant.families,
            file=file.path("./output/rtables",
                           paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "kraken2-dominant-families.tsv",sep="-")),
            row.names = F,sep = "\t")

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
write.table(ps.q.agg.dominant.families.nonbact,
            file=file.path("./output/rtables",
                           paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "kraken2-dominant-families-nonbact.tsv",sep="-")),
            row.names = F,sep = "\t")
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
write.table(ps.q.agg.dominant.genera,
            file=file.path("./output/rtables",
                           paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "kraken2-dominant-genera.tsv",sep="-")),
            row.names = F,sep = "\t")

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
write.table(ps.q.agg.dominant.genera.nonbact,
            file=file.path("./output/rtables",
                           paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "kraken2-dominant-genera-nonbact.tsv",sep="-")),
            row.names = F,sep = "\t")

