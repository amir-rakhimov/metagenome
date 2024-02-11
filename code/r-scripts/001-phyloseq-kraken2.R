library(tidyverse)
# install.packages(
#   "microViz",
#   repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
# )
library(microViz)
library(phyloseq)
asvlevel=FALSE
agglom.rank<-"Species"
metadatadir<-paste0("../amplicon_nmr/data/metadata/pooled-metadata/") # directory with metadata

combined.report.filename<-"output/kraken2_pipeline/kraken2_reports/20240122_combined_mpa_clean.tsv"

combined.report<-read.table(combined.report.filename, 
           header = T, 
           sep = "\t",
           fill=TRUE,
           quote = "", 
           comment.char = "@",
           na.strings = c("","NA"))
dim(combined.report)
head(combined.report)
combined.report[is.na(combined.report)]<-"NA"
colnames(combined.report)<-gsub("X20240122_","",colnames(combined.report))
colnames(combined.report)<-
  gsub("_wms_no_minimizer_data.k2report","",colnames(combined.report))

# Create metadata
metadata.filename<-paste0(metadatadir,
                          paste("filenames-single-pooled-raw-supercomp.tsv", sep = "-"))
custom.md<-read.table(metadata.filename, header = T)
colnames(custom.md)[1]<-"Sample" # set the first column name as Sample
custom.md<-custom.md%>%
  filter(Sample%in%colnames(combined.report))%>%
  column_to_rownames(var = "Sample")%>%
  select(-absolute.filepath)
# table(grepl('k__',combined.report$Kingdom))
# table(grepl('p__|NA',combined.report$Phylum))
# table(grepl('c__|NA',combined.report$Class))
# table(grepl('o__|NA',combined.report$Order))
# table(grepl('f__|NA',combined.report$Family))
# table(grepl('g__|NA',combined.report$Genus))
# table(grepl('s__|NA',combined.report$Species))


df.otus<-combined.report%>%
  unite("OTU",Kingdom:Species,sep = "|")%>%
  column_to_rownames("OTU")

df.taxa<-combined.report%>%
  select(Kingdom:Species)%>%
  unite("OTU",Kingdom:Species,sep = "|",remove = FALSE)%>%
  column_to_rownames("OTU")

df.taxa$Kingdom<-
  gsub("k__","",df.taxa$Kingdom)
df.taxa$Phylum<-
  gsub("p__","",df.taxa$Phylum)
df.taxa$Class<-
  gsub("c__","",df.taxa$Class)
df.taxa$Order<-
  gsub("o__","",df.taxa$Order)
df.taxa$Family<-
  gsub("f__","",df.taxa$Family)
df.taxa$Genus<-
  gsub("g__","",df.taxa$Genus)
df.taxa$Species<-
  gsub("s__","",df.taxa$Species)
df.taxa$Kingdom<-
  gsub("EukaryotaFungi","Fungi",df.taxa$Kingdom)
df.taxa$Kingdom<-
  gsub("EukaryotaMetazoa","Metazoa",df.taxa$Kingdom)
df.taxa$Kingdom<-
  gsub("EukaryotaViridiplantae","Viridiplantae",df.taxa$Kingdom)
df.taxa$Kingdom<-
  gsub("VirusesBamfordvirae","Viruses_Bamfordvirae",df.taxa$Kingdom)
df.taxa$Kingdom<-
  gsub("VirusesHelvetiavirae","Viruses_Helvetiavirae",df.taxa$Kingdom)
df.taxa$Kingdom<-
  gsub("VirusesHeunggongvirae","Viruses_Heunggongvirae",df.taxa$Kingdom)
df.taxa$Kingdom<-
  gsub("VirusesLoebvirae","Viruses_Loebvirae",df.taxa$Kingdom)
df.taxa$Kingdom<-
  gsub("VirusesOrthornavirae","Viruses_Orthornavirae",df.taxa$Kingdom)
df.taxa$Kingdom<-
  gsub("VirusesPararnavirae","Viruses_Pararnavirae",df.taxa$Kingdom)
df.taxa$Kingdom<-
  gsub("VirusesSangervirae","Viruses_Sangervirae",df.taxa$Kingdom)
df.taxa$Kingdom<-
  gsub("VirusesShotokuvirae","Viruses_Shotokuvirae",df.taxa$Kingdom)
df.taxa$Kingdom<-
  gsub("VirusesTrapavirae","Viruses_Trapavirae",df.taxa$Kingdom)
df.taxa$Kingdom<-
  gsub("VirusesZilligvirae","Viruses_Zilligvirae",df.taxa$Kingdom)

table(df.taxa$Kingdom)
# Convert into phyloseq taxonomyTable object
df.taxa<-tax_table(as.matrix(df.taxa))

df.taxa.pretty<-tax_fix(df.taxa)
df.taxa<-df.taxa.pretty
rm(df.taxa.pretty)

# Create a phyloseq object
ps.q<-phyloseq(otu_table(df.otus,taxa_are_rows = TRUE),
               tax_table(df.taxa),
               sample_data(custom.md))

if (asvlevel==TRUE){
  ps.q.agg<-ps.q %>%
    psmelt()  # transform the phyloseq object into an R dataframe
}else{
  ps.q.agg<-ps.q %>%
    tax_glom(agglom.rank,NArm = FALSE) %>% # agglomerate by agglom.rank
    psmelt()  # transform the phyloseq object into an R dataframe
}
# Remove entries with zero Abundance
ps.q.agg<-ps.q.agg%>%
  filter(Abundance!=0)


# Add relative abundance column: Abundance divided by total abundance in a sample
ps.q.agg<-ps.q.agg%>%
  group_by(class,Sample)%>%
  mutate(TotalSample=sum(Abundance))%>%
  group_by_at(c("class","Sample",agglom.rank))%>%
  mutate(RelativeAbundance=Abundance/TotalSample*100)


# Sanity check: is total relative abundance of each sample 100%?
ps.q.agg %>%
  group_by(Sample) %>% # Group by sample id
  summarise(RelativeAbundance = sum(RelativeAbundance)) %>% # Sum all abundances
  pull(RelativeAbundance) %>% # get only a vector of numbers
  `==`(100) %>% # compare each value to 100 (TRUE/FALSE)
  all() # get one TRUE/FALSE value

# Group the dataframe by classes (animal hosts)
classcol<-which(colnames(ps.q.agg) =="class")
# find a column by which we agglomerated the dataset
# for ASV, we actually use Species
# We find the column index because we can derive the previous taxonomic 
# rank with a simple agglom.rank.col-1.
if(agglom.rank=="OTU"){
  agglom.rank.col<-which(colnames(ps.q.agg) =="Species")
}else{
  agglom.rank.col<-which(colnames(ps.q.agg) ==agglom.rank)
}

# for each class, we take a agglom.rank, sum its abundances from all samples,
# then take a mean. This will be our MeanRelativeAbundance
if(asvlevel==TRUE){
  ps.q.agg<-ps.q.agg%>%
    group_by(class)%>% # group by class (animal host),
    mutate(TotalClass=sum(Abundance))%>%
    mutate(MeanRelativeAbundance = Abundance/TotalClass*100)
}else{ 
  ps.q.agg<-ps.q.agg%>%
    group_by(class)%>% # group by class,
    # compute MeanRelativeAbundance from Abundance 
    mutate(TotalClass=sum(Abundance))%>%
    mutate(MeanRelativeAbundance = Abundance/TotalClass*100)%>%
    select(-OTU)
}

ps.q.agg<-ps.q.agg%>%
  select(-TotalClass,-TotalSample)

objects.to.keep<-c("agglom.rank","ps.q.agg","asvlevel","custom.md")
objects.to.keep<-which(ls()%in%objects.to.keep)
rm(list = ls()[-objects.to.keep])
# Save the workspace
save.image(paste0("./output/rdafiles/",paste("kraken2",
                                             agglom.rank,
                                             "phyloseq-workspace.RData",sep = "-")))

sessionInfo()