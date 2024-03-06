library(tidyverse)
# install.packages(
#   "microViz",
#   repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
# )
library(microViz)
library(phyloseq)
agglom.rank<-"Species"
metadatadir<-paste0("../amplicon_nmr/data/metadata/pooled-metadata/") # directory with metadata

bracken_krona_txt.dir<-"output/kraken2_pipeline/bracken_krona_txt"
bracken_krona_txt.files<-
  list.files(bracken_krona_txt.dir)[grep("wms_table.tsv",list.files(bracken_krona_txt.dir))]

list.tables<-list()


for(filename in bracken_krona_txt.files){
  sample.base.name<-gsub("_wms_table.tsv","",filename)
  bracken_krona_table<-
    read.table(file.path(bracken_krona_txt.dir,filename), 
               header = T, 
               sep = "\t",
               fill=TRUE,
               quote = "", 
               comment.char = "@",
               na.strings = c("","NA"))
  bracken_krona_table<-bracken_krona_table%>%filter(Abundance!=0)
  colnames(bracken_krona_table)[which(colnames(bracken_krona_table)=="Abundance")]<-sample.base.name
  
  list.tables[[sample.base.name]]<-bracken_krona_table
  
}

red.merged<-list.tables %>%
  Reduce(function(dtf1,dtf2) left_join(dtf1,
                                       dtf2,
                                       by=c("Kingdom","Phylum",
                                            "Class","Order",
                                            "Family","Genus","Species")), .)
df.otus<-red.merged%>%
  unite("OTU",Kingdom:Species,sep = "|")%>%
  column_to_rownames("OTU")
df.taxa<-red.merged%>%
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

# Create metadata
metadata.filename<-paste0(metadatadir,
                          paste("filenames-single-pooled-raw-supercomp.tsv", sep = "-"))
custom.md<-read.table(metadata.filename, header = T)
colnames(custom.md)[1]<-"Sample" # set the first column name as Sample
custom.md<-custom.md%>%
  filter(Sample%in%colnames(red.merged))%>%
  column_to_rownames(var = "Sample")%>%
  select(-absolute.filepath)
custom.md$class<-as.factor(custom.md$class)
custom.md$animal<-as.factor(custom.md$animal)
custom.md$sex<-as.factor(custom.md$sex)
custom.md$birthday<-as.Date(custom.md$birthday)
custom.md$age<-year(as.period(interval(custom.md$birthday,now())))

# Create a phyloseq object
ps.q<-phyloseq(otu_table(df.otus,taxa_are_rows = TRUE),
               tax_table(df.taxa),
               sample_data(custom.md))

if (agglom.rank=="Species"|agglom.rank=="OTU"){
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
  summarise(sumRelativeAbundance = sum(RelativeAbundance)) %>% # Sum all abundances
  pull(sumRelativeAbundance) %>% # get only a vector of numbers
  `==`(100) %>% # compare each value to 100 (TRUE/FALSE)
  all() # get one TRUE/FALSE value

ps.q.agg %>%
  mutate(sumRelativeAbundance = sum(RelativeAbundance)) %>%
  ggplot(aes(x=Sample,y=sumRelativeAbundance))+
  geom_bar(stat="identity")

# Group the dataframe by classes (animal hosts)
classcol<-which(colnames(ps.q.agg) =="class")
# find a column by which we agglomerated the dataset
# We find the column index because we can derive the previous taxonomic 
# rank with a simple agglom.rank.col-1.
if(agglom.rank=="Species"|agglom.rank=="OTU"){
  agglom.rank.col<-which(colnames(ps.q.agg) =="Species")
}else{
  agglom.rank.col<-which(colnames(ps.q.agg) ==agglom.rank)
}

# for each class, we take a agglom.rank, sum its abundances from all samples,
# then take a mean. This will be our MeanRelativeAbundance
if(agglom.rank=="Species"|agglom.rank=="OTU"){
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

objects.to.keep<-c("agglom.rank","ps.q.agg","custom.md")
objects.to.keep<-which(ls()%in%objects.to.keep)
rm(list = ls()[-objects.to.keep])
# Save the workspace
save.image(paste0("./output/rdafiles/",paste("kraken2",
                                             agglom.rank,
                                             "phyloseq-workspace.RData",sep = "-")))

sessionInfo()