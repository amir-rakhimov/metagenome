# install.packages(c("tidyverse"))
# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# BiocManager::install("phyloseq")
# install.packages(
#   "microViz",
#   repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
# )
library(tidyverse)
library(microViz)
library(phyloseq)
rdafiles.directory<-"./output/rdafiles"
rtables.directory<-"./output/rtables"

agglom.rank<-"Species"
metadatadir<-paste0("../amplicon_nmr/data/metadata/pooled-metadata/") # directory with metadata
pipeline.name<-"kraken2"
pipeline.name<-"singlem"
if(pipeline.name=="kraken2"){
  # Kraken2
  combined.report.date_time<-"20240515_10_04_04"
  combined.report.filename<-file.path(rtables.directory,
                                      paste(combined.report.date_time,"combined_report.tsv",sep = "_"))
  # combined.report.filename<-file.path("output/kraken2_pipeline/kraken2_reports",
  #                                     paste(combined.report.date_time,"combined_mpa_clean.tsv",sep = "_"))
}else if (pipeline.name=="singlem"){
  # Singlem
  combined.report.date_time<-"20240929_23_31_55"
  combined.report.filename<-file.path(rtables.directory,
                                      paste(combined.report.date_time,"singlem_combined_profile.tsv",sep = "_"))
}

combined.report<-read.table(combined.report.filename, 
            header = T, 
            sep = "\t",
            fill=TRUE,
            quote = "", 
            comment.char = "@",
            na.strings = c("","NA"))
colnames(combined.report)<-gsub("X2D10","2D10",colnames(combined.report))
colnames(combined.report)<-gsub("X2D14","2D14",colnames(combined.report))

dim(combined.report)
head(combined.report)

# all.ranks<-c("Kingdom", "Phylum", "Class", "Order", "Family","Genus","Species")
# combined.report[,all.ranks][is.na(combined.report[,all.ranks])]<-"NA"
# combined.report[,-which(names(combined.report)%in%all.ranks)][is.na(combined.report[,-which(names(combined.report)%in%all.ranks)])]<-0

df.otus<-combined.report%>%
  unite("OTU",Kingdom:Species,sep = "|")%>%
  column_to_rownames("OTU")

df.otus[is.na(df.otus)]<-0

df.taxa<-combined.report%>%
  select(Kingdom:Species)%>%
  unite("OTU",Kingdom:Species,sep = "|",remove = FALSE)%>%
  column_to_rownames("OTU")

df.taxa$Kingdom<-
  gsub("k__","",df.taxa$Kingdom)
df.taxa$Kingdom<-
  gsub("d__","",df.taxa$Kingdom)
df.taxa$Phylum<-
  gsub("d__|p__","",df.taxa$Phylum)
df.taxa$Class<-
  gsub("d__|p__|c__","",df.taxa$Class)
df.taxa$Order<-
  gsub("d__|p__|c__|o__","",df.taxa$Order)
df.taxa$Family<-
  gsub("d__|p__|c__|o__|f__","",df.taxa$Family)
df.taxa$Genus<-
  gsub("d__|p__|c__|o__|f__|g__","",df.taxa$Genus)
df.taxa$Species<-
  gsub("d__|p__|c__|o__|f__|g__|s__","",df.taxa$Species)

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
  filter(Sample%in%colnames(combined.report))%>%
  column_to_rownames(var = "Sample")%>%
  select(-absolute.filepath)
custom.md$class<-as.factor(custom.md$class)
custom.md$animal<-as.factor(custom.md$animal)
custom.md$sex<-as.factor(custom.md$sex)
custom.md$birthday<-as.Date(custom.md$birthday)
# calculate age at the time when data was received
custom.md$age<-year(as.period(interval(custom.md$birthday,as.Date("2023-11-16"))))

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
# Remove human data
ps.q.agg<-ps.q.agg%>%
  filter(!get(agglom.rank)%in%c("Chordata","Mammalia","Primates",
                                "Hominidae","Homo","Homo_sapiens"))
  
# Remove entries with zero Abundance
ps.q.agg<-ps.q.agg%>%
  filter(Abundance!=0)
tax.levels<-c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
tax.levels<-tax.levels[1:which(agglom.rank==tax.levels)]

ps.q.agg<-ps.q.agg%>%
  select(all_of(c("OTU","Sample","class","Abundance",tax.levels)))

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
  group_by(Sample)%>%
  summarise(sumRelativeAbundance = sum(RelativeAbundance)) %>%
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
    group_by_at(c("class",agglom.rank))%>%
    mutate(TotalAgglomRank=sum(Abundance))%>%
    mutate(MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)
}else{ 
  ps.q.agg<-ps.q.agg%>%
    group_by(class)%>% # group by class,
    # compute MeanRelativeAbundance from Abundance 
    mutate(TotalClass=sum(Abundance))%>%
    group_by_at(c("class",agglom.rank))%>%
    mutate(TotalAgglomRank=sum(Abundance))%>%
    mutate(MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)%>%
    select(-OTU)
}

ps.q.agg<-ps.q.agg%>%
  select(-TotalClass,-TotalSample)%>%
  ungroup()

# Save the ASV table
write.table(ps.q.agg,
        file=file.path(rtables.directory,paste(
          paste(format(Sys.time(),format="%Y%m%d"),
                format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
          "phyloseq",pipeline.name,agglom.rank,
          "table.tsv",sep="-")),
        row.names = F,sep = "\t")
saveRDS(ps.q.agg,
            file=file.path(rdafiles.directory,paste(
              paste(format(Sys.time(),format="%Y%m%d"),
                    format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
              "phyloseq",pipeline.name,agglom.rank,
              "table.rds",sep="-")))

sessionInfo()
