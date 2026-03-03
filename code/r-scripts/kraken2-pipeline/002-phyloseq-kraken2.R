#' ---
#' output: 
#'   bookdown::html_document2:
#'      toc: true
#' ---
#' 
#' ```{r, setup 002-phyloseq-kraken2.R, include=FALSE}
#' knitr::opts_knit$set(root.dir = '/home/rakhimov/projects/metagenome')
#' ```
#' 
#+ echo=FALSE
# Convert Kraken2 output into a phyloseq object ####
#' # Convert Kraken2 output into a phyloseq object
#' 
#+ echo=FALSE
## Introduction ####
#'
#' ## Introduction
#' In this script, we will import the dataset from Kraken2 using 
#' phyloseq.
#' 

#+ echo=FALSE
## 1. Load necessary libraries. ####
#'
#' ## Load necessary libraries
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
#' Directories with input files:
rdafiles.directory<-"./output/rdafiles"
rtables.directory<-"./output/rtables"

#' Directory with metadata:
# metadatadir<-paste0("../amplicon_nmr/data/metadata/pooled-metadata/") 
metadatadir<-paste0("../amplicon_nmr/output/rdafiles") 

#+ echo=FALSE
## 2. Import the combined report from Kraken2. #### 
#'
#' ## Import the combined report from Kraken2.
combined.report.date_time<-"20240515_10_04_04"
combined.report.filename<-file.path(rtables.directory,
                                    paste(combined.report.date_time,"combined_report.tsv",sep = "_"))
combined.report<-read.table(combined.report.filename, 
            header = T, 
            sep = "\t",
            fill=TRUE,
            quote = "", 
            comment.char = "@",
            na.strings = c("","NA"))%>%
  as_tibble()

#' Convert into wide format:
combined.report<-combined.report%>%
  filter(Abundance!=0)%>%
  pivot_wider(names_from = Sample,
              values_from = Abundance)

dim(combined.report)
head(combined.report)

# all.ranks<-c("Kingdom", "Phylum", "Class", "Order", "Family","Genus","Species")
# combined.report[,all.ranks][is.na(combined.report[,all.ranks])]<-"NA"
# combined.report[,-which(names(combined.report)%in%all.ranks)][is.na(combined.report[,-which(names(combined.report)%in%all.ranks)])]<-0

#' Extract table with Species abundances where rows are Species and columns 
#' are samples. Here, we combine all the columns into one and call it "OTU".
df.otus<-combined.report%>%
  unite("OTU",Kingdom:Species,sep = "|")%>%
  column_to_rownames("OTU")

df.otus[is.na(df.otus)]<-0

#' Extract table with taxonomy where rows are OTUS and columns are taxonomy ranks:
#' Kingdom, ..., Species
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
head(df.taxa)
#+ echo=FALSE
## 3. Convert into phyloseq taxonomyTable object. ####
#'
#' ## Convert into phyloseq taxonomyTable object.
df.taxa<-tax_table(as.matrix(df.taxa))

df.taxa.pretty<-tax_fix(df.taxa)
df.taxa<-df.taxa.pretty
rm(df.taxa.pretty)

#' Add metadata
metadata.filename<-file.path(metadatadir,"custom.md.rds")
custom.md<-readRDS(metadata.filename)%>%
  filter(class =="NMR")%>%
  dplyr::select(Sample,class)
custom.md<-custom.md%>%
  filter(Sample%in%colnames(combined.report))

#+ echo=FALSE
## 4. Create a phyloseq object. ####
#'
#' ## Create a phyloseq object. 
ps.q<-phyloseq(otu_table(df.otus,taxa_are_rows = TRUE),
               tax_table(df.taxa),
               sample_data(custom.md))
#+ echo=FALSE
## 5. Convert the phyloseq object into a dataframe. ####
#'
#' ## Convert the phyloseq object into a dataframe.
ps.q.agg<-ps.q %>%
  psmelt() 
ps.q.agg.phylum<-ps.q %>%
  tax_glom("Phylum",NArm = FALSE) %>% # agglomerate by phylum
  psmelt()  # transform the phyloseq object into an R dataframe
ps.q.agg.family<-ps.q %>%
  tax_glom("Family",NArm = FALSE) %>% # agglomerate by family
  psmelt()  # transform the phyloseq object into an R dataframe
ps.q.agg.genus<-ps.q %>%
  tax_glom("Genus",NArm = FALSE) %>% # agglomerate by genus
  psmelt()  # transform the phyloseq object into an R dataframe
ps.list <- list("OTU" = ps.q.agg, 
                "Phylum" = ps.q.agg.phylum,
                "Family" = ps.q.agg.family, 
                "Genus" = ps.q.agg.genus)

for (ps.df.index in names(ps.list)){
  ps.df <- ps.list[[ps.df.index]]
  print(paste("Parsing data from the", ps.df.index, "table"))
  print(paste("Number of rows in the dataframe:", nrow(ps.df)))
  
  #' Remove human data:
  print("Removing rows with human data")
  ps.df<-ps.df%>%
    filter(!grepl(paste(c("Chordata","Mammalia","Primates",
                        "Hominidae","Homo","Homo_sapiens"),collapse = "|"),
                  get(ps.df.index)))
  # Remove entries with zero Abundance.
  print("Removing rows with zero Abundance")
  ps.df <- ps.df %>%
    dplyr::select(-sample_Sample)%>% # remove the duplicate column
    filter(Abundance!=0)
  print(paste("Number of rows in the filtered dataset:",nrow(ps.df)))
  ### 5.1 Number of samples in the filtered dataset. ####
  print(paste("Number of samples in the filtered dataset:"))
  ps.df%>%
    distinct(Sample)%>%
    tally()%>%
    print()
  ### 5.2 Number of features in the filtered dataset. ####
  print(paste("Number of features in the filtered dataset:"))
  ps.df%>%
    distinct(OTU)%>%
    tally()%>%
    print()
  ### 5.3 Total frequency in the filtered dataset. ####
  print(paste("Total frequency in the filtered dataset:"))
  ps.df%>%
    summarise(TotalAbundance=sum(Abundance))%>%
    print()
  ### 5.4 Summary statistics (min, median, max, quartiles) of the filtered dataset. ####
  print(paste("Summary statistics (min, median, max, quartiles) of the filtered dataset:"))
  ps.df%>%
    dplyr::select(Sample,Abundance)%>%
    group_by(Sample)%>%
    summarise(FrequencyPerSample=sum(Abundance))%>%
    dplyr::select(FrequencyPerSample)%>%
    summary()%>%
    print()
  ## 6. Add relative abundance column: Abundance divided by total abundance in a sample. ####
  ps.df<-ps.df%>%
    group_by(class,Sample)%>%
    mutate(TotalSample=sum(Abundance))%>%
    group_by_at(c("class","Sample",ps.df.index))%>%
    mutate(RelativeAbundance=Abundance/TotalSample*100)%>%
    ungroup()
  # Sanity check: is total relative abundance of each sample 100%?
  print(paste("Sanity check: is total relative abundance of each sample 100%?"))
  ps.df %>%
    group_by(Sample) %>% # Group by sample id
    summarise(sumRelativeAbundance = sum(RelativeAbundance)) %>% # Sum all abundances
    mutate(diff_from_100 = sumRelativeAbundance-100, # compare each value to 100
           is_different = as.logical(round(diff_from_100,digits = 10)))%>% 
    arrange(desc(diff_from_100))%>% # show the most deviating samples
    print()
  ps.df %>%
    group_by(Sample)%>%
    mutate(sumRelativeAbundance = sum(RelativeAbundance)) %>%
    ungroup()%>%
    distinct(Sample,.keep_all = T)%>%
    ggplot(aes(x=Sample,y=sumRelativeAbundance))+
    geom_bar(stat="identity")+
    coord_flip()
  print(last_plot())
  ### 6.1 Add mean relative abundance data. ####
  # We will group the dataset by three columns: class (animal host), 
  # two taxonomic ranks (e.g Genus, Family), and maybe OTU (actually ASV)
  # if we agglomerate at ASV level.
  
  # Group the dataframe by classes (animal hosts).
  # First, we calculate the library size per sample.
  # Then, inside each class, we take a agglom.rank, sum its abundances from all samples,
  # then take a mean. This will be our MeanRelativeAbundance.
  ps.df<-ps.df%>%
    group_by(class)%>% # group by class (animal host),
    mutate(TotalClass=sum(Abundance))%>%
    group_by_at(c("class",ps.df.index))%>%
    mutate(TotalAgglomRank=sum(Abundance))%>%
    mutate(MeanRelativeAbundance=TotalAgglomRank/TotalClass*100)%>%
    ungroup()
  if(ps.df.index != "OTU"){
    ps.df<-ps.df%>%
      dplyr::select(-OTU)
  }
  ps.df<-ps.df%>%
    ungroup()%>%
    dplyr::select(-TotalClass,-TotalSample,-TotalAgglomRank)
  ps.list[[ps.df.index]]<-ps.df
  
  # Save the tables in TSV format and as an RDS object
  # write.table(ps.q.agg,
  #             file=file.path(rtables.directory,paste(
  #               paste(format(Sys.time(),format="%Y%m%d"),
  #                     format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
  #               "phyloseq-kraken2",agglom.rank,
  #               "table.tsv",sep="-")),
  #             row.names = F,sep = "\t")
  # saveRDS(ps.q.agg,
  #         file=file.path(rdafiles.directory,paste(
  #           paste(format(Sys.time(),format="%Y%m%d"),
  #                 format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
  #           "phyloseq-kraken2",agglom.rank,
  #           "table.rds",sep="-")))
  
}

sessionInfo()
rm(list = ls(all=TRUE))
gc()