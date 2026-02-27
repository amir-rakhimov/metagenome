# Combine all the statistics from MAG assembly ####
library(tidyverse)
# TODO: most abundant taxa
# TODO: abundance of eukaryotic contigs
#TODO: abundance of MAGs with low GC%
# TODO: treponema MAGs

# Tables ####
# 1. Summary statistics from assembly (**Table 1**)
# 2. Most abundant taxa ?


# Figures: ####
# 1. Taxonomy barplot (**My Figure 1**) 
# 2. Abundance of eukaryotic contigs (**My Figure 2**) 
# 3. Number of unique taxa with low GC% and their abundance (**My Figure 3**)

# 1. Distribution of MAGs in the data using coverage  (**My Figure ?**) 



load(file="./output/rdafiles/mag-combine-all-data.Rdata")
#' Combine BLASTN and GTDB-tk results (on all contigs from MEGAHIT) 
#' with seqkit stats (on contigs and dereplicated bins)

# Stats
#' Seqkit stats on MEGAHIT data
seqkit.megahit.date_time<-"20250712_18_07_37"
#' Seqkit stats on dereplicated MAGs
seqkit.drep.date_time<-"20250712_18_07_37"
#' Seqkit stats on MetaBAT2 data
seqkit.metabat2.date_time<-"20250712_18_07_37"

# Taxonomy data
#' BLASTN output
blastn.date_time<-"20250701_05_24_10"
#' GTDB-tk output
gtdbtk.date_time<-"20250705_19_24_02"
# To clean up GTDB-tk genome names
metabat2.date_time<-"20250612_13_37_47"

# MAG and contig coverages
coverm.drep.date_time<-"20250803_05_07_53"
coverm.metabat2.date_time<-"20250803_05_07_53"
coverm.megahit.date_time<-"20250803_05_07_53"


# Directories with output
seqkit.output.dir<-"./output/mag_assembly/seqkit_output"
blastn.output.dir<-file.path("./output/mag_assembly/blastn_output",
                             blastn.date_time)
gtdbtk.output.dir<-file.path("output/mag_assembly/gtdbtk_output")
coverm.output.dir<-file.path("output/mag_assembly/coverm_output")

# Specify file names
seqkit.megahit.fname<-file.path(seqkit.output.dir,
                                    paste(seqkit.megahit.date_time,
                                          "megahit_combined_stats.tsv",sep="_"))
seqkit.drep.fname<-file.path(seqkit.output.dir,
                                 paste(seqkit.drep.date_time,
                                       "sa_95perc_combined_stats.tsv",sep="_"))
seqkit.metabat2.fname<-file.path(seqkit.output.dir,
                                    paste(seqkit.metabat2.date_time,
                                          "metabat2_combined_stats.tsv",sep="_"))

blastn.taxonomy.fname<-file.path(blastn.output.dir,
                                    paste(blastn.date_time,
                                          "blastn_all_samples.tsv",sep="_"))
gtdbtk.taxonomy.fname<-file.path(gtdbtk.output.dir,
                                    paste(gtdbtk.date_time,
                                          "sa_95perc_classification_combined.tsv",
                                          sep="_"))

coverm.drep.fname<-file.path(coverm.output.dir,
                            paste(coverm.drep.date_time,
                                  "all_samples_sa_95perc_coverm.tsv",
                                  sep="_"))
coverm.metabat2.fname<-file.path(coverm.output.dir,
                            paste(coverm.metabat2.date_time,
                                  "all_samples_bins_coverm.tsv",
                                  sep="_"))
coverm.megahit.fname<-file.path(coverm.output.dir,
                            paste(coverm.megahit.date_time,
                                  "all_samples_megahit_coverm.tsv",
                                  sep="_"))

image.formats<-c("tiff","png")
barplot.directory<-"./images/barplots/"
boxplot.directory<-"./images/boxplots/"


# Read the seqkit stats ####
# Seqkit on MEGAHIT [VERY LARGE]
seqkit.megahit.df<-data.table::fread(seqkit.megahit.fname,
                                     header = TRUE,
                              sep="\t")%>%
  as_tibble()

# Seqkit on MetaBAT2
seqkit.metabat2.df<-read.table(seqkit.metabat2.fname,header = TRUE,
                              sep="\t",comment.char = "")%>%
  as_tibble()

# Seqkit on dRep
seqkit.drep.df<-read.table(seqkit.drep.fname,header = TRUE,
                              sep="\t",comment.char = "")%>%
  as_tibble()


head(seqkit.megahit.df)
head(seqkit.metabat2.df)
head(seqkit.drep.df)
#' 3118148 contigs from MEGAHIT
nrow(seqkit.megahit.df)
#' 156182 contigs from MetaBAT2
nrow(seqkit.metabat2.df)
#' 11456 contigs from dRep
nrow(seqkit.drep.df)
summary(seqkit.megahit.df)
summary(seqkit.metabat2.df)
summary(seqkit.drep.df)


## Add drep clusters information ####
# We want to know which bins were assigned to which clusters, and which
# bins weren't because they had low quality.
# We need it to get a more precise estimate of mapped reads. 
# Read the data with all high-quality MAGs and their clusters
drep.clusters.fname<-paste0("./output/mag_assembly/drep_output/",
                            "20250627_19_56_42_sa_95perc/data_tables/Cdb.csv")
drep.clusters<-read.csv(drep.clusters.fname)%>%
  as_tibble()%>%
  mutate(genome=gsub(paste0(metabat2.date_time,"_"),"",genome),
         genome=gsub(".fa","",genome))%>%
  separate_wider_delim(genome,delim = "_bin.",names = c("sample","bin_id"))%>%
  mutate(bin_id=as.integer(bin_id))
# Read the data with representative genomes and their clusters
drep.winning_genomes.fname<-paste0("./output/mag_assembly/drep_output/",
                                   "20250627_19_56_42_sa_95perc/data_tables/Wdb.csv")
drep.winning_genomes<-read.csv(drep.winning_genomes.fname)%>%
  as_tibble()%>%
  mutate(genome=gsub(paste0(metabat2.date_time,"_"),"",genome),
         genome=gsub(".fa","",genome))%>%
  separate_wider_delim(genome,delim = "_bin.",names = c("sample","bin_id"))%>%
  mutate(bin_id=as.integer(bin_id))

# Combine drep MAGs with other checkm2 MAGs
drep.clusters.joined<-drep.winning_genomes%>%
  select(sample,bin_id,cluster)%>%
  rename("secondary_cluster" = "cluster",
         "drep_sample"= "sample",
         "drep_bin_id"= "bin_id")%>%
  full_join(drep.clusters[,c("sample","bin_id","secondary_cluster")],
            join_by("secondary_cluster"))%>%
  # Add a column that shows if a bin is representative or not.
  # Then, remove drep_sample and drep_bin_id because it's the same information
  mutate(is_drep_bin=ifelse(sample==drep_sample&bin_id==drep_bin_id,
                            TRUE,FALSE))%>%
  select(-drep_sample,-drep_bin_id)

rm(drep.clusters)
rm(drep.clusters.fname)
rm(drep.winning_genomes)
rm(drep.winning_genomes.fname)
### Join cluster information with seqkit stats ####
seqkit.drep.df<-seqkit.drep.df%>%
  left_join(drep.clusters.joined,by = join_by(sample, bin_id))

seqkit.metabat2.df<-seqkit.metabat2.df%>%
  left_join(drep.clusters.joined,by = join_by(sample, bin_id))
rm(drep.clusters.joined)

# Combine all dfs into one ####
seqkit.metabat2.df<-seqkit.metabat2.df%>%
  mutate(sample_contig=paste(sample,contig_id,sep = "_"),
         sample_bin_contig=paste(sample,bin_id,contig_id,sep = "_"))
seqkit.drep.df<-seqkit.drep.df%>%
  mutate(sample_contig=paste(sample,contig_id,sep = "_"),
         sample_bin_contig=paste(sample,bin_id,contig_id,sep = "_"))
seqkit.megahit.df<-seqkit.megahit.df%>%
  mutate(sample_contig=paste(sample,contig_id,sep = "_"))

# Get the common IDs between datasets to know whether a contig
# was assembled (metabat2 or metabat2_drep) or not (megahit).
# If it's metabat2, it is not representative. If it's metabat2_drep,
# it is representative.
megahit.metabat2.common.ids<-intersect(seqkit.megahit.df$sample_contig,
                                      seqkit.metabat2.df$sample_contig)
metabat2.drep.common.ids<-intersect(seqkit.metabat2.df$sample_bin_contig,
                                   seqkit.drep.df$sample_bin_contig)
megahit.drep.common.ids<-intersect(seqkit.megahit.df$sample_contig,
                                   seqkit.drep.df$sample_contig)

# All contigs from drep are found in metabat2
length(setdiff(seqkit.drep.df$sample_bin_contig,
        seqkit.metabat2.df$sample_bin_contig))
# All contigs from metabat2 are found in megahit
length(setdiff(seqkit.metabat2.df$sample_contig,
        seqkit.megahit.df$sample_contig))

## Combined df ####
seqkit.all.df<-seqkit.megahit.df%>%
  full_join(subset(seqkit.metabat2.df,select = -c(AT,GC,sample_contig,
                                                      contig_length)),
            by = join_by(sample, contig_id))%>%
  relocate(sample,contig_id,bin_id)%>%
  # If sample_contig is in the vector of IDs common between megahit and 
  # metabat2 IDs, it was used in the assembly by metabat2.
  # If sample_bin_contig is in the vector of IDs common between metabat2 and 
  # drep, it was selected as representative.
  mutate(source=ifelse(sample_contig%in%megahit.metabat2.common.ids,
                       "metabat2","megahit"),
         source=ifelse(sample_bin_contig%in%metabat2.drep.common.ids,
                       "metabat2_drep",source))%>%
  select(-sample_contig,-sample_bin_contig)

rm(megahit.drep.common.ids)
rm(megahit.metabat2.common.ids)
rm(metabat2.drep.common.ids)

# Barplot to count how many contigs were retained: most contigs weren't assembled
seqkit.all.df%>%
  count(sample,source)%>%
  pivot_wider(names_from = "source",
              values_from="n")

seqkit.all.df%>%
  ggplot(aes(x=sample,fill=source))+
  geom_bar(stat="count")


rm(seqkit.drep.df)
rm(seqkit.megahit.df)
rm(seqkit.metabat2.df)
gc()

# Add taxonomy ####
## Add GTDB-tk taxonomy (for dereplicated bins) ####
gtdbtk.taxonomy<-read.table(gtdbtk.taxonomy.fname,header = T,
                            sep="\t",comment.char = "")%>%
  as_tibble()%>%
  select(user_genome,classification) # remove unnecessary columns
head(gtdbtk.taxonomy)
# 128 bins
nrow(gtdbtk.taxonomy)

# Tidy up GTDB-tk taxonomy: substitute s__ and g__ with unclassified
# Separate sample_bin column into separate columns
gtdbtk.taxonomy<-gtdbtk.taxonomy %>%
  mutate(classification=gsub(";g__;s__",";unclassified;unclassified",
                             classification))%>%
  mutate(classification=gsub(";s__$",";unclassified",classification))%>%
  mutate(classification=gsub("d__|p__|c__|o__|f__|g__|s__","",classification)) %>%
  mutate(user_genome=gsub(paste0(metabat2.date_time,"_"),"",user_genome))%>%
  separate_wider_delim(user_genome,delim="_bin.",names=c("sample","bin_id"))%>%
  mutate(bin_id=as.integer(bin_id))

## Add BLAST taxonomy (for all contigs from MEGAHIT) ####
blastn.taxonomy<-read.table(blastn.taxonomy.fname,header = F,
                            sep="\t",comment.char = "")%>%
  as_tibble()
head(blastn.taxonomy)
# 1899 contigs
nrow(blastn.taxonomy)
colnames(blastn.taxonomy)<-c("sample","contig_id","silva_ref_id","perc_identity",
                             "alignment_length","bitscore","classification")

# Join stats with taxonomy ####
seqkit.with_taxonomy<-seqkit.all.df%>%
  full_join(gtdbtk.taxonomy,
            by = join_by(sample, bin_id))%>%
  full_join(blastn.taxonomy[,c("sample","contig_id","classification")],
            by = join_by(sample, contig_id))%>%
  rename("gtdbtk_result"="classification.x",
         "blastn_result"="classification.y")%>%
  ## Add column to show if a contig was classified, and which pipeline classified it ####
  mutate(tax.classification.type=ifelse(!is.na(gtdbtk_result)&!is.na(blastn_result),
                        "gtdbtk_and_blastn","unclassified"),
         tax.classification.type=ifelse(!is.na(gtdbtk_result)&is.na(blastn_result),
                        "gtdbtk_only",tax.classification.type),
         tax.classification.type=ifelse(is.na(gtdbtk_result)&!is.na(blastn_result),
                        "blastn_only",tax.classification.type))
rm(seqkit.all.df)


# Add CoverM data  ####
## CoverM data from dRep ####
coverm.drep.df<-read.table(coverm.drep.fname,header = T,
                           sep = "\t")%>%
  as_tibble()%>%
  mutate(genome=gsub(paste0(metabat2.date_time,"_"),"",genome))

## CoverM data from MetaBAT2 ####
coverm.metabat2.df<-read.table(coverm.metabat2.fname,header = T,
                               sep = "\t")%>%
  as_tibble()%>%
  mutate(genome=gsub(paste0(metabat2.date_time,"_"),"",genome))

## CoverM data from MEGAHIT (use data.table package)####
coverm.megahit.df<-data.table::fread(coverm.megahit.fname,header = T,
                                  sep = "\t")%>%
  as_tibble()%>%
  rename(sample_contig=contig)%>%
  separate_wider_delim(sample_contig,delim = "_",
                       names = c("sample","contig_id"),
                       too_many="merge" )

# How much relative abundance is classified in drep bins
coverm.drep.df%>%
  mutate(sample=ifelse(grepl("unmapped",genome),
                       gsub("unmapped_","",genome),
                       gsub("_bin.[0-9]*","",genome)))%>%
  filter(!grepl("unmapped",genome))%>%
  group_by(sample)%>%
  summarise(sum_relative_ab=sum(relative_abundance))

# How much relative abundance is classified in metabat2 bins
coverm.metabat2.df%>%
  mutate(sample=ifelse(grepl("unmapped",genome),
                       gsub("unmapped_","",genome),
                       gsub("_bin.[0-9]*","",genome)))%>%
  filter(!grepl("unmapped",genome))%>%
  group_by(sample)%>%
  summarise(sum_relative_ab=sum(relative_abundance))


# How much relative abundance is classified in metabat2 bins that belonged to
# drep clusters
seqkit.with_taxonomy%>%
  mutate(genome=paste(sample,bin_id,sep="_bin."))%>%
  distinct(sample,bin_id,genome,secondary_cluster)%>%
  left_join(coverm.metabat2.df,by=join_by("genome"))%>%
  filter(!is.na(secondary_cluster))%>%
  group_by(sample)%>%
  summarise(sum_relative_ab=sum(relative_abundance))

# Save image ####
save.image(file="./output/rdafiles/mag-combine-all-data.Rdata")

# Join relative abundance with taxonomy ####
## GTDB-tk + dRep ####
gtdbtk.coverm.drep<-seqkit.with_taxonomy%>%
  filter(!is.na(secondary_cluster))%>%
  select(sample,bin_id,secondary_cluster)%>%
  distinct()%>%
  right_join(gtdbtk.taxonomy)%>%
  mutate(genome=paste(sample,bin_id,sep="_bin."))%>%
  left_join(coverm.drep.df[,c("genome","relative_abundance")],
            by=join_by("genome"))%>%
  select(-genome)

## GTDB-tk + MetaBAT2 ####
gtdbtk.coverm.metabat2<-seqkit.with_taxonomy%>%
  distinct(sample,bin_id,secondary_cluster)%>%
  filter(!is.na(secondary_cluster))%>%
  full_join(gtdbtk.taxonomy,by = join_by(sample, bin_id))%>%
  group_by(secondary_cluster) %>%
  mutate(classification = classification[!is.na(classification)][1L])%>%
  ungroup%>%
  mutate(genome=paste(sample,bin_id,sep="_bin."))%>%
  left_join(coverm.metabat2.df[,c("genome","relative_abundance")],
            by=join_by("genome"))%>%
  select(-genome)

## BLASTN + MEGAHIT
blastn.coverm.megahit<-coverm.megahit.df%>%
  select(sample, contig_id,tpm)%>%
  inner_join(blastn.taxonomy,by = join_by(sample, contig_id))

rm(coverm.drep.df)
rm(coverm.metabat2.df)
rm(coverm.megahit.df)

# 
# gtdbtk.taxonomy%>%
#   distinct(sample,bin_id,classification)%>%
#   left_join(drep.winning_genomes[,c("sample","bin_id","cluster")])%>%
#   full_join(seqkit.metabat2.df,by=join_by(sample,bin_id))%>%
#   distinct(sample,bin_id,.keep_all = T)%>%View

# Most contigs weren't classified
# Actually, they probably were. We don't see it because we classified
# only dereplicated genomes. So, redundant contigs were ignored.
seqkit.with_taxonomy%>%
  ggplot(aes(x=sample,fill=tax.classification.type))+
  geom_bar(stat="count")


# Now, fill NAs for bins that were not used in GTDBtk but actually 
# belonged to a drep cluster. They were redundant but their cluster was 
# classified, so they were supposed to be classified, too.

# Fill NA by group (secondary_cluster). Group by secondary_cluster, 
# then select the first non-NA value [1L] and use that non-NA value to 
# fill NAs for other entries in the group 
seqkit.with_taxonomy%>%
  group_by(sample,bin_id,secondary_cluster) %>% 
  mutate(gtdbtk_result = gtdbtk_result[!is.na(gtdbtk_result)][1L])%>%
  ungroup%>%
  mutate(tax.classification.type=ifelse(!is.na(gtdbtk_result)&is.na(blastn_result),
                                        "gtdbtk_only",tax.classification.type))%>%
  # filter(tax.classification.type!="unclassified")%>%
  ggplot(aes(x=sample,fill=tax.classification.type))+
  geom_bar(stat="count")


# But among classified contigs, GTDB-tk classified more
seqkit.with_taxonomy%>%
  filter(tax.classification.type!="unclassified")%>%
  ggplot(aes(x=sample,fill=tax.classification.type))+
  geom_bar(stat="count")



# Now remove unclassified contigs to get a better view. Now there's more contigs
seqkit.with_taxonomy%>%
  group_by(sample,bin_id,secondary_cluster) %>% 
  mutate(gtdbtk_result = gtdbtk_result[!is.na(gtdbtk_result)][1L])%>%
  ungroup%>%
  mutate(tax.classification.type=ifelse(!is.na(gtdbtk_result)&is.na(blastn_result),
                                        "gtdbtk_only",tax.classification.type))%>%
  filter(tax.classification.type!="unclassified")%>%
  ggplot(aes(x=sample,fill=tax.classification.type))+
  geom_bar(stat="count")



# Save data ####
saveRDS(seqkit.with_taxonomy,file = "./output/rdafiles/seqkit-with_taxonomy.rds")
saveRDS(gtdbtk.coverm.drep,file = "./output/rdafiles/gtdbtk-coverm-drep.rds")
saveRDS(gtdbtk.coverm.metabat2,file = "./output/rdafiles/gtdbtk-coverm-metabat2.rds")
saveRDS(blastn.coverm.megahit,file = "./output/rdafiles/blastn-coverm-megahit.rds")

data.table::fwrite(seqkit.with_taxonomy,
            file = "./output/rtables/seqkit-with_taxonomy.tsv",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE,
            sep = '\t',
            na =NA)

write.table(gtdbtk.coverm.drep,
            file = "./output/rtables/gtdbtk-coverm-drep.tsv",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE,
            sep = '\t')

write.table(gtdbtk.coverm.metabat2,
            file = "./output/rtables/gtdbtk-coverm-metabat2.tsv",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE,
            sep = '\t')

write.table(blastn.coverm.megahit,
            file = "./output/rtables/blastn-coverm-megahit.tsv",
            col.names = TRUE,
            row.names = FALSE,
            quote = FALSE,
            sep = '\t')
