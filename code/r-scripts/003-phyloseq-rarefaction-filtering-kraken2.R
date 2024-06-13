library(tidyverse)
library(phyloseq)
library(vegan)
date_time<-"20240518_13_40_08"
agglom.rank<-"Species"
# Load the Workspace from phyloseq (output of 001-phyloseq-qiime2.R)
load(file.path("./output/rdafiles",paste(
  date_time,"kraken2",agglom.rank,
  "phyloseq-workspace.RData",sep = "-")))

# custom.levels<-c("NMR","B6mouse")
# better labels for facets
pretty.facet.labels<-c("NMR" = "*Heterocephalus glaber*"
)
custom.levels<-intersect(names(pretty.facet.labels),custom.md$class)

# If we're working with NMR only
custom.levels<-c("NMR")

# Import data ####
# For some reason, one family is duplicated
ps.q.df <-ps.q.agg%>%
  select(all_of(c("Sample","Abundance","class",agglom.rank)))%>%
  group_by_at(c("class","Sample",agglom.rank))%>%
  summarise(Abundance=sum(Abundance))%>%
  filter(Abundance!=0)%>%
  ungroup()

ps.q.df<-ps.q.df%>%
  filter(class %in% custom.levels,Abundance!=0)

# convert the data frame into wide format 
ps.q.df.wide<-ps.q.df%>%
  pivot_wider(names_from = agglom.rank, 
              values_from = "Abundance",
              values_fill = 0)%>%
  as.data.frame()

rownames(ps.q.df.wide)<-ps.q.df.wide$Sample
ps.q.df.wide<-ps.q.df.wide[,-c(1,2)]  

# Rarefaction ####
# find the smallest sample size
min.n_seqs.all<-ps.q.df%>%
  group_by(Sample)%>%
  summarize(n_seqs=sum(Abundance))%>%
  summarize(min=min(n_seqs))%>%
  pull(min)

# rarefied asv table with vegan
set.seed(1)
ps.q.df.rare<-rrarefy(ps.q.df.wide,sample=min.n_seqs.all)
ps.q.df.rare<-ps.q.df.rare%>%
  as_tibble(rownames="Sample")%>%
  pivot_longer(-Sample)%>%
  as.data.frame()%>%
  left_join(unique(ps.q.agg[,c("Sample","class","sex","birthday","age")]),
             by="Sample")
# rename the name column
# !!agglom.rank:= will evaluate the variable
ps.q.df.rare<-ps.q.df.rare%>%
  rename(!!agglom.rank:=name,Abundance=value)%>%
  filter(Abundance!=0)
write.table(ps.q.df.rare,
            file = file.path("./output/rtables",paste0(
              paste(
                paste(format(Sys.time(),format="%Y%m%d"),
                      format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                "kraken2-ps.q.df.rare-nonfiltered",agglom.rank,
                paste(custom.levels,collapse = '-'),sep = "-"),
              ".tsv")),
            row.names = F,
            sep = "\t")

# Filtering by prevalence ####
# for each host, calculate the number of samples an ASV was observed in 
# so the same ASV can be observed in multiple hosts.
# observed_samples shows number of samples per host
ps.q.sample_counts <- ps.q.agg%>%
  select(OTU,Sample,class,Abundance,Taxon,Genus,Family)%>%
  filter(Abundance!=0)%>%
  group_by(class,Family,Genus) %>% 
  summarize(observed_samples = n_distinct(Sample))

# calculate the total number of samples for each host
ps.q.host_counts <- ps.q.agg%>%
  select(OTU,Sample,class,Abundance,Taxon,Genus,Family)%>%
  filter(Abundance!=0)%>%
  group_by(class)%>%
  summarize(total_samples = n_distinct(Sample)) %>%
  ungroup()

# Join host counts with sample counts
# Create the PercentageSamples column: num of samples  from a specific host
# an ASV was observed in divided by the total number of samples in that host
# Then, we join our dataframe with mean relative abundance df (ps.agg.rel.pooled).
# Then, we add the Taxon column from ps.q.agg.abs
ps.q.prevalences <- ps.q.host_counts %>%
  left_join(ps.q.sample_counts, by = "class") %>%
  mutate(PercentageSamples = observed_samples/total_samples*100)%>% # percentage
  # of samples an ASV was observed in
  left_join(ps.q.agg[,c("class","Genus","Family","Taxon","MeanRelativeAbundance")],
            by=c("class","Genus","Family")) # add mean
  # relative abundances and the Taxon column

## filter by percentage of samples an ASV was observed in ####
ps.q.df.norare.filtered <-ps.q.df%>%
  left_join(ps.q.prevalences,by=c("class","Taxon"),relationship = "many-to-many")%>%
  filter(PercentageSamples>=10)
ps.q.df.rare.filtered<-ps.q.df.rare%>%
  left_join(ps.q.prevalences,by=c("class","Taxon"),relationship = "many-to-many")%>%
  filter(PercentageSamples>=10,MeanRelativeAbundance>=1)


write.table(ps.q.df.norare.filtered,
            file = paste0("./rtables/",authorname,"/ps.q.df.norare.filtered-",
                          paste(agglom.rank,custom.levels,collapse = '-'),".tsv"),
            row.names = F,
            sep = "\t")
write.table(ps.q.df.rare.filtered,
            file = paste0("./rtables/",authorname,"/ps.q.df.rare.filtered-",
                          paste(agglom.rank,custom.levels,collapse = '-'),".tsv"),
            row.names = F,
            sep = "\t")
