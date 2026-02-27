library(tidyverse)
library(phyloseq)
library(Maaslin2)
library(vegan)
# phyloseq.date_time<-"20240619_14_11_06"
pipeline.name<-"kraken2"
ps.q.df.preprocessed.date_time<-"20241004_15_12_22"

rdafiles.directory<-"./output/rdafiles"
rtables.directory<-"./output/rtables"
metadata.directory<-"../amplicon_nmr/output/rdafiles"
# choose what to compare
comparison<-"age"
# comparison<-"sex"
# comparison<-"strain"
# choose the host of interest
host<-"NMR"
# host<-"mice"
ref.level<-"agegroup0_10" # choose the reference level
# ref.level<-"F"
# this is for file names
host.labels<-c("NMR" = "*Heterocephalus glaber*")
custom.levels<-c("NMR")
agglom.rank<-"Species"
custom.md<-readRDS(file.path(metadata.directory,"custom.md.ages.rds"))%>%
  filter(sequencing_type == "Naked mole-rat whole metagenome sequencing")

rare.status<-"rare"
filter.status<-"nonfiltered"

# Import data ####
ps.q.df.preprocessed<-read.table(
  file.path(rtables.directory,paste0(
    paste(ps.q.df.preprocessed.date_time,
          pipeline.name,"ps.q.df.rare-nonfiltered",agglom.rank,
          paste(host,collapse = '-'),sep = "-"),".tsv")),
  header = T,sep = "\t")

ps.q.df.preprocessed$Sample<-as.factor(ps.q.df.preprocessed$Sample)
ps.q.df.preprocessed$class<-as.factor(ps.q.df.preprocessed$class)

# add age groups
ps.q.df.preprocessed <- ps.q.df.preprocessed%>%
  left_join(custom.md[,c("Sample", "agegroup")])

# we create these new levels because maaslin is itsy bitsy
unique_levels<-custom.md%>%
  distinct(agegroup)%>%
  arrange(agegroup)%>%
  pull()

# Creating custom levels ####
if (comparison=="age"){
  # names for levels are age groups
  pretty.level.names<-names(table(custom.md$old_agegroup))
  names(pretty.level.names)<-names(table(custom.md$agegroup))
  custom.levels<-names(pretty.level.names)
  
}else if (comparison=="sex"){
  pretty.level.names<-
    c("F" = "Females",
      "M" = "Males")
  custom.levels<-names(pretty.level.names)
  
}else if(comparison=="strain"){
  pretty.level.names<-
    c("B6mouse" = "B6 mouse",
      "MSMmouse" = "MSM/Ms mouse",
      "FVBNmouse" = "FVB/N mouse"
    )
  custom.levels<-intersect(names(pretty.level.names),custom.md$class)
  pretty.level.names<-pretty.level.names[which(names(pretty.level.names)%in%custom.levels)]
}
# Preparing the dataset ####
# filter the dataset
ps.q.df <-ps.q.df.preprocessed%>%
  dplyr::select(all_of(c("Sample","Abundance","class",agglom.rank)))%>%
  filter(Abundance!=0)
ps.q.df.wide<-ps.q.df%>%
  pivot_wider(names_from = agglom.rank, # or OTU
              values_from = "Abundance",
              values_fill = 0)%>%
  as.data.frame()
# colnames are OTUs and rownames are sample IDs
rownames(ps.q.df.wide)<-ps.q.df.wide$Sample
ps.q.df.wide<-subset(ps.q.df.wide,select = -Sample)
ps.q.df.wide<-subset(ps.q.df.wide,select = -class)


# 3.2 Running MaAsLin 2 ####
if (comparison=="age"){
  maaslin.reference<-paste("agegroup",ref.level,sep = ",")
  maaslin.comparison<-c("agegroup")
}else if(comparison=="sex"){
  maaslin.reference<-paste("sex","F",sep = ",")
  maaslin.comparison<-"sex"
}else if(comparison=="strain"){
  maaslin.reference<-paste("class",ref.level,sep = ",")
  maaslin.comparison<-"class"
}

# Add relations
relations<-read.table("../amplicon_nmr/data/metadata/pooled-metadata/nmr-relations.tsv",
                      header = T,
                      sep = "\t")
custom.md<-custom.md%>%
  left_join(relations,by="Sample")
rownames(custom.md)<-custom.md$Sample

set.seed(1)
maaslin.fit_data = 
  Maaslin2(input_data = ps.q.df.wide, 
           input_metadata = custom.md, 
           min_prevalence = 0,
           analysis_method = "LM",
           normalization = "TSS", #TSS for Kraken2, CLR for HUMANn3
           transform = "LOG",
           random_effects = c("relation"), 
           standardize = FALSE,
           output = file.path("./output/maaslin2",paste0(pipeline.name,"-output"),
                              rare.status,paste(
                                paste(format(Sys.time(),format="%Y%m%d"),
                                      format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                host,filter.status,agglom.rank,comparison,
                                paste(custom.levels,collapse = '-'),
                                "ref",ref.level,sep = "-")), 
           fixed_effects =maaslin.comparison,
           reference = maaslin.reference,
           max_significance = 0.05)

save.image(file.path("./output/rdafiles",paste(
  paste(format(Sys.time(),format="%Y%m%d"),
        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
  "maaslin2",pipeline.name,rare.status,filter.status,host,agglom.rank,
  comparison,paste(custom.levels,collapse = '-'),"ref",
  ref.level,"workspace.RData",sep="-")))
q()

# Load output for plotting ######
library(tidyverse)
maaslin.date_time<-"20240621_18_31_55"
if(pipeline.name=="kraken2"){
  ps.q.agg.date_time<-"20241003_13_52_43"
}else if (pipeline.name=="singlem"){
  ps.q.agg.date_time<-"20240929_23_33_37"
}
host<-"NMR"
comparison<-"age"
ref.level<-"agegroup0_10" # choose the reference level
# this is for file names
if(host=="NMR"){
  host.labels<-c("NMR" = "*Heterocephalus glaber*")
}else{
  host.labels<-
    c("B6mouse" = "B6 mouse",
      "MSMmouse" = "MSM/Ms mouse",
      "FVBNmouse" = "FVB/N mouse")
}
custom.levels<-c("agegroup10_16","agegroup0_10")
metric.labs<-c("agegroup0_10"="Young naked mole-rats",
               "agegroup10_16"="Old naked mole-rats")
agglom.rank<-"Species"
rare.status<-"rare"
filter.status<-"nonfiltered"

load(file.path("./output/rdafiles",paste(
  maaslin.date_time,
  "maaslin2-kraken2",rare.status,filter.status,host,agglom.rank,
  comparison,paste(sort(custom.levels),collapse = '-'),"ref",
  ref.level,"workspace.RData",sep="-")))
ps.q.agg<-readRDS(file.path(rdafiles.directory,
                            paste(ps.q.agg.date_time,"phyloseq-kraken2",
                                  agglom.rank,"table.rds",sep = "-")))
source("../amplicon_nmr/code/r-scripts/make_features_maaslin.R")
if(min(maaslin.fit_data$results$qval)<0.05){
  maaslin.signif.features<-maaslin.fit_data$results%>%
    filter(qval<0.05) # should be qval
}else{
  maaslin.signif.features<-maaslin.fit_data$results%>%
    arrange(qval)%>%
    head(n = 10) # if no significant results found
}



foo<-ps.q.agg
foo$maaslin<-foo$Species
foo<-make_features_maaslin(foo,"maaslin")
foo<-unique(foo[,c("maaslin","Species")])
maaslin.signif.features<-maaslin.signif.features%>%
  left_join(foo[,c("maaslin","Species")],by=c("feature"="maaslin"))%>%
  distinct()
rm(foo)
maaslin.signif.features$feature<-maaslin.signif.features$Species
maaslin.signif.features<-subset(maaslin.signif.features, select=-Species)

table(maaslin.signif.features$feature%in%ps.q.agg$Species)

# Calculate SD and mean if some samples don't have a certain taxon ####
sample.vector<-rep(0,length.out=nrow(custom.md))
names(sample.vector)<-custom.md$Sample
abundance.df<-ps.q.agg%>%
  ungroup()%>%
  filter(Species=="Staphylococcus_shinii")%>%
  distinct(Sample,.keep_all = T)%>%
  select(Sample,RelativeAbundance)
abundance.vector<-abundance.df$RelativeAbundance
names(abundance.vector)<-abundance.df$Sample
sample.vector
new.sample.vector<-ifelse(names(sample.vector)%in%names(abundance.vector),abundance.vector,0)
names(new.sample.vector)<-names(sample.vector)
new.sample.vector
mean(new.sample.vector)
sd(new.sample.vector)

# I want to add zero rows ###
ps.q.agg.relab<-ps.q.agg%>%ungroup()
for (sample in unique(ps.q.agg.relab$Sample)) {
  missing_species <- setdiff(maaslin.signif.features$feature, ps.q.agg.relab$Species[ps.q.agg.relab$Sample == sample])
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
  left_join(unique(custom.md),by="Sample")%>%
  fill(age, .direction = "down")%>%
  fill(agegroup, .direction = "down")#%>%
  # fill(old_agegroup, .direction = "down")


# for (i in seq(nrow(maaslin.signif.features))){
#   pretty.feature_name<-gsub("_"," ",maaslin.signif.features$feature[i])
#   pretty.feature_name<-paste0("<i>",pretty.feature_name,"</i>")
#   taxon.plot<-ps.q.agg.relab%>%
#     # left_join(custom.md,by="Sample")%>%
#     group_by(Sample)%>%
#     filter(Species==maaslin.signif.features$feature[i])%>%
#     # bind_rows(foo)%>%
#     ggplot(aes(x=Sample,y=RelativeAbundance))+
#     geom_bar(stat="identity")+
#     facet_grid(~factor(agegroup,
#                        levels=custom.levels),
#                scales = "free_x",
#                labeller = as_labeller(metric.labs))+
#     theme_bw()+
#     ggtitle(paste(pretty.feature_name,"relative abundance"))+
#     labs(y="Relative abundance (%)")+
#     theme(axis.line = element_blank(), 
#           strip.text.x = ggtext::element_markdown(size = 20),# the name of 
#           # each facet will be recognised as a markdown object, so we can
#           # add line breaks (cause host names are too long)
#           axis.text.x = element_text(size=30),# rotate 
#           # the x-axis labels by 45 degrees and shift to the right
#           axis.text.y = element_text(size=30), # size of y axis ticks
#           axis.title = element_text(size = 30), # size of axis names
#           plot.title = ggtext::element_markdown(size = 35), # the plot 
#           # title will be recognised as a markdown object, so we can
#           # add line breaks (cause host names are too long)
#           plot.caption = element_text(size=23),# size of plot caption
#           legend.text = element_text(size = 20),# size of legend text
#           legend.title = element_text(size = 25), # size of legend title
#           legend.position = "none")
#   for(image.format in c("png","tiff")){
#     ggsave(paste0("./images/barplots/",
#                   paste(paste(format(Sys.time(),format="%Y%m%d"),
#                               format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                         "barplot","NMR",maaslin.signif.features$feature[i],
#                         agglom.rank,sep = "-"),".",image.format),
#            plot=taxon.plot,
#            width = 6000,height = 3000,
#            units = "px",dpi=300,device = image.format)
#   }
# }

# Check the average relative abundance of significant taxa ####
# Separate by agegroup
ps.q.agg%>%
  left_join(custom.md)%>%
  group_by(agegroup)%>% # group by class (animal host),
  mutate(TotalAgegroup=sum(Abundance))%>%
  group_by_at(c("agegroup",agglom.rank))%>%
  mutate(TotalAgglomRankAge=sum(Abundance))%>%
  mutate(MeanRelativeAbundanceAgegroup=TotalAgglomRankAge/TotalAgegroup*100)%>%
  filter(Species%in%maaslin.signif.features$feature)%>%
  select(MeanRelativeAbundance,MeanRelativeAbundanceAgegroup)


# Plot differentially abundant species ####
gg.labs.name<-"Age group"

sample.levels<-custom.md%>%
  filter(Sample%in%ps.q.agg$Sample)%>%
  select(Sample,age)%>%
  arrange(age)%>%
  distinct()
sample.levels$Sample<-factor(sample.levels$Sample,
                             levels=sample.levels$Sample)

diff.abund.plot<-ps.q.agg.relab%>%
  filter(Species%in%maaslin.signif.features$feature)%>%
  left_join(maaslin.signif.features[,c("feature","qval")],
            by=c("Species"="feature"))%>%
  mutate(Species=gsub("_"," ",Species),
         Species=paste0("<i>",Species,"</i>"," (p = ",round(qval,digits = 3),")"))%>%
  mutate(Sample=factor(Sample,levels=sample.levels$Sample))%>%
  group_by(class,Species)%>%
  left_join(custom.md[,c("Sample","agegroup","age")])%>%
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
  coord_cartesian(expand = c("bottom" = FALSE))+
  scale_color_manual(breaks = pretty.level.names,
                     labels=unname(pretty.level.names))+
  scale_x_discrete(labels=pretty.level.names,
                   limits=sample.levels$Sample)+ 
  scale_fill_manual(labels=pretty.level.names)+
  scale_fill_viridis_d(option = "C")+
  # ggtitle(paste0("Relative abundances of differentially abundant species in different naked mole-rat age groups"))+
  theme(
    # plot.margin=unit(c(1,1,1,2), 'cm'),
        axis.title = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        axis.text.y = ggtext::element_markdown(size=10),
        axis.text.x = element_text(size=10),
        strip.text.x = ggtext::element_markdown(size=10),
        plot.title = element_text(size = 17),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.position = "right",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())


for (image.format in c("png","tiff")){
  if(nrow(maaslin.signif.features)==1){
    # diff.abund.plot<-diff.abund.plot+
    #   ggtitle(paste0("Relative abundance of differentially abundant species\nin different naked mole-rat age groups"))
  }
  ggsave(paste0("./images/barplots/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "diffabund-bacteria-nmr",
                      sep = "-"),".",image.format),
         plot=diff.abund.plot,
         width=7, height=4,units="in",
         # width = ifelse(nrow(maaslin.signif.features)==1,4000,7000),
         # height = ifelse(nrow(maaslin.signif.features)==1,2000,9000), units = "px",
         dpi=300,device = image.format)
}

