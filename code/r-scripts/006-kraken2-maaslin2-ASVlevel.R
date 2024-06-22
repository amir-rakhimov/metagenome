library(tidyverse)
library(phyloseq)
library(Maaslin2)
library(vegan)
# phyloseq.date_time<-"20240619_14_11_06"
ps.q.df.preprocessed.date_time<-"20240619_20_17_34"
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
if(host=="NMR"){
  host.labels<-c("NMR" = "*Heterocephalus glaber*")
}else{
  host.labels<-
    c("B6mouse" = "B6 mouse",
      "MSMmouse" = "MSM/Ms mouse",
      "FVBNmouse" = "FVB/N mouse")
}
custom.levels<-c("NMR")
agglom.rank<-"Species"
# Load the Workspace from phyloseq (output of 001-phyloseq-qiime2.R)
# load(file.path("./output/rdafiles",paste(
#   phyloseq.date_time,"kraken2",agglom.rank,
#   "phyloseq-workspace.RData",sep = "-")))
custom.md<-readRDS("../amplicon_nmr/output/rdafiles/custom.md.rds")

rare.status<-"rare"
filter.status<-"nonfiltered"

# Import data ####
ps.q.df.preprocessed<-read.table(
  file.path("./output/rtables",paste0(
    paste(ps.q.df.preprocessed.date_time,
          "kraken2-ps.q.df.rare-nonfiltered",agglom.rank,
          paste(host,collapse = '-'),sep = "-"),".tsv")),
  header = T,sep = "\t")

ps.q.df.preprocessed$Sample<-as.factor(ps.q.df.preprocessed$Sample)
ps.q.df.preprocessed$class<-as.factor(ps.q.df.preprocessed$class)
ps.q.df.preprocessed$sex<-as.factor(ps.q.df.preprocessed$sex)
ps.q.df.preprocessed$birthday<-as.Date(ps.q.df.preprocessed$birthday)

# Break the data into age groups
# age.breaks<-c(min(custom.md$age),10,max(custom.md$age))
age.breaks<-c(0,10,16)
if(host=="NMR"){
  # select nmr and add age groups
  ps.q.df.preprocessed<-ps.q.df.preprocessed%>%
    mutate(agegroup=cut(age, breaks = age.breaks, 
                         right = FALSE))
  # we create these new levels because maaslin is itsy bitsy
  unique_levels <- ps.q.df.preprocessed %>%
    ungroup()%>%
    distinct(agegroup)%>%
    arrange(agegroup) %>%
    mutate(new_agegroup = paste0("agegroup", agegroup))%>%
    mutate(new_agegroup = gsub("\\(|\\)|\\[|\\]","",new_agegroup))%>%
    mutate(new_agegroup = gsub("\\,","_",new_agegroup))
  ps.q.df.preprocessed <- ps.q.df.preprocessed %>%
    left_join(unique_levels, by = "agegroup")
  colnames(ps.q.df.preprocessed)[which(colnames(ps.q.df.preprocessed)=="agegroup")]<-"old_agegroup"
  colnames(ps.q.df.preprocessed)[which(colnames(ps.q.df.preprocessed)=="new_agegroup")]<-"agegroup"
  # add age group to metadata
  custom.md$Sample<-rownames(custom.md)
  custom.md<-custom.md%>% 
    filter(Sample %in% ps.q.df.preprocessed$Sample)%>%
    group_by(Sample)%>%
    mutate(birthday=as.Date(birthday))%>%
    mutate(age=year(as.period(interval(birthday,as.Date("2023-11-16")))))%>%
    left_join(unique(ps.q.df.preprocessed[,c("Sample","agegroup")]),by="Sample")%>%
    mutate(agegroup=as.factor(agegroup))%>%
    as.data.frame()
  rownames(custom.md)<-custom.md$Sample
}else if(host=="mice"){
  # select mice and add age groups: B6, old, or young
  # B6 are separate
  # mice born before 2020 are old
  # after 2023 are young
  custom.levels<-c("B6mouse",
                   "MSMmouse",
                   "FVBNmouse")
  ps.q.df.preprocessed<-ps.q.df.preprocessed%>%
    filter(class%in%custom.levels,Abundance!=0)
  ps.q.df.preprocessed$agegroup<-ifelse(ps.q.df.preprocessed$class=="B6mouse","B6",
                                         ifelse(grepl("2020",ps.q.df.preprocessed$birthday),"old","young"))
  # also to metadata
  custom.md<-custom.md%>%
    filter(class%in%custom.levels)
  custom.md$agegroup<-ifelse(custom.md$class=="B6mouse","B6",
                             ifelse(grepl("2020",custom.md$birthday),"old","young"))
}
# Creating custom levels ####
if (comparison=="age"){
  # names for levels are age groups
  pretty.level.names<-names(table(ps.q.df.preprocessed$old_agegroup))
  names(pretty.level.names)<-names(table(ps.q.df.preprocessed$agegroup))
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
  dplyr::select(all_of(c("Sample","Abundance","class","agegroup","sex",agglom.rank)))%>%
  filter(Abundance!=0)
ps.q.df.wide<-ps.q.df%>%
  pivot_wider(names_from = agglom.rank, # or OTU
              values_from = "Abundance",
              values_fill = 0)%>%
  as.data.frame()
# colnames are OTUs and rownames are sample IDs
rownames(ps.q.df.wide)<-ps.q.df.wide$Sample
ps.q.df.wide<-ps.q.df.wide[,-c(1:4)]  


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
           normalization = "TSS",
           transform = "LOG",
           analysis_method = "LM",
           random_effects = c("relation"), 
           standardize = FALSE,
           output = file.path("./output/maaslin2","kraken2-output",
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
  "maaslin2-kraken2",rare.status,filter.status,host,agglom.rank,
  comparison,paste(custom.levels,collapse = '-'),"ref",
  ref.level,"workspace.RData",sep="-")))
q()

# Load output for plotting ######
library(tidyverse)
maaslin.date_time<-"20240621_18_31_55"
ps.q.agg.date_time<-"20240620_12_30_49"
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
ps.q.agg<-readRDS(file.path("output/rdafiles",
                            paste(ps.q.agg.date_time,"phyloseq-kraken2",
                                  agglom.rank,"table.Rda",sep = "-")))
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
  left_join(unique(ps.q.df.preprocessed[,c("Sample","agegroup","old_agegroup")]),by="Sample")%>%
  fill(age, .direction = "down")%>%
  fill(agegroup, .direction = "down")%>%
  fill(old_agegroup, .direction = "down")


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


# Plot differentially abundant species ####
set.seed(1)
library(Polychrome)
custom.fill<-createPalette(length(custom.levels),
                           seedcolors = c("#EE2C2C","#5CACEE","#00CD66",
                                          "#FF8C00","#BF3EFF", "#00FFFF",
                                          "#FF6EB4","#00EE00","#EEC900",
                                          "#FFA07A"))
names(custom.fill)<-custom.levels
swatch(custom.fill)

gg.labs.name<-"Age group"



sample.levels<-ps.q.agg.relab%>%
  select(Sample,agegroup)%>%
  arrange(agegroup)%>%
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
  scale_color_manual(breaks = pretty.level.names,
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
  ggtitle(paste0("Relative abundance of differentially abundant species in different naked mole-rat age groups"))


for (image_format in c("png","tiff")){
  if(nrow(maaslin.signif.features)==1){
    diff.abund.plot<-diff.abund.plot+
      ggtitle(paste0("Relative abundance of differentially abundant species\nin different naked mole-rat age groups"))
  }
  ggsave(paste0("./images/barplots/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      "diffabund-bacteria-nmr",
                      sep = "-"),".",image_format),
         plot=diff.abund.plot,
         width = ifelse(nrow(maaslin.signif.features)==1,4000,7000),
         height = ifelse(nrow(maaslin.signif.features)==1,2000,9000),
         units = "px",dpi=300,device = image_format)
}

