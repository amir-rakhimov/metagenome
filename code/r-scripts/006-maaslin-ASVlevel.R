library(tidyverse)
library(phyloseq)
library(Maaslin2)
library(vegan)
# choose what to compare
comparison<-"age"
# comparison<-"sex"
# comparison<-"strain"
# choose the host of interest
host<-"NMR"
# host<-"mice"
ref.level<-"between2and10" # choose the reference level
# this is for file names
if(host=="NMR"){
  host.labels<-c("NMR" = "*Heterocephalus glaber*")
}else{
  host.labels<-
    c("B6mouse" = "B6 mouse",
      "MSMmouse" = "MSM/Ms mouse",
      "FVBNmouse" = "FVB/N mouse")
}
custom.levels<-c("NMR"
)
agglom.rank<-"Species"
load(paste0("./output/rdafiles/",paste("kraken2",
                                       agglom.rank,
                                       "phyloseq-workspace.RData",sep = "-")))
rare.status<-"rare"
filter.status<-"nonfiltered"

# Import data ####
ps.q.df.preprocessed<-read.table(paste0("./output/rtables/","kraken2-",
                                        "ps.q.df.rare.nonfiltered-",
                                        agglom.rank,"-",
                                        paste(custom.levels,collapse = '-'),
                                        ".tsv"),
                                 header = T,sep = "\t")
ps.q.df.preprocessed$Sample<-as.factor(ps.q.df.preprocessed$Sample)
ps.q.df.preprocessed$class<-as.factor(ps.q.df.preprocessed$class)
ps.q.df.preprocessed$sex<-as.factor(ps.q.df.preprocessed$sex)
ps.q.df.preprocessed$birthday<-as.Date(ps.q.df.preprocessed$birthday)

# Break the data into age groups
age.breaks<-c(min(custom.md$age),10,max(custom.md$age))
custom.md<-custom.md%>%
  mutate(age_group=cut(age, breaks = age.breaks, 
                       include.lowest = TRUE))%>%
  mutate(age_group=as.factor(age_group))
if(host=="NMR"){
  # select nmr and add age groups
  ps.q.df.preprocessed<-ps.q.df.preprocessed%>%
    mutate(age_group=cut(age, breaks = age.breaks, 
                        include.lowest = TRUE))
  # we create these new levels because maaslin is itsy bitsy
  unique_levels <- custom.md %>%
    distinct(age_group) %>%
    arrange(age_group) %>%
    mutate(new_age_group=gsub("\\[|\\]|\\(","",age_group))%>%
    separate(new_age_group,c("leftbound","rightbound"))%>%
    mutate(new_age_group = paste0("between", leftbound,"and",rightbound))%>%
    select(-leftbound,-rightbound)
  ps.q.df.preprocessed <- ps.q.df.preprocessed %>%
    left_join(unique_levels, by = "age_group")
  colnames(ps.q.df.preprocessed)[which(colnames(ps.q.df.preprocessed)=="age_group")]<-"old_age_group"
  colnames(ps.q.df.preprocessed)[which(colnames(ps.q.df.preprocessed)=="new_age_group")]<-"age_group"
  # add age group to metadata
  custom.md$Sample<-rownames(custom.md)
  # Extract unique levels from the original_vector
  unique_levels <- custom.md %>%
    distinct(age_group) %>%
    arrange(age_group) %>%
    mutate(new_age_group=gsub("\\[|\\]|\\(","",age_group))%>%
    separate(new_age_group,c("leftbound","rightbound"))%>%
    mutate(new_age_group = paste0("between", leftbound,"and",rightbound))%>%
    select(-leftbound,-rightbound)
  custom.md <- custom.md %>%
    left_join(unique_levels, by = "age_group")
  colnames(custom.md)[which(colnames(custom.md)=="age_group")]<-"old_age_group"
  colnames(custom.md)[which(colnames(custom.md)=="new_age_group")]<-"age_group"
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
  ps.q.df.preprocessed$age_group<-ifelse(ps.q.df.preprocessed$class=="B6mouse","B6",
                                         ifelse(grepl("2020",ps.q.df.preprocessed$birthday),"old","young"))
  # also to metadata
  custom.md<-custom.md%>%
    filter(class%in%custom.levels)
  custom.md$age_group<-ifelse(custom.md$class=="B6mouse","B6",
                             ifelse(grepl("2020",custom.md$birthday),"old","young"))
}
# Creating custom levels ####
if (comparison=="age"){
  # names for levels are age groups
  pretty.facet.labels<-names(table(ps.q.df.preprocessed$age_group))
  names(pretty.facet.labels)<-names(table(ps.q.df.preprocessed$age_group))
  custom.levels<-names(pretty.facet.labels)
  
}else if (comparison=="sex"){
  pretty.facet.labels<-
    c("F" = "Females",
      "M" = "Males")
  custom.levels<-names(pretty.facet.labels)
  
}else if(comparison=="strain"){
  pretty.facet.labels<-
    c("B6mouse" = "B6 mouse",
      "MSMmouse" = "MSM/Ms mouse",
      "FVBNmouse" = "FVB/N mouse"
    )
  custom.levels<-intersect(names(pretty.facet.labels),custom.md$class)
  pretty.facet.labels<-pretty.facet.labels[which(names(pretty.facet.labels)%in%custom.levels)]
}
# Preparing the dataset ####
# filter the dataset
ps.q.df <-ps.q.df.preprocessed%>%
  dplyr::select(all_of(c("Sample","Abundance","class","age_group","sex",agglom.rank)))%>%
  filter(Abundance!=0)
ps.q.df.maaslin.input.wide<-ps.q.df%>%
  pivot_wider(names_from = agglom.rank, # or OTU
              values_from = "Abundance",
              values_fill = 0)%>%
  as.data.frame()
# colnames are OTUs and rownames are sample IDs
rownames(ps.q.df.maaslin.input.wide)<-ps.q.df.maaslin.input.wide$Sample
ps.q.df.maaslin.input.wide<-ps.q.df.maaslin.input.wide[,-c(1:4)]  


# 3.2 Running MaAsLin 2 ####
if (comparison=="age"){
  maaslin.reference<-paste("age_group",ref.level,sep = ",")
  maaslin.comparison<-"age_group"
}else if(comparison=="sex"){
  maaslin.reference<-paste("sex","F",sep = ",")
  maaslin.comparison<-"sex"
}else if(comparison=="strain"){
  maaslin.reference<-paste("class",ref.level,sep = ",")
  maaslin.comparison<-"class"
}

set.seed(1)
maaslin.fit_data = 
  Maaslin2(input_data = ps.q.df.maaslin.input.wide, 
           input_metadata = custom.md, 
           min_prevalence = 0,
           normalization = "TSS",
           transform = "AST",
           analysis_method = "LM",
           random_effects = NULL,
           standardize = FALSE,
           output = paste0("./output/maaslin2/","kraken2","-output/",
                          rare.status,"/",paste(host,filter.status,agglom.rank,
                                                comparison,
                                                ref.level,sep="-")), 
           fixed_effects = maaslin.comparison,
           reference = maaslin.reference,
           max_significance = 0.05)


save.image(paste0("./output/rdafiles/",
                  paste("maaslin",rare.status,filter.status,host,agglom.rank,
                        comparison,ref.level,
                  "workspace.RData",sep="-")))
q()

load(paste0("./output/rdafiles/",
            paste("maaslin",rare.status,filter.status,host,agglom.rank,
                  comparison,ref.level,
                  "workspace.RData",sep="-")))
source("../amplicon_nmr/code/r-scripts/make_features_maaslin.R")
maaslin.signif.features<-maaslin.fit_data$results%>%
  filter(qval<0.05) # should be qval


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


for (i in seq(nrow(maaslin.signif.features))){
  taxon.plot<-ps.q.agg%>%
    left_join(custom.md,by="Sample")%>%
    group_by(Sample)%>%
    filter(Species==maaslin.signif.features$feature[i])%>%
    ggplot(aes(x=Sample,y=RelativeAbundance))+
    ylab("Relative Abundance (%)")+
    geom_bar(stat="identity")+
    ggtitle(paste(maaslin.signif.features$feature[i],"relative abundance"))+
    facet_grid(~factor(age_group,levels=rev(custom.levels)),scales = "free_x")+
    theme_bw()+
    theme(axis.line = element_blank(), 
          strip.text.x = ggtext::element_markdown(size = 20),# the name of 
          # each facet will be recognised as a markdown object, so we can
          # add line breaks (cause host names are too long)
          axis.text.x = element_text(size=20),# rotate 
          # the x-axis labels by 45 degrees and shift to the right
          axis.text.y = element_text(size=20), # size of y axis ticks
          axis.title = element_text(size = 20), # size of axis names
          plot.title = ggtext::element_markdown(size = 25), # the plot 
          # title will be recognised as a markdown object, so we can
          # add line breaks (cause host names are too long)
          plot.caption = element_text(size=23),# size of plot caption
          legend.text = element_text(size = 20),# size of legend text
          legend.title = element_text(size = 25), # size of legend title
          legend.position = "right")
  ggsave(paste0("./images/barplots/",
                paste(Sys.Date(),"barplot","NMR",maaslin.signif.features$feature[i],
                      agglom.rank,sep = "-"),".png"),
         plot=taxon.plot,
         width = 6000,height = 3000,
         units = "px",dpi=300,device = "png")
  ggsave(paste0("./images/barplots/",
                paste(Sys.Date(),"barplot","NMR",maaslin.signif.features$feature[i],
                      agglom.rank,sep = "-"),".tiff"),
         plot=taxon.plot,
         width = 6000,height = 3000,
         units = "px",dpi=300,device = "tiff")
}



