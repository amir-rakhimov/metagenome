#' ---
#' output: 
#'   bookdown::html_document2:
#'      toc: true
#' ---

#' ```{r, setup 007-diffabund-tests.R, include=FALSE}
#' knitr::opts_knit$set(root.dir = '/home/rakhimov/projects/metagenome')
#' ```
#+ echo=FALSE
# Differential microbial abundance tests with MaAsLin2 ####
#' # Differential microbial abundance tests with MaAsLin2
#+ echo=FALSE
## Introduction ####
#'
#' ## Introduction
#' This script performs differential microbial abundance tests on 
#' different hosts. We will use MaAsLin2.
#+ echo=FALSE
## 1. Load necessary libraries. ####
#'
#' ## Load necessary libraries.
# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# Current: Bioconductor version 3.20 (BiocManager 1.30.27), R 4.4.3 (2025-02-28 ucrt)
#' Maaslin2: 1.20.0
# BiocManager::install(c("Maaslin2","phyloseq"), version = "3.20")
# BiocManager::install("phyloseq", version = "3.17")
# install.packages(c("tidyverse","vegan"))
library(tidyverse)
library(phyloseq)
library(Maaslin2)
library(vegan)
#+ echo=FALSE
## 2. Specifying parameters and directory/file names. #### 
#'
#' ## Specifying parameters and directory/file names. 
#' Directories with input files:
rdafiles.directory<-"./output/rdafiles"
rtables.directory<-"./output/rtables"
metadata.directory<-"../amplicon_nmr/output/rdafiles"
#' Import the rarefied abundance table (RDS):
ps.q.df.preprocessed.date_time<-"20260227_17_58_35"
#' Choose what to compare:
comparison<-"age"
#' Choose the reference level:
ref.level<-"agegroup0_10" 
agglom.rank<-"Species"
rare.status<-"rare"
filter.status<-"nonfiltered"
#' Add metadata:
custom.md<-readRDS(file.path(metadata.directory,"custom.md.ages.rds"))%>%
  filter(sequencing_type == "Naked mole-rat whole metagenome sequencing")
#+ echo=FALSE
## 3. Import the rarefied dataset (RDS). ####
#'
#' ## 3. Import the rarefied dataset (RDS). 
ps.q.df.preprocessed<-
  readRDS(file = file.path(rdafiles.directory,
                           paste(ps.q.df.preprocessed.date_time,
                                 "kraken2","ps.q.df.rare-nonfiltered",agglom.rank,
                                 "NMR.rds",sep = "-")) )%>%
  as_tibble()%>%
  mutate(class = factor(class))

#' Add age groups:
ps.q.df.preprocessed <- ps.q.df.preprocessed%>%
  left_join(custom.md[,c("Sample", "agegroup")])

#+ echo=FALSE
## 4. Prepare the data for test. ####
#'
#' ## 4. Prepare the data for test. 
if (comparison=="age"){
  # names for levels are age groups
  pretty.level.names<-names(table(custom.md$old_agegroup))
  names(pretty.level.names)<-names(table(custom.md$agegroup))
  custom.levels<-names(pretty.level.names)
  
}else if (comparison=="sex"){
  pretty.level.names<-
    c("female" = "Females",
      "male" = "Males")
  custom.levels<-names(pretty.level.names)
  
}
#+ echo=FALSE
## 5. Filter the dataset and convert into a wide format. ####
#'
#' ## Filter the dataset and convert into a wide format.
ps.q.df <-ps.q.df.preprocessed%>%
  dplyr::select(all_of(c("Sample","Abundance","class",agglom.rank)))%>%
  filter(Abundance!=0)
ps.q.df.wide<-ps.q.df%>%
  pivot_wider(names_from = agglom.rank, # or OTU
              values_from = "Abundance",
              values_fill = 0)%>%
  as.data.frame()
#' Colnames are Species and rownames are sample IDs
rownames(ps.q.df.wide)<-ps.q.df.wide$Sample
ps.q.df.wide<-subset(ps.q.df.wide,select = -Sample)
ps.q.df.wide<-subset(ps.q.df.wide,select = -class)

#+ echo=FALSE
## 5. Run MaAsLin 2. ####
#'
#' ## Run MaAsLin 2.
if (comparison=="age"){
  maaslin.reference<-paste("agegroup",ref.level,sep = ",")
  maaslin.comparison<-c("agegroup")
}else if(comparison=="sex"){
  maaslin.reference<-paste("sex","F",sep = ",")
  maaslin.comparison<-"sex"
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
           output = file.path("./output/maaslin2",paste0("kraken2","-output"),
                              rare.status,paste(
                                paste(format(Sys.time(),format="%Y%m%d"),
                                      format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                "NMR",filter.status,agglom.rank,comparison,
                                paste(custom.levels,collapse = '-'),
                                "ref",ref.level,sep = "-")), 
           fixed_effects =maaslin.comparison,
           reference = maaslin.reference,
           max_significance = 0.05)
#' Save the fit_data object as an RDS file.
# saveRDS(maaslin.fit_data,
#         file = file.path("output/rdafiles",
#                          paste(
#                            paste(format(Sys.time(),format="%Y%m%d"),
#                                  format(Sys.time(),format = "%H_%M_%S"),
#                                  sep = "_"),
#                            "kraken2-maaslin.fit_data.age-NMR-Species-age-ref",
#                            ref.level,".rds",sep = "-")))

# save.image(file.path("./output/rdafiles",paste(
#   paste(format(Sys.time(),format="%Y%m%d"),
#         format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#   "maaslin2","kraken2",rare.status,filter.status,"NMR",agglom.rank,
#   comparison,paste(custom.levels,collapse = '-'),"ref",
#   ref.level,"workspace.RData",sep="-")))

#+ echo=FALSE
## 6. Plotting the output. ####
#'
#' ## Plotting the output.
# maaslin.date_time<-"20240621_18_31_55"
#' Load the unrarefied abundance table.
ps.q.agg.date_time<-"20241003_13_52_43"
metric.labs<-c("agegroup0_10"="Young naked mole-rats",
               "agegroup10_16"="Old naked mole-rats")
# load(file.path("./output/rdafiles",paste(
#   maaslin.date_time,
#   "maaslin2-kraken2",rare.status,filter.status,"NMR",agglom.rank,
#   comparison,paste(sort(custom.levels),collapse = '-'),"ref",
#   ref.level,"workspace.RData",sep="-")))
ps.q.agg<-readRDS(file.path(rdafiles.directory,
                            paste(ps.q.agg.date_time,"phyloseq-kraken2",
                                  agglom.rank,"table.rds",sep = "-")))
source("../amplicon_nmr/code/r-scripts/make_features_maaslin.R")
#' Extract features with qvalue < 0.05.
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

#' Plot differentially abundant species.
gg.labs.name<-"Age group"

sample.levels<-custom.md%>%
  filter(Sample%in%ps.q.agg$Sample)%>%
  select(Sample,age)%>%
  arrange(age)%>%
  distinct()%>%
  mutate(Sample=factor(Sample,
                       levels=Sample))

diff.abund.plot<-ps.q.agg%>%
  mutate(Sample=factor(Sample,levels=sample.levels$Sample))%>%
  filter(Species%in%maaslin.signif.features$feature)%>%
  left_join(maaslin.signif.features[,c("feature","qval")],
            by=c("Species"="feature"))%>%
  mutate(Species=gsub("_"," ",Species),
         Species=paste0("<i>",Species,"</i>"," (p = ",round(qval,digits = 3),")"))%>%
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
  scale_fill_viridis_d(option = "C")+
  # ggtitle(paste0("Relative abundances of differentially abundant species in different naked mole-rat age groups"))+
  theme(
    # plot.margin=unit(c(1,1,1,2), 'cm'),
        axis.title = element_text(size = 10),
        axis.title.y = element_text(size = 12),
        axis.text.y = ggtext::element_markdown(size=10),
        axis.text.x = element_text(size=10),
        strip.text.x = ggtext::element_markdown(size=10),
        plot.title = element_text(size = 13),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.position = "right",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

#+ fig.width=7, fig.height=4,
print(diff.abund.plot +
        ggtitle(paste0("Relative abundances of differentially abundant species in different naked mole-rat age groups")))
# for (image.format in c("png","tiff")){
#   if(nrow(maaslin.signif.features)==1){
#     # diff.abund.plot<-diff.abund.plot+
#     #   ggtitle(paste0("Relative abundance of differentially abundant species\nin different naked mole-rat age groups"))
#   }
#   ggsave(paste0("./images/barplots/",
#                 paste(paste(format(Sys.time(),format="%Y%m%d"),
#                             format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                       "diffabund-bacteria-nmr",
#                       sep = "-"),".",image.format),
#          plot=diff.abund.plot,
#          width=7, height=4,units="in",
#          # width = ifelse(nrow(maaslin.signif.features)==1,4000,7000),
#          # height = ifelse(nrow(maaslin.signif.features)==1,2000,9000), units = "px",
#          dpi=300,device = image.format)
# }

sessionInfo()
rm(list = ls(all=TRUE))
gc()