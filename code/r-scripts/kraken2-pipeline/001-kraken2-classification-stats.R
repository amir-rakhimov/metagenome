#' ---
#' output: 
#'   bookdown::html_document2:
#'      toc: true
#' ---
#' 
#' ```{r, setup 001-kraken2-classification-stats.R, include=FALSE}
#' knitr::opts_knit$set(root.dir = '/home/rakhimov/projects/metagenome')
#' ```
#' # Whole metagenome sequencing (shotgun sequencing) analysis {-}
#+ echo=FALSE
# Create a table with classification statistics ####
#' # Create a table with classification statistics
#' 
#+ echo=FALSE
## Introduction ####
#'
#' ## Introduction
#' This script creates a table that shows percentage of unclassified reads 
#' per sample. The number of unclassified reads was taken from the job report.
# install.packages(c("tidyverse"))
library(tidyverse)
source(here::here("config/R/config.R"))# config file with global variables
dir.create(community.composition.tables,recursive = TRUE)
dir.create(community.composition.rdafiles,recursive = TRUE)
dir.create(community.composition.figures,recursive = TRUE)
kraken2.results<-read.table(kraken2.unclassified.stats.fname,header = T)
kraken2.results<-kraken2.results%>%
  mutate(ClassifiedRate=Classified/Total*100)%>%
  mutate(meanClassifiedRate=round(mean(ClassifiedRate)))%>%
  mutate(SDClassifiedRate=round(sd(ClassifiedRate),3))
knitr::kable(kraken2.results, format = "simple")%>%
  print()

classification.stats.fname <- file.path(community.composition.tables,
                                        "kraken2-classification-stats.tsv")
if (! file.exists(classification.stats.fname)){
  write.table(kraken2.results, classification.stats.fname,
              row.names = F,quote = F,sep = "\t")
}
sessionInfo()
rm(list =setdiff(ls(all.names = TRUE), c("markdown.dir")))
gc()