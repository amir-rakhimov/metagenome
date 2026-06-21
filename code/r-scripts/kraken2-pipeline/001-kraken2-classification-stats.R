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
kraken2.results<-read.table(kraken2.unclassified.stats,header = T)
kraken2.results<-kraken2.results%>%
  mutate(ClassifiedRate=Classified/Total*100)%>%
  mutate(meanClassifiedRate=round(mean(ClassifiedRate)))%>%
  mutate(SDClassifiedRate=round(sd(ClassifiedRate),3))
knitr::kable(kraken2.results, format = "simple")%>%
  print()

# write.table(kraken2.results,file.path("./output/rtables",
#                       paste(date_time,"kraken2-classification-stats.tsv",
#                             sep="_")),row.names = F,quote = F,sep = "\t")
