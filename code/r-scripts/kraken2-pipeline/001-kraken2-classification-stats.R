#' ---
#' output: 
#'   bookdown::html_document2:
#'      toc: true
#' ---
#' 
#' ```{r, setup 001-kraken2-classification-stats.R, include=FALSE}
#' knitr::opts_knit$set(root.dir = '/home/rakhimov/projects/metagenome')
#' ```
#' ```{r, echo = FALSE}
#' # For showing images, tables, etc: Use global path
#' # knitr::spin("code/r-scripts/kraken2-pipeline/001-kraken2-classification-stats.R", 
#' #             knit = FALSE)
#' # file.rename("code/r-scripts/kraken2-pipeline/001-kraken2-classification-stats.Rmd", 
#' #             "markdown/001-kraken2-classification-stats.Rmd")
#' # rmarkdown::render('./markdown/001-kraken2-classification-stats.Rmd', 
#' #                   'html_document',
#' #                   knit_root_dir="/home/rakhimov/projects/metagenome/")
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
date_time<-"20240409_17_32_40"
kraken2.results<-read.table(file.path("./output/kraken2_pipeline",
                                     paste(date_time,"unclassified_reads.tsv",
                                            sep="_")),header = T)
kraken2.results<-kraken2.results%>%
  mutate(ClassifiedRate=Classified/Total*100)%>%
  mutate(meanClassifiedRate=round(mean(ClassifiedRate)))%>%
  mutate(SDClassifiedRate=round(sd(ClassifiedRate),3))
kraken2.results
# write.table(kraken2.results,file.path("./output/rtables",
#                       paste(date_time,"kraken2-classification-stats.tsv",
#                             sep="_")),row.names = F,quote = F,sep = "\t")
