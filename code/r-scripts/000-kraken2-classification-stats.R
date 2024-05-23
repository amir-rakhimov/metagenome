library(tidyverse)
date_time<-"20240409_17_32_40"
kraken2.results<-read.table(file.path("./output/kraken2_pipeline",
                                     paste(date_time,"unclassified_reads.tsv",
                                            sep="_")),header = T)
kraken2.results<-kraken2.results%>%
  mutate(ClassifiedRate=Classified/Total*100)%>%
  mutate(meanClassifiedRate=round(mean(ClassifiedRate)))%>%
  mutate(SDClassifiedRate=round(sd(ClassifiedRate),3))
write.table(kraken2.results,file.path("./output/rtables",
                      paste(date_time,"kraken2-classification-stats.tsv",
                            sep="_")),row.names = F,quote = F,sep = "\t")
