library(tidyverse)
combined.stats.fname<-"./output/mag_assembly/seqkit_output/20250211_22_55_27_all_samples_stats.tsv"
combined.stats<-read.table(combined.stats.fname,header = TRUE,
                           sep="\t",comment.char = "")
mean.gc.table<-combined.stats%>%
  group_by(Sample,Bin)%>%
  summarise(mean_gc=mean(GC))%>%
  arrange(mean_gc)

combined.stats%>%
  left_join(mean.gc.table)%>%
  filter(contig_length>=50000)%>%
  filter(AT>70)%>%
  group_by(Sample,Bin)%>%
  unique()%>%
  tally%>%
  group_by(Sample)%>%
  unique()%>%
  tally()%>%
  summarise(ss=sum(n))
