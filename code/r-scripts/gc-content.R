library(tidyverse)
date_time="20250303_17_54_26"
combined.stats.fname<-paste0("./output/mag_assembly/seqkit_output/",date_time,
                             "_all_samples_stats.tsv")
combined.stats<-read.table(combined.stats.fname,header = TRUE,
                           sep="\t",comment.char = "")
mean.gc_at.table<-combined.stats%>%
  group_by(Sample,Bin)%>%
  summarise(mean_gc=mean(GC),
            mean_at=mean(AT))%>%
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
