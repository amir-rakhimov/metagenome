library(tidyverse)
#' Discard MAGs with >10% contamination
#' Three categories:
#' 1. High quality MAGs: >=90% complete, <=5% contamination
#' 2. Medium quality MAGs: <100% but >=50% complete, <10% contamination
#' 3. Low quality MAGs: <50% complete, <10% contamination
#' 
combined.mag_stats.dir<-"./output/mag_assembly/seqkit_output"
combined.mag_stats.date_time<-"20250430_15_23_58"
combined.mag_stats.fname<-file.path(combined.mag_stats.dir,
                             paste(combined.mag_stats.date_time,
                             "all_samples_stats.tsv",sep="_"))
checkm2.output.date_time<-"20250501_07_53_55"
checkm2.output.dir.all<-"./output/mag_assembly/checkm2_output"
metabat2.output.date_time<-"20250417_23_21_29"
#' MAG stats analysis
combined.mag_stats<-read.table(combined.mag_stats.fname,header = TRUE,
                           sep="\t",comment.char = "")
head(combined.mag_stats)
#' Number of contigs: 156252
combined.mag_stats%>%
  distinct(contig_id)%>%
  nrow

#'Calculate assembly status of each sample:
#'Total length of contigs by sample, longest contig, shortest contig,
#'average contig length, number of bins per sample,
#'number of contigs per sample,
#'total length of contigs in all samples combined.
combined.mag_stats%>%
  group_by(Sample)%>%
  mutate(total_len=sum(contig_length),
         longest_contig_len=max(contig_length),
         shortest_contig_len=min(contig_length),
         average_len=mean(contig_length),
         n_bin=n_distinct(Bin),
         n_contig=n_distinct(contig_id))%>%
  ungroup()%>%
  distinct(Sample,.keep_all = T)%>%
  select(-Bin,-contig_id,-contig_length,-GC,-AT)%>%
  mutate(total_len_all=sum(total_len))

#'Number of total bins (MAGs): 1347
combined.mag_stats%>%
  distinct(Sample,Bin)%>%
  nrow

#'Length distribution
ggplot(combined.mag_stats,aes(x=contig_length))+
  geom_histogram()

#'Number of contigs per bin
contigs.per.bin<-combined.mag_stats%>%
  group_by(Sample,Bin)%>%
  summarise(n_contig_per_bin=n_distinct(contig_id))%>%
  group_by(Sample)%>%
  summarise(mean_n_contig=mean(n_contig_per_bin))



#' CheckM2 output analysis
#' match the number at the start of sample name (2 in 2D10)
#' match the letters (D in 2D10)
#' match numbers in the sample name (10 in 2D10)
#' match a letter at the end (b in Y51b)
checkm2.output.sample.dirs<-list.files(checkm2.output.dir.all)[grep(paste0(
  checkm2.output.date_time,"_[0-9]*[A-Z]+[0-9]+[a-z]*_checkm2$"),
  list.files(checkm2.output.dir.all))]
checkm2.output.sample.dirs.contents<-list.files(file.path(checkm2.output.dir.all,
                                                          checkm2.output.sample.dirs),
                                                full.names = T)
checkm2.quality.reports.paths<-checkm2.output.sample.dirs.contents[grep(
  "quality_report\\.tsv",checkm2.output.sample.dirs.contents)]

checkm2.quality.reports.combined<-do.call(rbind, 
                                          lapply(checkm2.quality.reports.paths, 
                                                 read.table, 
                                                 header = TRUE, 
                                                 sep = "\t"))%>%
  as_tibble()
head(checkm2.quality.reports.combined)
#' MAGs with >10% contamination: 45 out of 1347 MAGs
contaminated.mags<-checkm2.quality.reports.combined%>%
  filter(Contamination>10)

#' High-quality MAGs (>=90% complete, <=5% contamination): 313
high.quality.mags<-checkm2.quality.reports.combined%>%
  filter(Completeness>=90,
         Contamination<=5)
#' Medium-quality MAGs (<100% but >=50% complete, <10% contamination):509
medium.quality.mags<-checkm2.quality.reports.combined%>%
  filter(Completeness<100,
         Completeness>=50,
         Contamination<10,
         !Name %in%high.quality.mags$Name)
#' Low-quality MAGs (<50% complete, <10% contamination): 480
low.quality.mags<-checkm2.quality.reports.combined%>%
  filter(Completeness<50,
         Contamination<10)
selected.mags<-rbind(high.quality.mags,medium.quality.mags,low.quality.mags)%>%distinct(Name)%>%pull(Name)
#' Sanity check: are there any MAGs we didn't include that aren't contaminated?
checkm2.quality.reports.combined%>%
  filter(Name %in%setdiff(checkm2.quality.reports.combined$Name,selected.mags),
         !Name%in%contaminated.mags$Name)%>%
  pull(Name)

#' Work with high-quality MAGs now
high.quality.mags.bins<-high.quality.mags%>%
  mutate(Bin=str_match(Name,"bin\\.[0-9]+"),
         Bin=gsub("bin\\.","",Bin),
         Bin=as.numeric(Bin),
         Sample=gsub(paste0(metabat2.output.date_time,"_"),"",Name),
         Sample=gsub("_bin\\.[0-9]+","",Sample))%>%
  select(Name,Bin,Sample)
high.quality.mags.bins%>%
  head

# Filter to keep contigs from high-quality MAGs
semi_join(combined.mag_stats,high.quality.mags.bins,by=c("Sample","Bin"))%>%View
