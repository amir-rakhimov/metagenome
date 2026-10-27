---
output: 
  bookdown::html_document2:
     toc: true
---




# Gene annotation plots.





## Introduction
Plot gene annotation data.





## Load necessary libraries.


``` r
# install.packages(c("tidyverse"))
library(tidyverse)
```

```
## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
## ✔ dplyr     1.1.4     ✔ readr     2.1.5
## ✔ forcats   1.0.0     ✔ stringr   1.5.1
## ✔ ggplot2   4.0.0     ✔ tibble    3.2.1
## ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
## ✔ purrr     1.0.4     
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ purrr::%||%()   masks base::%||%()
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
```

``` r
library(gplots)
```

```
## 
## ---------------------
## gplots 3.3.0 loaded:
##   * Use citation('gplots') for citation info.
##   * Homepage: https://talgalili.github.io/gplots/
##   * Report issues: https://github.com/talgalili/gplots/issues
##   * Ask questions: https://stackoverflow.com/questions/tagged/gplots
##   * Suppress this message with: suppressPackageStartupMessages(library(gplots))
## ---------------------
## 
## 
## Attaching package: 'gplots'
## 
## The following object is masked from 'package:stats':
## 
##     lowess
```

``` r
library(stringr)
```



## Import datasets.


``` r
source(here::here("config/R/config.R"))# config file with global variables
source(here::here("config/R/themes.R"))# config file with themes
```

Mapping between PROKKA, MMseqs2, samples, contigs, and
HTseq quantification datasets:


``` r
gene.annotation.df <- readRDS(file = gene.annotation.df.fname.rds)
```

KOfamScan output (KEGG Ortholog classification):


``` r
kofamscan.df <- readRDS(file = kofamscan.df.fname.rds)
```

dbCAN output (CAZyme classification):


``` r
dbcan.df <- readRDS(file = dbcan.df.fname.rds)
```

Representative gene lengths:


``` r
representative.gene_lengths <- readRDS(file = representative.gene_lengths.fname.rds)
```

BLASTP classification of CAZymes:


``` r
cazymes.blast.df <- read.table(file = cazymes.blast.df.fname,
                               sep = "\t",header = F,fill=T)%>%
  as_tibble()

colnames(cazymes.blast.df)<-c("qacc","sacc", "qseqid","sseqid", 
                              "pident", "length", "mismatch", "gapopen", 
                              "qstart", "qend", "sstart", "send",
                              "evalue", "bitscore", "sframe","uniparc",
                              "caz_db","Kingdom","Species")
```

Clean the dataset by removing the "fasta" suffix, removing qseqid and sseqid,
renaming qacc, and renaming some taxa:


``` r
cazymes.blast.df<-cazymes.blast.df%>%
  mutate(caz_db=gsub(".*fasta\\|","",caz_db))%>%
  select(-qseqid,-sseqid)%>%
  rename("locus_tag_cluster"="qacc")%>%
  mutate(Kingdom=ifelse(Kingdom=="","unmapped",Kingdom),
         Species=ifelse(Kingdom=="unmapped","unmapped",Species))
head(cazymes.blast.df)%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## locus_tag_cluster   sacc            pident   length   mismatch   gapopen   qstart   qend   sstart   send   evalue   bitscore   sframe  uniparc      caz_db   Kingdom    Species                              
## ------------------  ------------  --------  -------  ---------  --------  -------  -----  -------  -----  -------  ---------  -------  -----------  -------  ---------  -------------------------------------
## BPLGMAPI_00136      seq01436884     25.417      480        327        13        3    461       11    480        0        181        1  ACQ54708.1   GT4      Bacteria   Clostridium botulinum Ba4 str. 657   
## BPLGMAPI_00138      seq02318547     34.979      243        155         2        3    245       39    278        0        158        1  ULO03512.1   GT2      Bacteria   Campylobacter sp. RM12651            
## BPLGMAPI_00168      seq00968406    100.000      484          0         0        1    484        1    484        0       1012        1  QUT33056.1   GH57     Bacteria   Bacteroides xylanisolvens CL11T00C03 
## BPLGMAPI_00169      seq00250542    100.000      421          0         0        1    421        1    421        0        882        1  UVP26469.1   GT4      Bacteria   Bacteroides xylanisolvens BFG-514    
## BPLGMAPI_00170      seq00626770    100.000      645          0         0        1    645        1    645        0       1355        1  QUT27492.1   GH133    Bacteria   Bacteroides xylanisolvens CL11T00C41 
## BPLGMAPI_00353      seq01585385     44.444      657        313        12        1    613        1    649        0        528        1  CCO12633.1   GH84     Bacteria   Carnobacterium maltaromaticum LMA28
```

3764752 CDS:


``` r
gene.annotation.df%>%
  filter(gene_type=="CDS")%>%
  distinct(locus_tag, .keep_all =T )%>%
  nrow
```

```
## [1] 3764752
```

Min, max, mean gene lengths


``` r
gene.annotation.df%>%
  filter(gene_type=="CDS")%>%
  distinct(locus_tag_cluster, .keep_all =T )%>%
  select(gene_length)%>%
  summary()%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
##       gene_length    
## ---  ----------------
##      Min.   :  102.0 
##      1st Qu.:  252.0 
##      Median :  405.0 
##      Mean   :  582.5 
##      3rd Qu.:  714.0 
##      Max.   :74937.0
```



## Number of unique genes by annotation in total (kofamscan, dbcan, all).
Inspired by https://www.nature.com/articles/s42003-021-02827-2


``` r
gene.annotation.tool.list <- list("dbCAN" = unique(dbcan.df$locus_tag_cluster),
                                  "KofamScan" = unique(kofamscan.df$locus_tag_cluster),
                                  "all" = unique(gene.annotation.df$locus_tag_cluster))
gene.annotation.venn <- venn(gene.annotation.tool.list, show.plot = FALSE)
plot(gene.annotation.venn)
```

<img src="010-mag-gene-annotation-plots_files/figure-html/unnamed-chunk-16-1.png" alt="" width="672" />

``` r
attr(gene.annotation.venn, "intersections")%>% 
  summary()
```

```
##                     Length  Class  Mode     
## KofamScan:all        353652 -none- character
## dbCAN:KofamScan:all   22261 -none- character
## all                 1524557 -none- character
## dbCAN:all             22387 -none- character
```

``` r
gene.annotation.venn %>%
  head()
```

```
##         num dbCAN KofamScan all
## 000       0     0         0   0
## 001 1524557     0         0   1
## 010       0     0         1   0
## 011  353652     0         1   1
## 100       0     1         0   0
## 101   22387     1         0   1
```

Get all sets (both shared and unique to each tool)


``` r
uniq.gene.counts <- attr(gene.annotation.venn, "intersections")
```

Convert list into dataframe


``` r
uniq.gene.counts <- stack(uniq.gene.counts)%>%
  rename("locus_tag_cluster" = "values",
         "annotation_type" = "ind")%>%
  count(annotation_type, sort=TRUE)%>%
  mutate(annotation_type = gsub(":all", "", annotation_type),
         annotation_type = gsub("^all$", "None", annotation_type),
         annotation_type = gsub(":", "+", annotation_type))

head(uniq.gene.counts)%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## annotation_type          n
## ----------------  --------
## None               1524557
## KofamScan           353652
## dbCAN                22387
## dbCAN+KofamScan      22261
```



## The barplot of all genes for each annotation (4 bars).


``` r
uniq.gene.counts.bp <- uniq.gene.counts%>%
  mutate(annotation_type=factor(annotation_type,
                                levels=sort(uniq.gene.counts$annotation_type)))%>%
  ggplot(aes(x=annotation_type,y=n,fill=annotation_type))+
  geom_bar(stat="identity")+
  geom_text(aes(label=n),vjust = -0.3, nudge_y = -0.5,
            size=8)+ # add counts of genes
  scale_fill_viridis_d(option = "C")+
  labs(y="Number of unique genes",x="",
       fill="Annotation type")+
  theme_bw()+
  coord_cartesian(ylim = c(0,max(uniq.gene.counts$n)+200000),
                  expand = F)
```

Save the plot as rds object


``` r
saveRDS(uniq.gene.counts.bp, file = file.path(mag.manuscript.figures, 
                                                  "uniq.gene.counts.bp.rds"))

uniq.gene.counts.bp <- uniq.gene.counts.bp +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=20), # size of y axis ticks
        axis.title = element_text(size = 20), # size of axis names
        plot.title = element_text(size = 25), # size of plot title
        plot.caption = element_text(size=23), # size of plot caption
        legend.text = element_text(size = 20), # size of legend text
        legend.title = element_text(size = 25), # size of legend title
        legend.position = "inside", # legend inside the plot
        legend.position.inside=c(0.25,0.8),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) # legend on the left top
```

``` r
print(uniq.gene.counts.bp)
```

<img src="010-mag-gene-annotation-plots_files/figure-html/unnamed-chunk-22-1.png" alt="" width="768" />

``` r
for(image.format in image.formats){
  ggsave(filename = paste("number-of-genes-in-total", image.format, sep = "."),
         path = mag.figures,
         plot=uniq.gene.counts.bp,
         width=8, height=6,units="in",
         dpi=300,device = image.format)
}

rm(gene.annotation.tool.list)
rm(gene.annotation.venn)
rm(uniq.gene.counts)
rm(uniq.gene.counts.bp)
gc()
```

```
##            used  (Mb) gc trigger  (Mb)  max used  (Mb)
## Ncells  5867363 313.4    8845847 472.5   8845847 472.5
## Vcells 41942901 320.0  108229319 825.8 108226962 825.8
```



## Most abundant genes in total with and without filtering by prevalence.
Find non-redundant genes (locus_tag_clusters) in each sample.


``` r
non_red_genes.by_sample <- gene.annotation.df %>%
  distinct(sample_id, locus_tag_cluster)
```

Group by KOs:


``` r
non_red_genes.by_sample.ko <- non_red_genes.by_sample %>%
  left_join(kofamscan.df, by = "locus_tag_cluster")%>%
  group_by(sample_id, ko)%>%
  count(ko, name = "n_locus_tag_clusters", sort = TRUE) %>% # when you count, you only get sample_id, ko, n_locus_tag_cluster
  ungroup()%>%
  left_join(kofamscan.df[,c("ko", "ko_definition")]%>%distinct(), by = "ko")%>% # add KO definitions
  arrange(desc(n_locus_tag_clusters))%>%
  mutate(ko = ifelse(is.na(ko),"KO unavailable", ko),
         ko_definition = ifelse(is.na(ko_definition),"KO unavailable",ko_definition)) # clean NAs
```

Group by CAZymes:


``` r
non_red_genes.by_sample.caz_subclass <- non_red_genes.by_sample %>%
  left_join(dbcan.df[,c("locus_tag_cluster","caz_with_aa_coord",
                        "caz_subclass")]%>%distinct(),
            relationship = "many-to-many", 
            by = "locus_tag_cluster")%>%
  group_by(sample_id, caz_subclass)%>%
  count(caz_subclass, name = "n_locus_tag_clusters", sort = TRUE)%>%
  ungroup()%>%
  left_join(dbcan.df[,c("caz_subclass", "caz_class")]%>%distinct())%>%
  arrange(desc(n_locus_tag_clusters))
```

```
## Joining with `by = join_by(caz_subclass)`
```

Showing the first few rows:


``` r
head(non_red_genes.by_sample.ko)%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## sample_id   ko                n_locus_tag_clusters  ko_definition  
## ----------  ---------------  ---------------------  ---------------
## Y51b        KO unavailable                  275671  KO unavailable 
## Y66b        KO unavailable                  273214  KO unavailable 
## G18         KO unavailable                  267520  KO unavailable 
## H15         KO unavailable                  261886  KO unavailable 
## H3          KO unavailable                  249607  KO unavailable 
## 2D14        KO unavailable                  237742  KO unavailable
```

``` r
head(non_red_genes.by_sample.caz_subclass)%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## sample_id   caz_subclass    n_locus_tag_clusters  caz_class 
## ----------  -------------  ---------------------  ----------
## Y51b        NA                            369950  NA        
## Y66b        NA                            367548  NA        
## H15         NA                            362373  NA        
## G18         NA                            358524  NA        
## H3          NA                            338610  NA        
## 2D14        NA                            324005  NA
```

Count prevalence of KOs (the number of samples it was observed in):


``` r
ko.prevalence <- non_red_genes.by_sample.ko %>%
  distinct(sample_id, ko )%>%
  count(ko, name = "prevalence", sort = TRUE)
```

Count prevalence of CAZyme subclasses (the number of samples it was observed in):


``` r
caz_subclass.prevalence <- non_red_genes.by_sample.caz_subclass %>%
  filter(!is.na(caz_subclass))%>%
  distinct(sample_id, caz_subclass)%>%
  count(caz_subclass, name = "prevalence", sort = TRUE)
```

Showing the first few rows:


``` r
head(ko.prevalence)%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## ko        prevalence
## -------  -----------
## K00003            11
## K00005            11
## K00008            11
## K00009            11
## K00010            11
## K00012            11
```

``` r
head(caz_subclass.prevalence)%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## caz_subclass    prevalence
## -------------  -----------
## AA4                     11
## AA7                     11
## CBM0                    11
## CBM13                   11
## CBM2                    11
## CBM20                   11
```



### Plot KOs with the highest number of non-redundant genes.


``` r
top.ko.initial.boxplot <- non_red_genes.by_sample.ko %>%
  filter(!is.na(ko), ko!= "KO unavailable")%>% # we have to remove NAs because they'd skew the plot
  group_by(ko)%>%
  summarise(median_n_locus_tag_clusters = median(n_locus_tag_clusters), # calculate median and mean
            mean_n_locus_tag_clusters = mean(n_locus_tag_clusters))%>%
  arrange(desc(median_n_locus_tag_clusters))%>%
  slice_head(n = 15)%>% # get the first 15
  left_join(non_red_genes.by_sample.ko)%>% # add the mean and median to the full dataset (keep only rows of 15 KOs)
  mutate(ko = factor(ko, levels = unique(ko)),  # sort the names for plot
         ko_definition = factor(ko_definition, levels = unique(ko_definition)))%>%
  ggplot(aes(x = ko, y = n_locus_tag_clusters, fill = ko_definition))+
  geom_boxplot()+
  geom_jitter() +
  # geom_text(subset(top.ko,
  #                  locus_tag_cluster%in%top15.ko.initial$locus_tag_cluster) %>%
  #             distinct(locus_tag_cluster,.keep_all = TRUE),
  #           inherit.aes = TRUE,
  #           mapping=aes(x=locus_tag_cluster,
  #                       # y=tpm,
  #                       y = n_locus_tags,
  #                       label=rep_gene_len),
  #           vjust = -0.3, nudge_y = -0.5,
  #           size=6)+ # add counts of genes
  theme_bw()+
  scale_fill_viridis_d(option = "C",
                       alpha = 0.5)+
  labs(x="KO ID",
       # y="Median TPM",
       y="Number of unique genes per sample",
       fill="KO annotation",
       title = "Top 15 KOs by median number of unique gene sequences",
       subtitle = "Only non-redundant sequences were used")+
  theme(#plot.margin=unit(c(1,1,1,1.5), 'cm'),
        axis.text.x = element_text(size=10,angle = 45,hjust=1),
        axis.text.y = element_text(size=10), # size of y axis ticks
        axis.title = element_text(size = 10), # size of axis names
        plot.title = element_text(size = 15), # size of plot title
        plot.caption = element_text(size=13), # size of plot caption
        legend.text = element_text(size = 10), # size of legend text
        legend.title = element_text(size = 15), # size of legend title
        legend.position = "right" ,
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())# legend on the right
```

```
## Joining with `by = join_by(ko)`
```

``` r
print(top.ko.initial.boxplot)
```

<img src="010-mag-gene-annotation-plots_files/figure-html/unnamed-chunk-33-1.png" alt="" width="1056" />

``` r
for(image.format in image.formats){
  ggsave(filename = paste("top-ko-initial", image.format, sep = "."),
         path = mag.figures,
         plot=top.ko.initial.boxplot,
         width=10, height=5,units="in",
         dpi=300,device = image.format)
}
```



#### Keep KOs with prevalence >2 samples out of 11.


``` r
top.ko.filtered.boxplot <- ko.prevalence %>%
  filter(prevalence > 2)%>%
  left_join(non_red_genes.by_sample.ko, by =join_by(ko))%>%
  filter(!is.na(ko), ko_definition != "KO unavailable")%>%
  group_by(ko)%>%
  summarise(median_n_locus_tag_clusters = median(n_locus_tag_clusters), # calculate median and mean
            mean_n_locus_tag_clusters = mean(n_locus_tag_clusters))%>%
  arrange(desc(median_n_locus_tag_clusters))%>%
  slice_head(n = 15)%>% # get the first 15
  left_join(non_red_genes.by_sample.ko)%>% # add the mean and median to the full dataset (keep only rows of 15 KOs)
  mutate(ko = factor(ko, levels = unique(ko)),  # sort the names for plot
         ko_definition = factor(ko_definition, levels = unique(ko_definition)))%>%
  ggplot(aes(x = ko, y = n_locus_tag_clusters, fill = ko_definition))+
  geom_boxplot()+
  geom_jitter() +
  # geom_text(subset(top.ko,
  #                  locus_tag_cluster%in%top15.genes.filtered$locus_tag_cluster)%>%
  #             distinct(locus_tag_cluster,.keep_all = TRUE),
  #           inherit.aes = TRUE,
  #           mapping=aes(x=locus_tag_cluster,
  #                       # y=tpm,
  #                       y = n_locus_tags,
  #                       label=rep_gene_len),
  #           vjust = -0.3, nudge_y = -0.5,
  #           size=3)+ # add counts of genes
  theme_bw()+
  scale_fill_viridis_d(option = "C",
                       alpha = 0.5)+
  labs(x="KO ID",
       # y="Median TPM",
       y="Number of unique genes per sample",
       title = "Top 15 KOs by median number of unique gene sequences",
       fill="KO annotation",
       subtitle = "Only prevalent KOs found in at least 2 out of 11 samples are shown here")+
  theme(
    # plot.margin=unit(c(0.2,0.2,0.2,0.5), 'cm'),
    axis.text.x = element_text(size=10,angle = 45,hjust=1),
    axis.text.y = element_text(size=10), # size of y axis ticks
    axis.title = element_text(size = 10), # size of axis names
    plot.title = element_text(size = 15), # size of plot title
    plot.caption = element_text(size=13), # size of plot caption
    legend.text = element_text(size = 10), # size of legend text
    legend.title = element_text(size = 15), # size of legend title
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.position = "right" )
```

```
## Joining with `by = join_by(ko)`
```

``` r
print(top.ko.filtered.boxplot)
```

<img src="010-mag-gene-annotation-plots_files/figure-html/unnamed-chunk-36-1.png" alt="" width="960" />

``` r
for(image.format in image.formats){
  ggsave(filename = paste("top-ko-filtered-prevalence", image.format,sep = "."),
         path = mag.figures,
         plot=top.ko.filtered.boxplot,
         width=10, height=5,units="in",
         dpi=300,device = image.format)
}
```



### Plot CAZyme subclasses with the highest number of non-redundant genes.


``` r
top.cazyme.initial.boxplot <-  non_red_genes.by_sample.caz_subclass %>%
  filter(!is.na(caz_subclass ))%>% # we have to remove NAs because they'd skew the plot
  group_by(caz_subclass )%>%
  summarise(median_n_locus_tag_clusters = median(n_locus_tag_clusters), # calculate median and mean
            mean_n_locus_tag_clusters = mean(n_locus_tag_clusters))%>%
  arrange(desc(median_n_locus_tag_clusters))%>%
  slice_head(n = 15)%>% # get the first 15
  left_join(non_red_genes.by_sample.caz_subclass)%>% # add the mean and median to the full dataset (keep only rows of 15 KOs)
  mutate(caz_subclass  = factor(caz_subclass , levels = unique(caz_subclass )),  # sort the names for plot
         caz_class = factor(caz_class, levels = unique(caz_class)))%>%
  ggplot(aes(x = caz_subclass , y = n_locus_tag_clusters, fill = caz_class))+
  geom_boxplot()+
  geom_jitter(cex = 0.8)+
  labs(x="CAZyme subclass",
       # y="Median TPM",
       y="Number of unique genes per sample",
       title = "Top 15 CAZyme subclasses by median number of unique gene sequences",
       fill="CAZyme class",
       subtitle = "Some genes were assigned to multiple subclasses because they have multiple modules;\ntherefore, some subclasses may seem more abundant")+
  theme_bw()+
  scale_fill_viridis_d(option = "C",
                       alpha = 0.5)+
  # coord_cartesian(expand = FALSE)+
  theme(axis.text.x = element_text(size=12,angle=45,vjust = 1,hjust = 1),
                      axis.text.y = element_text(size=10), # size of y axis ticks
                      axis.title = element_text(size = 10), # size of axis names
                      plot.title = element_text(size = 15), # size of plot title
                      plot.caption = element_text(size=13), # size of plot caption
                      legend.text = element_text(size = 10), # size of legend text
                      legend.title = element_text(size = 15), # size of legend title
                      panel.grid.minor = element_blank(),
                      panel.grid.major = element_blank(),
                      legend.position = "right") # legend on the right
```

```
## Joining with `by = join_by(caz_subclass)`
```

``` r
print(top.cazyme.initial.boxplot)
```

<img src="010-mag-gene-annotation-plots_files/figure-html/unnamed-chunk-39-1.png" alt="" width="960" />

``` r
for(image.format in image.formats){
  ggsave(filename =paste( "top-cazymes-initial", image.format, sep = "."),
         path = mag.figures,
         plot=top.cazyme.initial.boxplot,
         width=10, height=5,units="in",
         dpi=300,device = image.format)
}
```



#### Keep CAZyme subclasses with prevalence >2 samples out of 11.


``` r
top15.caz_subclasses.filtered <- caz_subclass.prevalence %>%
  filter(prevalence > 2)%>%
  left_join(non_red_genes.by_sample.caz_subclass, by =join_by(caz_subclass))%>%
  filter(!is.na(caz_subclass))%>%
  group_by(caz_subclass)%>%
  summarise(median_n_locus_tag_clusters = median(n_locus_tag_clusters), # calculate median and mean
            mean_n_locus_tag_clusters = mean(n_locus_tag_clusters))%>%
  arrange(desc(median_n_locus_tag_clusters))%>%
  slice_head(n = 15)%>% # get the first 15
  left_join(non_red_genes.by_sample.caz_subclass)%>% # add the mean and median to the full dataset (keep only rows of 15 KOs)
  mutate(caz_subclass = factor(caz_subclass, levels = unique(caz_subclass)),  # sort the names for plot
         caz_class = factor(caz_class, levels = unique(caz_class)))
```

```
## Joining with `by = join_by(caz_subclass)`
```

``` r
top.cazyme.filtered.boxplot <- top15.caz_subclasses.filtered %>%
  ggplot(aes(x = caz_subclass, y = n_locus_tag_clusters, fill = caz_class))+
  geom_boxplot()+
  geom_jitter() +
  # geom_text(subset(top.ko,
  #                  locus_tag_cluster%in%top15.genes.filtered$locus_tag_cluster)%>%
  #             distinct(locus_tag_cluster,.keep_all = TRUE),
  #           inherit.aes = TRUE,
  #           mapping=aes(x=locus_tag_cluster,
  #                       # y=tpm,
  #                       y = n_locus_tags,
  #                       label=rep_gene_len),
  #           vjust = -0.3, nudge_y = -0.5,
  #           size=3)+ # add counts of genes
  theme_bw()+
  scale_fill_viridis_d(option = "C",
                       alpha = 0.5)+
  labs(x="CAZyme subclass",
       # y="Median TPM",
       y="Number of unique genes per sample",
       title = "Top 15 CAZyme subclasses by median number of unique gene sequences",
       fill="CAZyme class",
       subtitle = "Only prevalent subclasses found in at least 2 out of 11 samples are shown here")+
  theme(
    # plot.margin=unit(c(0.2,0.2,0.2,0.5), 'cm'),
    axis.text.x = element_text(size=10,angle = 45,hjust=1),
    axis.text.y = element_text(size=10), # size of y axis ticks
    axis.title = element_text(size = 10), # size of axis names
    plot.title = element_text(size = 15), # size of plot title
    plot.caption = element_text(size=13), # size of plot caption
    legend.text = element_text(size = 10), # size of legend text
    legend.title = element_text(size = 15), # size of legend title
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    legend.position = "right" )
```

``` r
print(top.cazyme.filtered.boxplot)
```

<img src="010-mag-gene-annotation-plots_files/figure-html/unnamed-chunk-42-1.png" alt="" width="960" />

``` r
for(image.format in image.formats){
  ggsave(filename = paste("top-cazymes-filtered-prevalence", image.format,sep = "."),
         path = mag.figures,
         plot=top.cazyme.filtered.boxplot,
         width=10, height=5,units="in",
         dpi=300,device = image.format)
}
```



## Summary of the BLASTP results.


``` r
cazymes.blast.df%>%
  count(Kingdom,sort=TRUE)
```

```
## # A tibble: 6 × 2
##   Kingdom          n
##   <chr>        <int>
## 1 Bacteria     40958
## 2 Eukaryota     3010
## 3 Archaea        191
## 4 unmapped       169
## 5 Viruses        166
## 6 unclassified   154
```

Find species with the highest number of cazymes


``` r
top.cazyme_species<-cazymes.blast.df%>%
  add_count(Species,sort=TRUE)%>%
  distinct(Species,Kingdom,n)%>%
  rename(n_cazymes=n)
head(top.cazyme_species)
```

```
## # A tibble: 6 × 3
##   Species                                            Kingdom   n_cazymes
##   <chr>                                              <chr>         <int>
## 1 Lactiplantibacillus plantarum subsp. plantarum P-8 Bacteria       1160
## 2 Carnobacterium maltaromaticum LMA28                Bacteria        816
## 3 Arabidopsis arenosa                                Eukaryota       696
## 4 Parabacteroides goldsteinii BFG-241                Bacteria        632
## 5 Brassica rapa                                      Eukaryota       609
## 6 Muribaculum gordoncarteri TLL-A4                   Bacteria        573
```

``` r
top.cazyme_species.fname <- file.path(mag.tables, "cazymes-blastp-taxonomy-top-species.tsv")
if(!file.exists(top.cazyme_species.fname)){
  write.table(top.cazyme_species,
              file=top.cazyme_species.fname,
              sep="\t",col.names = TRUE,row.names = FALSE)
  
}
```

Check how many species the top 15 CAZyme subclasses were found in BLASTP 
(showing the first ten):


``` r
cazymes.blast.df%>%
  filter(caz_db %in% top15.caz_subclasses.filtered$caz_subclass)%>%
  rename("caz_subclass" = "caz_db")%>%
  group_by(caz_subclass)%>%
  count(Species, name = "n_cazymes", sort = TRUE)%>%
  ungroup() %>%
  slice_head(n = 10)%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## caz_subclass   Species                                                                    n_cazymes
## -------------  ------------------------------------------------------------------------  ----------
## GT2            Lactiplantibacillus plantarum subsp. plantarum P-8                               579
## GT2            Salmonella enterica subsp. enterica serovar Enteritidis str. EC20121176          359
## GH1            Carnobacterium maltaromaticum LMA28                                              340
## GH28           Arabidopsis arenosa                                                              274
## GT4            Lactiplantibacillus plantarum subsp. plantarum P-8                               252
## GH28           Vitis vinifera                                                                   237
## GT4            Ardenticatena Cfx-K                                                              211
## GH23           Lactiplantibacillus plantarum subsp. plantarum P-8                               196
## GT2            Sphingomonas bisphenolicum AO1                                                   196
## CE4            Volvariella volvacea V23                                                         190
```

Only for GH subclass:


``` r
cazymes.blast.df%>%
  rename("caz_subclass" = "caz_db")%>%
  left_join(non_red_genes.by_sample.caz_subclass[,c("caz_subclass", "caz_class")]%>%distinct(),
                          by = join_by("caz_subclass"))%>%
  filter(caz_class=="GH")%>%
  group_by(caz_subclass)%>%
  add_count(Species)%>%
  ungroup()%>%
  # filter(caz_subclass%in%c("GH45","GH11","GH44"))%>%
  distinct(caz_subclass,Species,Kingdom)%>%
  count(Species,sort=TRUE)%>%
  head()%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## Species                                  n
## -------------------------------------  ---
## Bacteroides faecium CBA7301             50
## Bacteroides cellulosilyticus BFG-250    40
## Parabacteroides goldsteinii BFG-241     38
## Alistipes senegalensis JC50             34
## Bacteroides sp. BFG-257                 34
## Lachnospiraceae bacterium CE91-St56     31
```

Save data.


``` r
rm(top.ko.initial.boxplot)
rm(top.ko.filtered.boxplot)
rm(top.cazyme.initial.boxplot)
rm(top.cazyme.filtered.boxplot)
rm(top15.caz_subclasses.filtered)
rm(cazymes.blast.df)
rm(representative.gene_lengths)
gc()
```

```
##            used  (Mb) gc trigger  (Mb)  max used  (Mb)
## Ncells  5897986 315.0    8845847 472.5   8845847 472.5
## Vcells 46110591 351.8  108229319 825.8 108226962 825.8
```

``` r
# save.image("./output/rdafiles/gene-annotation-workspace.RData")
```



## Find core CAZymes and KOs and compare with iMGMC.


``` r
# load("./output/rdafiles/gene-annotation-workspace.RData")
core.ko <- ko.prevalence %>%
  left_join(kofamscan.df)%>%
  rename("n_samples"="prevalence")%>%
  filter(!is.na(ko),
         n_samples>5)%>%
  filter(!is.na(ko))%>%
  distinct(ko,ko_definition)
```

```
## Joining with `by = join_by(ko)`
```

``` r
head(core.ko)%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## ko       ko_definition                                                                           
## -------  ----------------------------------------------------------------------------------------
## K00003   homoserine dehydrogenase [EC:1.1.1.3]                                                   
## K00005   glycerol dehydrogenase [EC:1.1.1.6]                                                     
## K00008   L-iditol 2-dehydrogenase [EC:1.1.1.14]                                                  
## K00009   mannitol-1-phosphate 5-dehydrogenase [EC:1.1.1.17]                                      
## K00010   myo-inositol 2-dehydrogenase / D-chiro-inositol 1-dehydrogenase [EC:1.1.1.18 1.1.1.369] 
## K00012   UDPglucose 6-dehydrogenase [EC:1.1.1.22]
```

``` r
nrow(core.ko)
```

```
## [1] 3810
```

``` r
core.caz_subclasses <- caz_subclass.prevalence%>%
  left_join(dbcan.df, by = "caz_subclass")%>%
  rename("n_samples"="prevalence")%>%
  filter(n_samples>5)%>%
  select(caz_subclass,caz_class)%>%
  filter(!is.na(caz_subclass))%>%
  distinct(caz_subclass,.keep_all = TRUE)
head(core.caz_subclasses)%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## caz_subclass   caz_class 
## -------------  ----------
## AA4            AA        
## AA7            AA        
## CBM0           CBM       
## CBM13          CBM       
## CBM2           CBM       
## CBM20          CBM
```

``` r
nrow(core.caz_subclasses)
```

```
## [1] 258
```

``` r
core.ko.fname <- file.path(mag.tables, "core-kos-by-n_samples.tsv")
core.caz_subclasses.fname <- file.path(mag.tables, "core-caz_subclasses-by-n_samples.tsv")
if(!file.exists(core.ko.fname)){
  write.table(core.ko,
              file=core.ko.fname,
              sep = "\t",
              row.names = F)
}
if(!file.exists(core.caz_subclasses.fname)){
  write.table(core.caz_subclasses,
              file=core.caz_subclasses.fname,
              sep = "\t",
              row.names = F)
}

core.ko.names<-core.ko%>%
  distinct(ko)%>%
  pull(ko)
core.caz_subclasses.names<-core.caz_subclasses%>%
  distinct(caz_subclass)%>%
  pull(caz_subclass)
```



### Compare core CAZymes.
Mouse core gut data are found here:
https://zenodo.org/records/3631711/files/iMGMC_map_functionality.tar.gz?download=1
https://pmc.ncbi.nlm.nih.gov/articles/PMC7059117/
Read the iMGMC core CAZyme dataset:


``` r
imgmc.caz<-read.table(imgmc.cazy.fname,
                      header = FALSE, sep = "\t")%>%
  as_tibble()
```

Tidy up the data:


``` r
imgmc.caz.vector<-imgmc.caz%>%
  rename("caz"="V2")%>%
  distinct(caz)%>%
  mutate(caz_split=sub("[A-Z0-9\\._]*\\|","",caz),
         caz_split=gsub("_[0-9]+","",caz_split),
         caz_split=gsub("\\+","|",caz_split))%>%
  pull(caz_split)
```

Get unique subclasses:


``` r
imgmc.caz.str<-paste(imgmc.caz.vector,collapse = "|")
imgmc.caz.all<-str_split(imgmc.caz.vector, "\\|")
imgmc.caz.all<-unlist(imgmc.caz.all,use.names = F)
imgmc.caz.all<-unique(imgmc.caz.all)
```

Find which core NMR subclasses are absent in core mouse data.


``` r
nmr.imgmc.core.caz.venn <- venn(list ("NMR" = core.caz_subclasses.names,
                                      "iMGMC" = imgmc.caz.all))
```

<img src="010-mag-gene-annotation-plots_files/figure-html/unnamed-chunk-55-1.png" alt="" width="672" />

``` r
nmr.core.caz.vector <- attr(nmr.imgmc.core.caz.venn, "intersections")$NMR
length(nmr.core.caz.vector)
```

```
## [1] 50
```

Save the core NMR subclasses.


``` r
nmr.core.caz.df <- nmr.core.caz.vector%>%
  as_tibble()%>%
  rename("caz_subclass"="value")%>%
  inner_join(core.caz_subclasses)
```

```
## Joining with `by = join_by(caz_subclass)`
```

``` r
head(nmr.core.caz.df)
```

```
## # A tibble: 6 × 2
##   caz_subclass caz_class
##   <chr>        <chr>    
## 1 CBM88        CBM      
## 2 CBM91        CBM      
## 3 CBM93        CBM      
## 4 CBM97        CBM      
## 5 CE17         CE       
## 6 CE19         CE
```

``` r
nmr.core.caz.df.fname <- file.path(mag.tables, "core-caz_subclasses-nmr_specific.tsv")
if(!file.exists(nmr.core.caz.df.fname)){
  write.table(nmr.core.caz.df,
              file = nmr.core.caz.df.fname,
              sep = "\t",
              row.names = F)
  
}
```



### Compare core KOs.


``` r
imgmc.ko<-read.table(imgmc.ko.fname,
                     header = FALSE, sep = "\t")%>%
  as_tibble()

imgmc.ko.vector<-imgmc.ko%>%
  rename("ko"="V2")%>%
  distinct(ko)%>%
  pull(ko)
```

Find the core NMR-specific KOs:


``` r
nmr.imgmc.core.ko.venn <- venn(list ("NMR" = core.ko.names,
                                     "iMGMC" = imgmc.ko.vector))
```

<img src="010-mag-gene-annotation-plots_files/figure-html/unnamed-chunk-59-1.png" alt="" width="672" />

``` r
nmr.core.ko.vector<-attr(nmr.imgmc.core.ko.venn, "intersections")$NMR
length(nmr.core.ko.vector)
```

```
## [1] 555
```

Save the dataset:


``` r
nmr.core.ko.df<-core.ko%>%
  distinct(ko,.keep_all = TRUE)%>%
  filter(ko%in%nmr.core.ko.vector)%>%
  select(ko,ko_definition)
head(nmr.core.ko.df)
```

```
## # A tibble: 6 × 2
##   ko     ko_definition                                                          
##   <chr>  <chr>                                                                  
## 1 K00078 dihydrodiol dehydrogenase / D-xylose 1-dehydrogenase (NADP) [EC:1.3.1.…
## 2 K00188 2-oxoisovalerate ferredoxin oxidoreductase delta subunit [EC:1.2.7.7]  
## 3 K00200 formylmethanofuran dehydrogenase subunit A [EC:1.2.7.12]               
## 4 K00201 formylmethanofuran dehydrogenase subunit B [EC:1.2.7.12]               
## 5 K00203 formylmethanofuran dehydrogenase subunit D [EC:1.2.7.12]               
## 6 K00319 methylenetetrahydromethanopterin dehydrogenase [EC:1.5.98.1]
```

``` r
nmr.core.ko.df.fname <- file.path(mag.tables, "core-kos-nmr_specific.tsv")
if(!file.exists(nmr.core.ko.df.fname)){
  write.table(nmr.core.ko.df,
              file = nmr.core.ko.df.fname,
              sep = "\t",
              row.names = FALSE)
  
}
```



## Match gene annotation with taxonomic classification from GTDB-tk.
Also look for a few taxa that were found in 16S data.


``` r
seqkit.all.df <- readRDS(file = seqkit.all.df.fname.rds)
gtdbtk.coverm.drep <- readRDS(file = gtdbtk.coverm.drep.fname.rds)
gtdbtk.coverm.metabat2 <- readRDS(file = gtdbtk.coverm.metabat2.fname.rds)
blastn.coverm.megahit <- readRDS(file = blastn.coverm.megahit.fname.rds)
gtdbtk.taxonomy <- readRDS(gtdbtk.taxonomy.clean.fname.rds)
```

Find p-251-o5 in BLASTN classification of contigs: 14 contigs


``` r
blastn.coverm.megahit%>%
  filter(grepl("p-251-o5",classification))%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## sample_id   contig_id             tpm  silva_ref_id       perc_identity   alignment_length   bitscore  classification                                                                               
## ----------  ------------  -----------  ----------------  --------------  -----------------  ---------  ---------------------------------------------------------------------------------------------
## 2D14        k141_82886     176.739960  EU776984.1.1390           99.496               1390       2527  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;p-251-o5;Incertae Sedis;uncultured bacterium 
## G14         k141_160344     41.026450  EU776984.1.1390           99.841               1259       2313  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;p-251-o5;Incertae Sedis;uncultured bacterium 
## H15         k141_234799      9.134913  EU463139.1.1400           90.553                561        734  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;p-251-o5;Incertae Sedis;uncultured bacterium 
## H15         k141_108011    623.028400  EU776984.1.1390          100.000                726       1341  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;p-251-o5;Incertae Sedis;uncultured bacterium 
## H21         k141_106317     17.535395  EU776984.1.1390           99.712               1390       2543  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;p-251-o5;Incertae Sedis;uncultured bacterium 
## H3          k141_110992    298.576140  EU776984.1.1390           99.877                812       1493  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;p-251-o5;Incertae Sedis;uncultured bacterium 
## H3          k141_45731     353.773500  EU463139.1.1400           99.722                719       1317  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;p-251-o5;Incertae Sedis;uncultured bacterium 
## H4          k141_185402     50.323930  EU776984.1.1390           99.501               1002       1821  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;p-251-o5;Incertae Sedis;uncultured bacterium 
## H4          k141_181003    205.511730  EU463139.1.1400           99.622                529        966  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;p-251-o5;Incertae Sedis;uncultured bacterium 
## O15         k141_110667    380.984560  EU776984.1.1390           93.037                517        750  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;p-251-o5;Incertae Sedis;uncultured bacterium 
## O15         k141_246933    838.466700  EU463139.1.1400           98.571                560        989  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;p-251-o5;Incertae Sedis;uncultured bacterium 
## O15         k141_130609    128.941960  EU463139.1.1400           92.721               1154       1646  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;p-251-o5;Incertae Sedis;uncultured bacterium 
## Y51b        k141_275532    277.199770  EU463139.1.1400           99.857                700       1288  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;p-251-o5;Incertae Sedis;uncultured bacterium 
## Y66b        k141_244384    336.055820  EU776984.1.1390           99.640               1390       2538  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;p-251-o5;Incertae Sedis;uncultured bacterium
```

Prevotellaceae UCG-003: 32 contigs


``` r
blastn.coverm.megahit%>%
  filter(grepl("Prevotellaceae UCG-003",classification))%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## sample_id   contig_id             tpm  silva_ref_id       perc_identity   alignment_length   bitscore  classification                                                                                                   
## ----------  ------------  -----------  ----------------  --------------  -----------------  ---------  -----------------------------------------------------------------------------------------------------------------
## 2D10        k141_186065    193.503900  EU463094.1.1400           99.284               1117       2019  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotellaceae UCG-003;uncultured bacterium       
## 2D14        k141_28550     192.559890  EU777021.1.1401           98.405                627       1101  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotellaceae UCG-003;uncultured bacterium       
## 2D14        k141_237558     90.355470  EU850578.1.1452           90.388                593        778  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotellaceae UCG-003;uncultured rumen bacterium 
## 2D14        k141_283071    261.068500  EU777021.1.1401           99.172                725       1306  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotellaceae UCG-003;uncultured bacterium       
## 2D14        k141_71293      31.835516  AB824397.1.1521           91.246               1051       1430  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotellaceae UCG-003;uncultured bacterium       
## G14         k141_32626     296.800350  EU777021.1.1401           99.001               1401       2508  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotellaceae UCG-003;uncultured bacterium       
## G14         k141_11103      37.553677  AB824397.1.1521           90.092               1302       1685  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotellaceae UCG-003;uncultured bacterium       
## G14         k141_203929      4.582249  AB555283.1.1490           90.724               1423       1890  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotellaceae UCG-003;uncultured rumen bacterium 
## G18         k141_317875    201.604630  AB821632.1.1521           96.789                654       1092  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotellaceae UCG-003;uncultured bacterium       
## G18         k141_72113     285.300570  EU463094.1.1400           99.007               1007       1805  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotellaceae UCG-003;uncultured bacterium       
## H15         k141_28147     154.333620  EU777021.1.1401           99.511                613       1116  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotellaceae UCG-003;uncultured bacterium       
## H15         k141_214830     21.270971  AB821632.1.1521           90.529               1531       2008  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotellaceae UCG-003;uncultured bacterium       
## H15         k141_246170     75.723830  AB824397.1.1521           91.728                677        939  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotellaceae UCG-003;uncultured bacterium       
## H15         k141_153696    141.432160  EU463094.1.1400           94.945               1266       1977  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotellaceae UCG-003;uncultured bacterium       
## H21         k141_9720      386.549320  EU777021.1.1401           99.459                555       1009  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotellaceae UCG-003;uncultured bacterium       
## H21         k141_9722      195.177490  EU463094.1.1400           90.110               1001       1279  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotellaceae UCG-003;uncultured bacterium       
## H21         k141_151573     40.777860  AB824397.1.1521           91.728                677        939  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotellaceae UCG-003;uncultured bacterium       
## H21         k141_67575      24.039602  EU842573.1.1505           92.224               1196       1688  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotellaceae UCG-003;uncultured bacterium       
## H21         k141_132121     42.348760  KJ874004.1.1320           92.161                523        737  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotellaceae UCG-003;uncultured bacterium       
## H3          k141_220988    141.488950  EU463094.1.1400           96.159               1406       2287  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotellaceae UCG-003;uncultured bacterium       
## H3          k141_24360      20.713750  AB821632.1.1521           90.594               1531       2013  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotellaceae UCG-003;uncultured bacterium       
## H3          k141_156189    229.004210  EU777021.1.1401           99.438                534        970  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotellaceae UCG-003;uncultured bacterium       
## H4          k141_191576     26.339390  AB615139.1.1491           90.017                591        765  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotellaceae UCG-003;uncultured rumen bacterium 
## H4          k141_223439      1.634680  AB615139.1.1491           93.802                726       1092  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotellaceae UCG-003;uncultured rumen bacterium 
## H4          k141_93847     253.524230  EU463094.1.1400           97.502               1401       2392  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotellaceae UCG-003;uncultured bacterium       
## O15         k141_172131    361.914600  EU777021.1.1401           99.492                591       1075  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotellaceae UCG-003;uncultured bacterium       
## O15         k141_111578    125.135000  EU777021.1.1401           93.665               1405       2093  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotellaceae UCG-003;uncultured bacterium       
## Y51b        k141_35801      94.221950  EU463094.1.1400           99.143               1400       2519  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotellaceae UCG-003;uncultured bacterium       
## Y51b        k141_208982     15.751941  AB824397.1.1521           92.672                655        942  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotellaceae UCG-003;uncultured bacterium       
## Y66b        k141_343895     32.357647  AB824397.1.1521           92.240                683        966  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotellaceae UCG-003;uncultured bacterium       
## Y66b        k141_88138      17.835249  EU871423.1.1411           90.145                761        985  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotellaceae UCG-003;uncultured rumen bacterium 
## Y66b        k141_56393      98.050100  EU463094.1.1400           99.000               1400       2508  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotellaceae UCG-003;uncultured bacterium
```

Treponema and Fibrobacter in representative MAGs (GTDB-tk): 8 and 1 MAGs


``` r
gtdbtk.taxonomy%>%
  filter(grepl("Treponema|Fibrobacter",Genus))%>%
  arrange(Genus)%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## sample_id    bin_id  secondary_cluster   Kingdom    Phylum           Class           Order             Family             Genus         Species           
## ----------  -------  ------------------  ---------  ---------------  --------------  ----------------  -----------------  ------------  ------------------
## O15              37  108_1               Bacteria   Fibrobacterota   Fibrobacteria   Fibrobacterales   Fibrobacteraceae   Fibrobacter   Fibrobacter Genus 
## 2D10             76  105_0               Bacteria   Spirochaetota    Spirochaetia    Treponematales    Treponemataceae    Treponema_D   Treponema_D Genus 
## 2D14             78  106_0               Bacteria   Spirochaetota    Spirochaetia    Treponematales    Treponemataceae    Treponema_D   Treponema_D Genus 
## G18              27  113_1               Bacteria   Spirochaetota    Spirochaetia    Treponematales    Treponemataceae    Treponema_D   Treponema_D Genus 
## H15              23  112_1               Bacteria   Spirochaetota    Spirochaetia    Treponematales    Treponemataceae    Treponema_D   Treponema_D Genus 
## H3               10  111_1               Bacteria   Spirochaetota    Spirochaetia    Treponematales    Treponemataceae    Treponema_D   Treponema_D Genus 
## H4              108  119_1               Bacteria   Spirochaetota    Spirochaetia    Treponematales    Treponemataceae    Treponema_D   Treponema_D Genus 
## O15              14  118_1               Bacteria   Spirochaetota    Spirochaetia    Treponematales    Treponemataceae    Treponema_D   Treponema_D Genus 
## Y66b             56  120_0               Bacteria   Spirochaetota    Spirochaetia    Treponematales    Treponemataceae    Treponema_D   Treponema_D Genus
```

Erysipelotrichaceae and Prevotellaceae: 5 MAGs


``` r
gtdbtk.taxonomy%>%
  filter(grepl("Erysipelotrichaceae|Prevotellaceae",Family))%>%
  arrange(Genus)%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## sample_id    bin_id  secondary_cluster   Kingdom    Phylum      Class     Order                Family                Genus            Species              
## ----------  -------  ------------------  ---------  ----------  --------  -------------------  --------------------  ---------------  ---------------------
## O15              22  18_0                Bacteria   Bacillota   Bacilli   Erysipelotrichales   Erysipelotrichaceae   Allobaculum      Allobaculum Genus    
## 2D14             91  107_0               Bacteria   Bacillota   Bacilli   Erysipelotrichales   Erysipelotrichaceae   Colivicinus      Colivicinus Genus    
## H21             115  109_1               Bacteria   Bacillota   Bacilli   Erysipelotrichales   Erysipelotrichaceae   Colivicinus      Colivicinus Genus    
## O15              10  17_0                Bacteria   Bacillota   Bacilli   Erysipelotrichales   Erysipelotrichaceae   Faecalibaculum   Faecalibaculum Genus 
## G18              52  2_0                 Bacteria   Bacillota   Bacilli   Erysipelotrichales   Erysipelotrichaceae   Ileibacterium    Ileibacterium Genus
```



### Which individual MAGs have the highest number of unique CAZyme subclasses.
Join samples and contigs from classified bins with taxonomy classification,
then with gene annotation for each contig


``` r
mag.classification.per_contig <- seqkit.all.df %>%
  filter(is_representative_bin == TRUE)%>% # Only representative 128 bins
  mutate(bin_id = as.integer(bin_id))%>%
  left_join(gtdbtk.taxonomy, by = join_by(sample_id, bin_id, secondary_cluster))
```

Combine distinct locus_tag_clusters with gene annotation per contig and add
CAZyme data:


``` r
mag.classification.caz <- mag.classification.per_contig %>%
  left_join(gene.annotation.df,by = join_by(sample_id, contig_id))%>%
  # distinct(secondary_cluster, locus_tag_cluster)%>% 
  left_join(dbcan.df,by = join_by(locus_tag_cluster),
            relationship = "many-to-many") # Add CAZymes
```

Count unique CAZymes and CAZyme subclasses in each MAG:


``` r
mag.classification.caz %>% 
  group_by(secondary_cluster)%>%
  summarise(n_uniq_caz_in_mag=n_distinct(caz_with_aa_coord ),
            n_uniq_caz_subclass_in_mag=n_distinct(caz_subclass),
            n_uniq_cas_class_in_mag=n_distinct(caz_class))%>%
  arrange(desc(n_uniq_caz_subclass_in_mag))%>%
  left_join(gtdbtk.taxonomy, by = join_by(secondary_cluster))%>%
  ggplot(aes (x = Phylum, y = n_uniq_caz_subclass_in_mag))+
  geom_boxplot()
```

<img src="010-mag-gene-annotation-plots_files/figure-html/unnamed-chunk-70-1.png" alt="" width="672" />

Find CAZyme subclasses shared by representative MAGs (presence/absence):


``` r
mag_family.caz <- mag.classification.caz %>%
  filter(!is.na(caz_subclass))%>%
  distinct(Family, caz_subclass)
mag_genus.caz <- mag.classification.caz %>%
  filter(!is.na(caz_subclass))%>%
  distinct(Genus, caz_subclass)
```

First, convert the data into wide format: the values are genus names, so
we can collapse them into strings.


``` r
caz.membership.genus <- mag_genus.caz%>%
  pivot_wider(names_from = Genus,
              values_from = Genus,
              values_fill = NA)%>%
  as.data.frame()%>%
  column_to_rownames("caz_subclass")%>% 
  select(order(colnames(.)))%>%
  unite(data = ., 
        col ="membership", 
        sep = ";",
        remove = FALSE,
        na.rm = TRUE) %>% # Concatenate genus names into a string, remove NA so we don't have a ton of semicolons
  mutate(across(!membership, ~replace(., !is.na(.), "1")), # Replace non-NA values with 1
         across(!membership, ~replace(., is.na(.), "0")))%>%# Replace non-NA values with 0
  unite(col ="membership_vector",  
        -"membership",
        sep = "",
        remove = TRUE,
        na.rm = TRUE) %>%# Concatenate 1 and 0 except the column with concatenated names
  mutate(n_taxa = str_count(membership_vector, pattern = "1")) %>% # count "1" characters
  select(membership, membership_vector, n_taxa)%>%
  as_tibble(rownames = "caz_subclass")
#same code 
caz.membership.family <- mag_family.caz%>%
  pivot_wider(names_from = Family,
              values_from = Family,
              values_fill = NA)%>%
  as.data.frame()%>%
  column_to_rownames("caz_subclass")%>% 
  select(order(colnames(.)))%>%
  unite(data = ., 
        col ="membership", 
        sep = ";",
        remove = FALSE,
        na.rm = TRUE) %>% # Concatenate genus names into a string, remove NA so we don't have a ton of semicolons
  mutate(across(!membership, ~replace(., !is.na(.), "1")), # Replace non-NA values with 1
         across(!membership, ~replace(., is.na(.), "0")))%>%# Replace non-NA values with 0
  unite(col ="membership_vector",  
        -"membership",
        sep = "",
        remove = TRUE,
        na.rm = TRUE) %>%# Concatenate 1 and 0 except the column with concatenated names
  mutate(n_taxa = str_count(membership_vector, pattern = "1")) %>% # count "1" characters
  select(membership, membership_vector, n_taxa)%>%
  as_tibble(rownames = "caz_subclass")

# At MAG level
caz.membership.mag <- mag.classification.caz %>%
  filter(!is.na(caz_subclass))%>%
  distinct(secondary_cluster, caz_subclass)%>%
  pivot_wider(names_from = secondary_cluster,
              values_from = secondary_cluster,
              values_fill = NA)%>%
  as.data.frame()%>%
  column_to_rownames("caz_subclass")%>% 
  select(order(colnames(.)))%>%
  unite(data = ., 
        col ="membership", 
        sep = ";",
        remove = FALSE,
        na.rm = TRUE) %>% # Concatenate IDs into a string, remove NA so we don't have a ton of semicolons
  mutate(across(!membership, ~replace(., !is.na(.), "1")), # Replace non-NA values with 1
         across(!membership, ~replace(., is.na(.), "0")))%>%# Replace non-NA values with 0
  unite(col ="membership_vector",  
        -"membership",
        sep = "",
        remove = TRUE,
        na.rm = TRUE) %>%# Concatenate 1 and 0 except the column with concatenated names
  mutate(n_taxa = str_count(membership_vector, pattern = "1")) %>% # count "1" characters
  select(membership, membership_vector, n_taxa)%>%
  as_tibble(rownames = "caz_subclass")

treponema.ids<-gtdbtk.taxonomy%>%
  filter(grepl("Treponema", Genus))%>%
  pull(secondary_cluster)
# CAZ subclasses found in ALL Treponema MAGs
caz.membership.mag %>%
  mutate(n_matches = str_count(membership, pattern = paste(treponema.ids,collapse = "|"))  )%>%
  # filter(n_matches ==length(treponema.ids)) # ALL mags
  filter(n_matches >0) # ANY mags
```

```
## # A tibble: 90 × 5
##    caz_subclass membership                    membership_vector n_taxa n_matches
##    <chr>        <chr>                         <chr>              <int>     <int>
##  1 GH43         102_0;107_0;108_1;11_1;110_0… 0000100001101111…     52         4
##  2 GT2          1_1;10_1;100_1;101_1;102_0;1… 1111111111111111…    128         8
##  3 PL10         120_0;32_1;35_1;36_0;59_1;60… 0000000000000000…      6         1
##  4 CE8          101_1;103_1;104_1;108_1;110_… 0001011000100100…     92         1
##  5 GH51         108_1;11_1;111_1;118_1;119_1… 0000000000101010…     40         3
##  6 GH23         1_1;10_1;100_1;101_1;103_1;1… 1111011111111111…    119         8
##  7 GH73         1_1;10_1;101_1;105_0;106_0;1… 1101000111011011…     68         8
##  8 GT4          1_1;10_1;100_1;101_1;102_0;1… 1111111111111111…    128         8
##  9 GH1          1_1;10_1;100_1;101_1;102_0;1… 1111111111111111…    121         8
## 10 GH3          100_1;102_0;103_1;104_1;105_… 0010111111101111…    107         8
## # ℹ 80 more rows
```

``` r
caz.membership.genus %>%
  filter(grepl("Treponema", membership))
```

```
## # A tibble: 90 × 4
##    caz_subclass membership                              membership_vector n_taxa
##    <chr>        <chr>                                   <chr>              <int>
##  1 GH43         Acetatifactor;Acutalibacteraceae Famil… 1010101010001000…     28
##  2 GT2          Acetatifactor;Acutalibacter;Acutalibac… 1111111111111111…     70
##  3 PL10         Bacteroides;Cryptobacteroides;Prevotel… 0000000000001000…      4
##  4 CE8          Acetatifactor;Acutalibacter;Acutalibac… 1111101110011101…     50
##  5 GH51         Acetatifactor;Acutalibacter;Acutalibac… 1110010010001000…     19
##  6 GH23         Acetatifactor;Acutalibacter;Acutalibac… 1111111111001111…     65
##  7 GH73         Acetatifactor;Anaerostipes;Anaerovorac… 1000001100011100…     38
##  8 GT4          Acetatifactor;Acutalibacter;Acutalibac… 1111111111111111…     70
##  9 GH1          Acetatifactor;Acutalibacter;Acutalibac… 1111111111111101…     68
## 10 GH3          Acetatifactor;Acutalibacter;Acutalibac… 1111111111101110…     59
## # ℹ 80 more rows
```

``` r
  # filter(membership == "Treponema_D")

caz.membership.mag %>%
  filter(n_taxa == nrow(gtdbtk.taxonomy)) # 2 subclasses
```

```
## # A tibble: 2 × 4
##   caz_subclass membership                               membership_vector n_taxa
##   <chr>        <chr>                                    <chr>              <int>
## 1 GT2          1_1;10_1;100_1;101_1;102_0;103_1;104_1;… 1111111111111111…    128
## 2 GT4          1_1;10_1;100_1;101_1;102_0;103_1;104_1;… 1111111111111111…    128
```

``` r
caz.membership.mag%>%
  filter(n_taxa ==1)# 28
```

```
## # A tibble: 28 × 4
##    caz_subclass membership membership_vector                              n_taxa
##    <chr>        <chr>      <chr>                                           <int>
##  1 GH111        26_0       000000000000000000000000000000000000000000000…      1
##  2 PL22         68_0       000000000000000000000000000000000000000000000…      1
##  3 GH119        42_1       000000000000000000000000000000000000000000000…      1
##  4 GH180        68_0       000000000000000000000000000000000000000000000…      1
##  5 GH54         32_1       000000000000000000000000000000000000000000000…      1
##  6 CBM85        83_0       000000000000000000000000000000000000000000000…      1
##  7 2.4.1.12     126_0      000000000000000000000000000000010000000000000…      1
##  8 GT99         24_1       000000000000000000000000000000000000000000001…      1
##  9 GT114        112_1      000000000000000100000000000000000000000000000…      1
## 10 3.2.1.21     56_1       000000000000000000000000000000000000000000000…      1
## # ℹ 18 more rows
```

``` r
caz.membership.mag%>%
  filter(n_taxa ==1)%>%
  left_join(gtdbtk.taxonomy[,c("secondary_cluster","Family","Genus","Species")],
            by = join_by("membership" == "secondary_cluster"))%>%
  count(Genus)
```

```
## # A tibble: 17 × 2
##    Genus                       n
##    <chr>                   <int>
##  1 Acetatifactor               2
##  2 Adiutrix                    1
##  3 Aristaeellaceae Family      2
##  4 Atopobiaceae Family         1
##  5 Bacteroides                 2
##  6 CAG-196                     1
##  7 Cryptobacteroides           1
##  8 Eubacterium_G               4
##  9 Fibrobacter                 4
## 10 Fimivicinus                 1
## 11 JAGBWK01                    1
## 12 Oliverpabstia               1
## 13 Parasutterella              1
## 14 Phascolarctobacterium       1
## 15 Prevotella                  1
## 16 Selenomonadaceae Family     2
## 17 Treponema_D                 2
```

``` r
caz.membership.genus %>%
  filter(n_taxa == n_distinct(gtdbtk.taxonomy$Genus)) # 2 subclasses, same as above
```

```
## # A tibble: 2 × 4
##   caz_subclass membership                               membership_vector n_taxa
##   <chr>        <chr>                                    <chr>              <int>
## 1 GT2          Acetatifactor;Acutalibacter;Acutalibact… 1111111111111111…     70
## 2 GT4          Acetatifactor;Acutalibacter;Acutalibact… 1111111111111111…     70
```

``` r
caz.membership.genus %>%
  filter(n_taxa ==1)  # 42 because we're looking at genera
```

```
## # A tibble: 42 × 4
##    caz_subclass membership             membership_vector                  n_taxa
##    <chr>        <chr>                  <chr>                               <int>
##  1 CBM72        Prevotella             000000000000000000000000000000000…      1
##  2 GH111        Fimivicinus            000000000000000000000000000000001…      1
##  3 PL22         Aristaeellaceae Family 000000001000000000000000000000000…      1
##  4 GH37         Cryptobacteroides      000000000000000000000001000000000…      1
##  5 GH119        JAGBWK01               000000000000000000000000000000000…      1
##  6 GH180        Aristaeellaceae Family 000000001000000000000000000000000…      1
##  7 GH54         Cryptobacteroides      000000000000000000000001000000000…      1
##  8 PL15         Bacteroides            000000000000100000000000000000000…      1
##  9 CBM85        Oliverpabstia          000000000000000000000000000000000…      1
## 10 2.4.1.12     Parasutterella         000000000000000000000000000000000…      1
## # ℹ 32 more rows
```

``` r
caz.membership.genus %>%
  left_join(dbcan.df[,c("caz_subclass", "caz_class")]%>%distinct())%>%
  group_by(caz_class)%>%
  mutate(mean_n_taxa = mean(n_taxa))%>%
  ungroup()%>%
  arrange(desc(mean_n_taxa))%>%
  mutate(caz_class = factor(caz_class, levels = unique(caz_class)))%>%
  ggplot(aes(x= caz_class, y = n_taxa)) +
  geom_boxplot()
```

```
## Joining with `by = join_by(caz_subclass)`
```

<img src="010-mag-gene-annotation-plots_files/figure-html/unnamed-chunk-72-1.png" alt="" width="672" />

``` r
caz.membership.genus %>%
  left_join(dbcan.df[,c("caz_subclass", "caz_class")]%>%distinct())%>%
  arrange(desc(n_taxa))%>%
  mutate(caz_subclass = factor(caz_subclass, levels = caz_subclass))%>%
  ggplot(aes(x= caz_subclass, y = n_taxa)) +
  geom_bar(stat = "identity")+
  theme(axis.text.x = element_text(angle = 45))
```

```
## Joining with `by = join_by(caz_subclass)`
```

<img src="010-mag-gene-annotation-plots_files/figure-html/unnamed-chunk-72-2.png" alt="" width="672" />

``` r
# caz.membership.mag %>%
caz.membership.genus %>%
  filter(caz_subclass %in% nmr.core.caz.vector)%>%
  left_join(dbcan.df[,c("caz_subclass", "caz_class")]%>%distinct())%>%
  arrange(desc(n_taxa))%>%
  mutate(caz_subclass = factor(caz_subclass, levels = caz_subclass))%>%
  ggplot(aes(x= caz_subclass, y = n_taxa)) +
  geom_bar(stat = "identity")+
  theme(axis.text.x = element_text(angle = 45))
```

```
## Joining with `by = join_by(caz_subclass)`
```

<img src="010-mag-gene-annotation-plots_files/figure-html/unnamed-chunk-72-3.png" alt="" width="672" />

``` r
# CORE CAZY ARE NOT ALWAYS WITHIN MAGS
dbcan.df%>%
  filter(caz_subclass == "GH151")%>%
  left_join(gene.annotation.df)%>%
  select(sample_id,contig_id)%>%
  left_join(seqkit.all.df)
```

```
## Joining with `by = join_by(locus_tag_cluster)`
## Joining with `by = join_by(sample_id, contig_id)`
```

```
## # A tibble: 16 × 8
##    sample_id contig_id   bin_id      contig_length    GC    AT secondary_cluster
##    <chr>     <chr>       <chr>               <int> <dbl> <dbl> <chr>            
##  1 2D14      k141_169621 34                   4953  42.3  57.7 low_quality_mag  
##  2 2D14      k141_144722 103                  6958  58.4  41.6 low_quality_mag  
##  3 G18       k141_332438 97                  29207  62.3  37.7 low_quality_mag  
##  4 G18       k141_332438 97                  29207  62.3  37.7 low_quality_mag  
##  5 G18       k141_165808 unassembled          9329  58.5  41.5 unassembled      
##  6 H21       k141_91704  52                   3246  48.9  51.1 low_quality_mag  
##  7 H3        k141_12190  23                   8327  49.3  50.7 80_1             
##  8 H4        k141_120195 101                  8359  49.4  50.6 80_1             
##  9 Y51b      k141_311715 23                  11416  49.7  50.3 80_1             
## 10 Y66b      k141_242225 89                   8215  49.2  50.8 80_1             
## 11 Y66b      k141_24878  unassembled          2483  41.5  58.5 unassembled      
## 12 H15       k141_128781 unassembled          1227  47.3  52.7 unassembled      
## 13 H3        k141_206756 unassembled           830  39.2  60.8 unassembled      
## 14 O15       k141_171284 29                   2974  44.8  55.2 low_quality_mag  
## 15 G18       k141_332115 unassembled          1953  37.5  62.5 unassembled      
## 16 2D14      k141_100074 unassembled          1332  43.3  56.7 unassembled      
## # ℹ 1 more variable: is_representative_bin <chr>
```

``` r
# FIND ALL GENES IN ALL MAGS, NOT ONLY DREP
seqkit.all.df%>%
  filter(!secondary_cluster %in% c("unassembled", "low_quality_mag"))%>%
  left_join(gene.annotation.df)%>%
  distinct(locus_tag_cluster)
```

```
## Joining with `by = join_by(sample_id, contig_id)`
```

```
## # A tibble: 279,247 × 1
##    locus_tag_cluster
##    <chr>            
##  1 PLACIMPO_186011  
##  2 LBKEGPMD_360542  
##  3 GPCOLHEP_00317   
##  4 GPCOLHEP_00318   
##  5 GPCOLHEP_00319   
##  6 GPCOLHEP_00320   
##  7 GPCOLHEP_00322   
##  8 GPCOLHEP_00323   
##  9 GPCOLHEP_00324   
## 10 GPCOLHEP_00325   
## # ℹ 279,237 more rows
```

``` r
mag_family.caz %>%  
  filter( grepl("Treponemataceae", Family)) %>%
  nrow()
```

```
## [1] 90
```

``` r
mag_family.caz %>%  
  filter( grepl("Erysipelotrichaceae",Family))%>%
  nrow()
```

```
## [1] 64
```

``` r
mag_family.caz %>%
  filter( grepl("Treponemataceae", Family) & grepl("Erysipelotrichaceae",Family))%>%
  nrow()
```

```
## [1] 0
```

For Treponema and Fibrobacter:


``` r
mag.classification.caz %>%
  filter(!is.na(caz_subclass))%>%
  # left_join(gtdbtk.taxonomy, by = join_by(secondary_cluster))%>%
  group_by(Genus)%>%
  summarise(n_uniq_caz_in_taxon=n_distinct(caz_with_aa_coord),
            n_uniq_caz_subclass_in_taxon=n_distinct(caz_subclass),
            n_uniq_cas_class_in_taxon=n_distinct(caz_class))%>%
  arrange(desc(n_uniq_caz_subclass_in_taxon))%>%
  filter(grepl("Treponema|Fibrobacter",Genus))%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## Genus          n_uniq_caz_in_taxon   n_uniq_caz_subclass_in_taxon   n_uniq_cas_class_in_taxon
## ------------  --------------------  -----------------------------  --------------------------
## Treponema_D                    483                             90                           5
## Fibrobacter                    135                             65                           5
```

For Prevotellaceae (named Bacteroidaceae in GTDB-tk) and Erysipelotrichaceae:


``` r
mag.classification.caz %>%
  filter(!is.na(caz_subclass))%>%
  # left_join(gtdbtk.taxonomy, by = join_by(secondary_cluster))%>%
  group_by(Family)%>%
  summarise(n_uniq_caz_in_taxon=n_distinct(caz_with_aa_coord),
            n_uniq_caz_subclass_in_taxon=n_distinct(caz_subclass),
            n_uniq_cas_class_in_taxon=n_distinct(caz_class))%>%
  arrange(desc(n_uniq_caz_subclass_in_taxon))%>%
  filter(grepl("Erysipelotrichaceae|Bacteroidaceae",Family))%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## Family                 n_uniq_caz_in_taxon   n_uniq_caz_subclass_in_taxon   n_uniq_cas_class_in_taxon
## --------------------  --------------------  -----------------------------  --------------------------
## Bacteroidaceae                        2344                            187                           6
## Erysipelotrichaceae                    323                             64                           5
```

Break down each MAG by CAZyme subclasses: get the highest number of CAZymes 
for each subclass in each MAG.


``` r
caz_per_family <- mag.classification.caz %>%
  filter(!is.na(caz_subclass))%>%
  group_by(Family,caz_subclass)%>%
  summarise(n_uniq_caz_in_taxon=n_distinct(caz_with_aa_coord))%>%
  arrange(Family,desc(n_uniq_caz_in_taxon))%>%
  ungroup()
```

```
## `summarise()` has grouped output by 'Family'. You can override using the
## `.groups` argument.
```

``` r
caz_per_genus <- mag.classification.caz%>%
  filter(!is.na(caz_subclass))%>%
  group_by(Genus,caz_subclass)%>%
  summarise(n_uniq_caz_in_taxon=n_distinct(caz_with_aa_coord))%>%
  arrange(Genus,desc(n_uniq_caz_in_taxon))%>%
  ungroup()
```

```
## `summarise()` has grouped output by 'Genus'. You can override using the
## `.groups` argument.
```

``` r
treponema.caz <- caz_per_family %>%
  filter(grepl("Treponema",Family))
head(treponema.caz)%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## Family            caz_subclass    n_uniq_caz_in_taxon
## ----------------  -------------  --------------------
## Treponemataceae   GH13                             63
## Treponemataceae   GT4                              50
## Treponemataceae   GT2                              37
## Treponemataceae   GH3                              25
## Treponemataceae   GH5                              17
## Treponemataceae   GH23                             16
```

``` r
fibrobacter.caz <- caz_per_genus %>%
  filter(grepl("Fibrobacter",Genus))
head(fibrobacter.caz)%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## Genus         caz_subclass    n_uniq_caz_in_taxon
## ------------  -------------  --------------------
## Fibrobacter   GT4                              17
## Fibrobacter   GT2                              15
## Fibrobacter   GH5                              11
## Fibrobacter   CBM6                              8
## Fibrobacter   GH13                              6
## Fibrobacter   GH8                               5
```

``` r
erysipelotrichaceae.caz <- caz_per_family%>%
  filter(grepl("Erysipelotrichaceae",Family))
head(erysipelotrichaceae.caz)%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## Family                caz_subclass    n_uniq_caz_in_taxon
## --------------------  -------------  --------------------
## Erysipelotrichaceae   GT2                              38
## Erysipelotrichaceae   GH13                             35
## Erysipelotrichaceae   GT4                              34
## Erysipelotrichaceae   GH1                              30
## Erysipelotrichaceae   CE1                              12
## Erysipelotrichaceae   GH77                             12
```

``` r
prevotella.caz <- caz_per_genus%>%
  filter(grepl("Prevotella",Genus))
head(prevotella.caz)%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## Genus        caz_subclass    n_uniq_caz_in_taxon
## -----------  -------------  --------------------
## Prevotella   GT2                              75
## Prevotella   GH43                             64
## Prevotella   GH2                              58
## Prevotella   GT4                              55
## Prevotella   GH3                              39
## Prevotella   GH28                             36
```

How many CAZyme subclasses are found in Treponema but not Fibrobacter: 48


``` r
length(setdiff(treponema.caz$caz_subclass,fibrobacter.caz$caz_subclass))
```

```
## [1] 48
```

Reverse: 23


``` r
length(setdiff(fibrobacter.caz$caz_subclass,treponema.caz$caz_subclass))
```

```
## [1] 23
```

Alternative way of showing this:


``` r
caz_subclasses.venn <- venn(list ("Treponema" = unique(treponema.caz$caz_subclass),
                                     "Fibrobacter" = unique(fibrobacter.caz$caz_subclass),
                                  "Erysipelotrichaceae" = unique(erysipelotrichaceae.caz$caz_subclass),
                                  "Prevotella" = unique(prevotella.caz$caz_subclass)))
```

<img src="010-mag-gene-annotation-plots_files/figure-html/unnamed-chunk-78-1.png" alt="" width="672" />

Get all sets (both shared and unique to each taxon)


``` r
caz_subclasses.by_set <- attr(caz_subclasses.venn, "intersections")
```

Convert list into dataframe


``` r
caz_subclasses.by_set <- stack(caz_subclasses.by_set)%>%
  rename("caz_subclass" = "values",
         "Taxon" = "ind")

caz_subclasses.by_set %>%
  filter(grepl("Treponema", Taxon) & ! grepl("Fibrobacter", Taxon))%>%
  nrow()
```

```
## [1] 48
```

``` r
caz_subclasses.by_set %>%
  filter(grepl("Fibrobacter", Taxon)  & ! grepl("Treponema", Taxon))%>%
  nrow()
```

```
## [1] 23
```

``` r
caz_subclasses.by_set%>%
  filter(grepl("Erysipelotrichaceae", Taxon))
```

```
##    caz_subclass                                                Taxon
## 1         GH170                                  Erysipelotrichaceae
## 2           AA7                                  Erysipelotrichaceae
## 3         CBM66                                  Erysipelotrichaceae
## 4         GH110                                  Erysipelotrichaceae
## 5          CE13                                  Erysipelotrichaceae
## 6         GT116                                  Erysipelotrichaceae
## 7          GH42                        Treponema:Erysipelotrichaceae
## 8           CE9                        Treponema:Erysipelotrichaceae
## 9           GH4                        Treponema:Erysipelotrichaceae
## 10        GH179                       Erysipelotrichaceae:Prevotella
## 11         GH29                       Erysipelotrichaceae:Prevotella
## 12        CBM32                       Erysipelotrichaceae:Prevotella
## 13        GH163                       Erysipelotrichaceae:Prevotella
## 14        GH177                       Erysipelotrichaceae:Prevotella
## 15         GH20                       Erysipelotrichaceae:Prevotella
## 16         GH89                       Erysipelotrichaceae:Prevotella
## 17         GH92                       Erysipelotrichaceae:Prevotella
## 18         GT14                       Erysipelotrichaceae:Prevotella
## 19         GT47                       Erysipelotrichaceae:Prevotella
## 20         GT58                       Erysipelotrichaceae:Prevotella
## 21        CBM51                      Fibrobacter:Erysipelotrichaceae
## 22          GH0                      Fibrobacter:Erysipelotrichaceae
## 23        GH130             Treponema:Erysipelotrichaceae:Prevotella
## 24         GH31             Treponema:Erysipelotrichaceae:Prevotella
## 25          CE4             Treponema:Erysipelotrichaceae:Prevotella
## 26         GH32             Treponema:Erysipelotrichaceae:Prevotella
## 27         GH36             Treponema:Erysipelotrichaceae:Prevotella
## 28         GH73             Treponema:Erysipelotrichaceae:Prevotella
## 29        GT111             Treponema:Erysipelotrichaceae:Prevotella
## 30         CE20             Treponema:Erysipelotrichaceae:Prevotella
## 31        GH172             Treponema:Erysipelotrichaceae:Prevotella
## 32         GH38             Treponema:Erysipelotrichaceae:Prevotella
## 33         GH78             Treponema:Erysipelotrichaceae:Prevotella
## 34         GH95             Treponema:Erysipelotrichaceae:Prevotella
## 35        CBM48            Treponema:Fibrobacter:Erysipelotrichaceae
## 36        GH152            Treponema:Fibrobacter:Erysipelotrichaceae
## 37          CE1           Fibrobacter:Erysipelotrichaceae:Prevotella
## 38          CE6           Fibrobacter:Erysipelotrichaceae:Prevotella
## 39         GH84           Fibrobacter:Erysipelotrichaceae:Prevotella
## 40         GT30           Fibrobacter:Erysipelotrichaceae:Prevotella
## 41         GH13 Treponema:Fibrobacter:Erysipelotrichaceae:Prevotella
## 42          GT4 Treponema:Fibrobacter:Erysipelotrichaceae:Prevotella
## 43          GT2 Treponema:Fibrobacter:Erysipelotrichaceae:Prevotella
## 44          GH3 Treponema:Fibrobacter:Erysipelotrichaceae:Prevotella
## 45          GH5 Treponema:Fibrobacter:Erysipelotrichaceae:Prevotella
## 46         GH23 Treponema:Fibrobacter:Erysipelotrichaceae:Prevotella
## 47         GT28 Treponema:Fibrobacter:Erysipelotrichaceae:Prevotella
## 48         GH77 Treponema:Fibrobacter:Erysipelotrichaceae:Prevotella
## 49          GT5 Treponema:Fibrobacter:Erysipelotrichaceae:Prevotella
## 50         GH43 Treponema:Fibrobacter:Erysipelotrichaceae:Prevotella
## 51         GT51 Treponema:Fibrobacter:Erysipelotrichaceae:Prevotella
## 52          GT8 Treponema:Fibrobacter:Erysipelotrichaceae:Prevotella
## 53        CBM50 Treponema:Fibrobacter:Erysipelotrichaceae:Prevotella
## 54         GH25 Treponema:Fibrobacter:Erysipelotrichaceae:Prevotella
## 55          GH2 Treponema:Fibrobacter:Erysipelotrichaceae:Prevotella
## 56         GH51 Treponema:Fibrobacter:Erysipelotrichaceae:Prevotella
## 57          GH1 Treponema:Fibrobacter:Erysipelotrichaceae:Prevotella
## 58         GT35 Treponema:Fibrobacter:Erysipelotrichaceae:Prevotella
## 59         CE12 Treponema:Fibrobacter:Erysipelotrichaceae:Prevotella
## 60         GH28 Treponema:Fibrobacter:Erysipelotrichaceae:Prevotella
## 61          GT0 Treponema:Fibrobacter:Erysipelotrichaceae:Prevotella
## 62          GT1 Treponema:Fibrobacter:Erysipelotrichaceae:Prevotella
## 63        GT106 Treponema:Fibrobacter:Erysipelotrichaceae:Prevotella
## 64         GT83 Treponema:Fibrobacter:Erysipelotrichaceae:Prevotella
```



## Summary statistics of all MAGs including gene annotation.
Calculate genome size, contig number, GC%, CDS number, 
CheckM2 completeness and contamination


``` r
metabat2.date_time<-"20250612_13_37_47"
checkm2.high_quality_mags<-read.table(checkm2.high_quality_mags.fname,
                                      sep="\t",header = T)%>%
  as_tibble()%>%
  rename("genome"="Name")%>%
  mutate(genome=gsub(paste0(metabat2.date_time,"_"),"",genome))%>%
  separate_wider_delim(genome,delim = "_bin.",names = c("sample","bin_id"))%>%
  mutate(bin_id=as.integer(bin_id))

checkm2.all_mags<-read.table(checkm2.all_mags.fname,
                                      sep="\t",header = T)%>%
  as_tibble()%>%
  rename("genome"="Name")%>%
  mutate(genome=gsub(paste0(metabat2.date_time,"_"),"",genome))%>%
  separate_wider_delim(genome,delim = "_bin.",names = c("sample","bin_id"))%>%
  mutate(bin_id=as.integer(bin_id))

mag.stats.with_annotation.df <- seqkit.all.df %>%
  filter(is_representative_bin !="unassembled" ) %>%
  select(sample_id,contig_id,bin_id,contig_length,secondary_cluster,GC)%>%
  filter(!is.na(secondary_cluster))%>%
  group_by(sample_id, bin_id)%>%
  mutate(genome_size = sum(contig_length),
         mag_gc=mean(GC),
         n_contig = n_distinct(contig_id))%>%
  ungroup()%>%
  full_join(gene.annotation.df[,c("locus_tag_cluster","locus_tag","sample_id","contig_id",
                                  "gene_length","gene_type")],
            by = join_by(sample_id, contig_id))%>%
  filter(!is.na(secondary_cluster))%>% # remove genes that mapped to unassembled contigs
  group_by(sample_id, bin_id)%>%
  mutate(n_cds = n_distinct(locus_tag, na.rm = TRUE)) %>% # add gene annotations
  ungroup()%>%
  distinct(sample_id, bin_id, .keep_all = TRUE)%>%
  select(-contig_length, -contig_id, -GC, -locus_tag_cluster, -locus_tag,
         -gene_length, -gene_type)%>%
  mutate(bin_id = as.integer(bin_id))%>%
  left_join(checkm2.all_mags, #checkm2.high_quality_mags, 
            by = join_by("sample_id" == "sample", bin_id))%>%
  select(#-Completeness_Model_Used, -Translation_Table_Used,-Additional_Notes,
         -Genome_Size,-GC_Content, -Total_Contigs)%>%
  left_join(subset(gtdbtk.taxonomy, select = -c(sample_id, bin_id)),
            by = "secondary_cluster")%>%
  arrange(secondary_cluster)
head(mag.stats.with_annotation.df)
```

```
## # A tibble: 6 × 24
##   sample_id bin_id secondary_cluster genome_size mag_gc n_contig n_cds
##   <chr>      <int> <chr>                   <int>  <dbl>    <int> <int>
## 1 G14           92 100_1                 1212152   47.0       24  1082
## 2 Y66b          38 100_1                 1138829   47.6      112  1016
## 3 G18           28 101_1                 1728111   34.8       51  1348
## 4 H15           71 101_1                 1679416   34.3       91  1313
## 5 Y51b          64 101_1                 1740550   34.6       65  1345
## 6 H15           94 102_0                 1516872   27.7       55  1270
## # ℹ 17 more variables: Completeness <dbl>, Contamination <dbl>,
## #   Completeness_Model_Used <chr>, Translation_Table_Used <int>,
## #   Coding_Density <dbl>, Contig_N50 <int>, Average_Gene_Length <dbl>,
## #   Total_Coding_Sequences <int>, Max_Contig_Length <int>,
## #   Additional_Notes <chr>, Kingdom <chr>, Phylum <chr>, Class <chr>,
## #   Order <chr>, Family <chr>, Genus <chr>, Species <chr>
```

``` r
nrow(mag.stats.with_annotation.df)
```

```
## [1] 1341
```

``` r
mag.stats.with_annotation.df.fname <- file.path(mag.tables, "mag-stats-with-n_cds-new-all_mag.tsv")
if(!file.exists(mag.stats.with_annotation.df.fname)){
  write.table(mag.stats.with_annotation.df,
              file= mag.stats.with_annotation.df.fname,
              sep = "\t", row.names = FALSE)
  
}
sessionInfo()
```

```
## R version 4.4.3 (2025-02-28 ucrt)
## Platform: x86_64-w64-mingw32/x64
## Running under: Windows 11 x64 (build 22631)
## 
## Matrix products: default
## 
## 
## locale:
## [1] LC_COLLATE=English_United States.utf8 
## [2] LC_CTYPE=English_United States.utf8   
## [3] LC_MONETARY=English_United States.utf8
## [4] LC_NUMERIC=C                          
## [5] LC_TIME=English_United States.utf8    
## 
## time zone: Asia/Tokyo
## tzcode source: internal
## 
## attached base packages:
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] gplots_3.3.0    lubridate_1.9.4 forcats_1.0.0   stringr_1.5.1  
##  [5] dplyr_1.1.4     purrr_1.0.4     readr_2.1.5     tidyr_1.3.1    
##  [9] tibble_3.2.1    ggplot2_4.0.0   tidyverse_2.0.0
## 
## loaded via a namespace (and not attached):
##  [1] utf8_1.2.4         sass_0.4.10        generics_0.1.4     xml2_1.3.8        
##  [5] bitops_1.0-9       KernSmooth_2.23-26 gtools_3.9.5       stringi_1.8.4     
##  [9] hms_1.1.4          digest_0.6.37      magrittr_2.0.3     caTools_1.18.3    
## [13] evaluate_1.0.5     grid_4.4.3         timechange_0.3.0   RColorBrewer_1.1-3
## [17] bookdown_0.46      fastmap_1.2.0      rprojroot_2.1.1    jsonlite_2.0.0    
## [21] ggtext_0.1.2       viridisLite_0.4.3  scales_1.4.0       textshaping_1.0.3 
## [25] jquerylib_0.1.4    cli_3.6.4          rlang_1.1.5        withr_3.0.2       
## [29] cachem_1.1.0       yaml_2.3.12        otel_0.2.0         tools_4.4.3       
## [33] tzdb_0.5.0         here_1.0.2         vctrs_0.6.5        R6_2.6.1          
## [37] lifecycle_1.0.5    ragg_1.5.0         pkgconfig_2.0.3    pillar_1.11.1     
## [41] bslib_0.10.0       gtable_0.3.6       Rcpp_1.0.14        glue_1.8.0        
## [45] systemfonts_1.2.3  xfun_0.56          tidyselect_1.2.1   rstudioapi_0.18.0 
## [49] knitr_1.51         farver_2.1.2       htmltools_0.5.8.1  labeling_0.4.3    
## [53] rmarkdown_2.30     compiler_4.4.3     S7_0.2.0           gridtext_0.1.5
```

``` r
rm(list =setdiff(ls(all.names = TRUE), c("markdown.dir")))
gc()
```

```
##           used (Mb) gc trigger   (Mb)  max used   (Mb)
## Ncells 1747581 93.4   12402401  662.4  15281418  816.2
## Vcells 7230678 55.2  149849169 1143.3 171760064 1310.5
```

