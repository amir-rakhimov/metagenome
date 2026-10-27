---
output: 
  bookdown::html_document2:
     toc: true
---





# Analysing MAG data.





## Introduction
In this script, we will explore the MAG dataset using tidyverse and phyloseq.






## Load necessary libraries.


``` r
# install.packages(c("tidyverse"))
# BiocManager::install("phyloseq")
# install.packages(
#   "microViz",
#   repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
# )
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
library(phyloseq)
library(microViz)
```

```
## microViz version 0.12.1 - Copyright (C) 2021-2024 David Barnett
## ! Website: https://david-barnett.github.io/microViz
## ✔ Useful?  For citation details, run: `citation("microViz")`
## ✖ Silence? `suppressPackageStartupMessages(library(microViz))`
```



## Import the datasets.


``` r
source(here::here("config/R/config.R"))# config file with global variables
source(here::here("config/R/themes.R"))# config file with themes
```

Seqkit information on every contig and bin:


``` r
seqkit.all.df <- readRDS(file = seqkit.all.df.fname.rds)
```

Taxonomy classification and relative abundance of each representative bin from dRep


``` r
gtdbtk.coverm.drep <- readRDS(file = gtdbtk.coverm.drep.fname.rds)
```

Taxonomy classification and relative abundance of each bin from MetaBAT2


``` r
gtdbtk.coverm.metabat2 <- readRDS(file = gtdbtk.coverm.metabat2.fname.rds)
```

Taxonomy classification and relative abundance of each contig from BLASTN


``` r
blastn.coverm.megahit <- readRDS(file = blastn.coverm.megahit.fname.rds)
```



## Basic stats
On MEGAHIT data (contigs)


``` r
megahit.summary.stats <- seqkit.all.df %>% #seqkit.with_taxonomy%>%
  summarise(total_len=sum(contig_length),
         longest_contig_len=max(contig_length),
         shortest_contig_len=min(contig_length),
         average_len=mean(contig_length),
         n_contig=n_distinct(sample_id,contig_id))
megahit.summary.stats %>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
##   total_len   longest_contig_len   shortest_contig_len   average_len   n_contig
## -----------  -------------------  --------------------  ------------  ---------
##  4272286137               821016                   200      1370.136    3118148
```

On MetaBAT2 data (bins)


``` r
metabat2.summary.stats <- seqkit.all.df %>% # seqkit.with_taxonomy%>%
  filter(!secondary_cluster %in% c("unassembled"))%>%
  summarise(total_len=sum(contig_length),
         longest_contig_len=max(contig_length),
         shortest_contig_len=min(contig_length),
         average_len=mean(contig_length),
         n_contig=n_distinct(sample_id,contig_id)) 
metabat2.summary.stats%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
##   total_len   longest_contig_len   shortest_contig_len   average_len   n_contig
## -----------  -------------------  --------------------  ------------  ---------
##  2079334251               821016                  2500      13313.53     156182
```

On dRep data (high-quality bins)


``` r
drep.summary.stats <- seqkit.all.df %>% # seqkit.with_taxonomy%>%
  filter(is_representative_bin == TRUE)%>%
  summarise(total_len=sum(contig_length),
         longest_contig_len=max(contig_length),
         shortest_contig_len=min(contig_length),
         average_len=mean(contig_length),
         n_contig=n_distinct(sample_id,contig_id))
drep.summary.stats%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
##  total_len   longest_contig_len   shortest_contig_len   average_len   n_contig
## ----------  -------------------  --------------------  ------------  ---------
##  300172073               674940                  2500      26202.17      11456
```

Combine stats into one Summary statistics from assembly.


``` r
final.summary.stats<-rbind(megahit.summary.stats,
                           metabat2.summary.stats,
                           drep.summary.stats)%>%
  mutate(source=c("megahit","metabat2","drep"))%>%
  mutate(average_len=round(average_len))
final.summary.stats%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
##   total_len   longest_contig_len   shortest_contig_len   average_len   n_contig  source   
## -----------  -------------------  --------------------  ------------  ---------  ---------
##  4272286137               821016                   200          1370    3118148  megahit  
##  2079334251               821016                  2500         13314     156182  metabat2 
##   300172073               674940                  2500         26202      11456  drep
```

``` r
write.table(final.summary.stats,
            file = file.path(mag.tables, "mag-summary-stats.tsv"),
            sep = "\t",
            quote = FALSE,
            row.names = F)
rm(final.summary.stats)
rm(drep.summary.stats)
rm(metabat2.summary.stats)
rm(megahit.summary.stats)
```

Number of representative high-quality bins (MAGs): 128.


``` r
seqkit.all.df %>%
  filter(is_representative_bin == TRUE)%>%
  distinct(sample_id,bin_id)%>%
  nrow
```

```
## [1] 128
```

Length distribution:


``` r
seqkit.all.df %>%
  filter(is_representative_bin == TRUE)%>%
  ggplot(aes(x=contig_length))+
  geom_histogram()
```

```
## `stat_bin()` using `bins = 30`. Pick better value `binwidth`.
```

<img src="007-summary-stats-mag_files/figure-html/unnamed-chunk-17-1.png" alt="" width="672" />

Number of contigs per bin:


``` r
seqkit.all.df %>%
  filter(is_representative_bin == TRUE)%>%
  group_by(sample_id,bin_id)%>%
  summarise(n_contig_per_bin=n_distinct(contig_id))%>%
  group_by(sample_id)%>%
  summarise(mean_n_contig=mean(n_contig_per_bin))
```

```
## `summarise()` has grouped output by 'sample_id'. You can override using the
## `.groups` argument.
```

```
## # A tibble: 11 × 2
##    sample_id mean_n_contig
##    <chr>             <dbl>
##  1 2D10               74.6
##  2 2D14               94.7
##  3 G14                47.3
##  4 G18                94.9
##  5 H15                87.2
##  6 H21               107  
##  7 H3                 79.2
##  8 H4                 88.4
##  9 O15                94.3
## 10 Y51b              112. 
## 11 Y66b              108
```



## Fix the unclassified cells using tax_fix from microViz.
Need to convert the taxonomy table into a phyloseq tax_table


``` r
gtdbtk.drep.tax_fix <- gtdbtk.coverm.drep%>%
  mutate(sample_bin_cluster=paste(sample_id,"bin",bin_id,"cluster",secondary_cluster,sep = "_"))%>%
  select(classification,sample_bin_cluster)%>%
  rename(OTU=sample_bin_cluster)%>%
  separate_wider_delim(delim = ";",cols = classification,
                       names =  c("Kingdom","Phylum","Class","Order",
                                  "Family","Genus","Species"))%>%
  as.data.frame()%>%
  column_to_rownames("OTU")%>%
  as.matrix()
head(gtdbtk.drep.tax_fix)%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
##                             Kingdom    Phylum          Class          Order                Family               Genus               Species      
## --------------------------  ---------  --------------  -------------  -------------------  -------------------  ------------------  -------------
## 2D10_bin_42_cluster_69_0    Bacteria   Bacillota       Clostridia     Christensenellales   Aristaeellaceae      unclassified        unclassified 
## 2D10_bin_26_cluster_60_0    Bacteria   Bacteroidota    Bacteroidia    Bacteroidales        Bacteroidaceae       Prevotella          unclassified 
## 2D10_bin_62_cluster_28_1    Bacteria   Bacteroidota    Bacteroidia    Bacteroidales        UBA932               Cryptobacteroides   unclassified 
## 2D10_bin_39_cluster_54_1    Bacteria   Bacteroidota    Bacteroidia    Bacteroidales        Bacteroidaceae       Prevotella          unclassified 
## 2D10_bin_71_cluster_92_0    Bacteria   Bacillota       Clostridia     Oscillospirales      Acutalibacteraceae   Ruminococcoides     unclassified 
## 2D10_bin_76_cluster_105_0   Bacteria   Spirochaetota   Spirochaetia   Treponematales       Treponemataceae      Treponema_D         unclassified
```

``` r
gtdbtk.drep.tax_fix <- tax_table(gtdbtk.drep.tax_fix)
gtdbtk.drep.tax_fix <- tax_fix(gtdbtk.drep.tax_fix,
                                    unknowns = "unclassified")
gtdbtk.drep.tax_fix <- gtdbtk.drep.tax_fix%>%
  as.data.frame()%>%
  rownames_to_column(var="sample_bin_cluster")%>%
  as_tibble()%>%
  separate_wider_delim(sample_bin_cluster,delim = "_bin_",names = c("sample_id","bin_cluster"))%>%
  separate_wider_delim(bin_cluster,delim = "_cluster_",names = c("bin_id","secondary_cluster"))%>%
  mutate(bin_id=as.integer(bin_id))
head(gtdbtk.drep.tax_fix)%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## sample_id    bin_id  secondary_cluster   Kingdom    Phylum          Class          Order                Family               Genus                    Species                 
## ----------  -------  ------------------  ---------  --------------  -------------  -------------------  -------------------  -----------------------  ------------------------
## 2D10             42  69_0                Bacteria   Bacillota       Clostridia     Christensenellales   Aristaeellaceae      Aristaeellaceae Family   Aristaeellaceae Family  
## 2D10             26  60_0                Bacteria   Bacteroidota    Bacteroidia    Bacteroidales        Bacteroidaceae       Prevotella               Prevotella Genus        
## 2D10             62  28_1                Bacteria   Bacteroidota    Bacteroidia    Bacteroidales        UBA932               Cryptobacteroides        Cryptobacteroides Genus 
## 2D10             39  54_1                Bacteria   Bacteroidota    Bacteroidia    Bacteroidales        Bacteroidaceae       Prevotella               Prevotella Genus        
## 2D10             71  92_0                Bacteria   Bacillota       Clostridia     Oscillospirales      Acutalibacteraceae   Ruminococcoides          Ruminococcoides Genus   
## 2D10             76  105_0               Bacteria   Spirochaetota   Spirochaetia   Treponematales       Treponemataceae      Treponema_D              Treponema_D Genus
```

``` r
if(!file.exists(gtdbtk.taxonomy.clean.fname.rds) & !file.exists(gtdbtk.taxonomy.clean.fname.tsv)){
  saveRDS(gtdbtk.drep.tax_fix,
          file = gtdbtk.taxonomy.clean.fname.rds)
  write.table(file = gtdbtk.taxonomy.clean.fname.tsv,
              gtdbtk.drep.tax_fix,
              sep = '\t',
              row.names = F,
              col.names = T)
}
```



## Classification rate in each sample.


``` r
gtdbtk.drep.classification.rate <- gtdbtk.coverm.drep%>%
  group_by(sample_id)%>%
  summarise(sum_relative_abundance=sum(relative_abundance))%>%
  mutate(mean_classified=mean(sum_relative_abundance),
         sd_classified=sd(sum_relative_abundance))
gtdbtk.drep.classification.rate %>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## sample_id    sum_relative_abundance   mean_classified   sd_classified
## ----------  -----------------------  ----------------  --------------
## 2D10                      12.820675          8.248448        5.035587
## 2D14                       6.236914          8.248448        5.035587
## G14                        3.823877          8.248448        5.035587
## G18                       20.342875          8.248448        5.035587
## H15                        8.839494          8.248448        5.035587
## H21                        7.355897          8.248448        5.035587
## H3                         4.261603          8.248448        5.035587
## H4                        11.841826          8.248448        5.035587
## O15                        5.884497          8.248448        5.035587
## Y51b                       5.531590          8.248448        5.035587
## Y66b                       3.793681          8.248448        5.035587
```

``` r
gtdbtk.drep.classification.rate.fname <- file.path(mag.tables, "gtdbtk-drep-classification-rate.tsv")
write.table(gtdbtk.drep.classification.rate,
            file = gtdbtk.drep.classification.rate.fname,
            sep = "\t",
            col.names = T,
            row.names = F)
```



## Number of phyla, families, genera, species.
In GTDB-tk:


``` r
gtdbtk.drep.n_taxa <- gtdbtk.drep.tax_fix%>%
  summarise(PhylaPerHost=n_distinct(Phylum),
            FamiliesPerHost=n_distinct(Family),
            GeneraPerHost=n_distinct(Genus),
            SpeciesPerHost=n_distinct(secondary_cluster))
gtdbtk.drep.n_taxa %>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
##  PhylaPerHost   FamiliesPerHost   GeneraPerHost   SpeciesPerHost
## -------------  ----------------  --------------  ---------------
##            11                34              70              128
```

``` r
gtdbtk.drep.n_taxa.fname <- file.path(mag.tables, "mag-n_taxa-gtdbtk-drep.tsv")
if(!file.exists(gtdbtk.drep.n_taxa.fname)){
  write.table(gtdbtk.drep.n_taxa,
              file = gtdbtk.drep.n_taxa.fname,
              sep = "\t",
              col.names = T,
              row.names = F)
}
```

Recalculate relative abundances as relative abundance/sum(relative 
abundances with sample).


``` r
gtdbtk.coverm.drep.rescaled <- gtdbtk.coverm.drep%>%
  full_join(gtdbtk.drep.tax_fix,
            by = join_by(sample_id,bin_id,secondary_cluster))%>%
  group_by(sample_id)%>%
  mutate(old_relative_abundance=relative_abundance,
         relative_abundance=old_relative_abundance/sum(old_relative_abundance)*100)%>%
  ungroup
```

Dominant phyla (average relative abundance) and families within those phyla.
Rescaled is new, raw is old.


``` r
get_dominant_taxa_mag <- function(df, tax.rank){
  if(tax.rank == "Phylum"){
    tax.ranks.to_keep <- "Phylum"
  }else if(tax.rank == "Family"){
    tax.ranks.to_keep <- c( "Phylum", "Family")
  }else if(tax.rank == "Genus"){
    tax.ranks.to_keep <- c("Phylum", "Family", "Genus")
  }
  cols.to_keep <-c(tax.ranks.to_keep,
                   "min_old","max_old","mean_old","sd_old","n_mags",
                   "min_new","max_new","mean_new","sd_new")
  
  dominant.taxa <- df %>%
    group_by_at(tax.rank)%>%
    mutate(min_old=min(old_relative_abundance),
           max_old=max(old_relative_abundance),
           mean_old=mean(old_relative_abundance),
           sd_old=sd(old_relative_abundance),
           n_mags = n(),
           min_new=min(relative_abundance),
           max_new=max(relative_abundance),
           mean_new=mean(relative_abundance),
           sd_new=sd(relative_abundance))%>%
    ungroup%>%
    distinct(dplyr::across(all_of(cols.to_keep)))%>%
    arrange(-n_mags)
  return(dominant.taxa)
}

dominant.phyla.gtdbtk.drep <- 
  get_dominant_taxa_mag(gtdbtk.coverm.drep.rescaled,"Phylum")
  
dominant.phyla.gtdbtk.drep %>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## Phylum                  min_old     max_old    mean_old      sd_old   n_mags     min_new     max_new    mean_new      sd_new
## -------------------  ----------  ----------  ----------  ----------  -------  ----------  ----------  ----------  ----------
## Bacillota             0.1794458   5.8206810   0.5684516   0.7756941       57    1.398169   49.153576    7.067678    7.174763
## Bacteroidota          0.2089209   5.7907840   1.1421252   1.2664750       35    1.644123   61.088411   12.518065   11.573996
## Spirochaetota         0.1611320   0.5982314   0.3612197   0.1229951       10    1.607420   10.800039    5.343509    3.062228
## Actinomycetota        0.2446178   1.0590161   0.5122039   0.2764809        8    1.575957   14.396832    6.428019    3.722847
## Desulfobacterota_I    0.2877855   0.7019396   0.5107476   0.1847721        4    3.912310   11.550017    8.127407    3.797980
## Desulfobacterota      0.5666781   0.9980773   0.8053690   0.2193448        3    9.630016   13.650206   11.523779    2.020169
## Cyanobacteriota       0.3211731   0.7295291   0.5516531   0.2091986        3    7.536440   11.696957   10.052383    2.212903
## Pseudomonadota        0.2101144   0.2555319   0.2384205   0.0246916        3    1.227040    6.735724    4.500434    2.897351
## Methanobacteriota     0.7963774   0.9507980   0.8735877   0.1091919        2    4.673862   14.396899    9.535381    6.875225
## Elusimicrobiota       0.4302215   0.6607080   0.5454648   0.1629785        2    2.114851    5.579444    3.847147    2.449837
## Fibrobacterota        0.9795533   0.9795533   0.9795533          NA        1   16.646338   16.646338   16.646338          NA
```

``` r
dominant.families.gtdbtk.drep<-
  get_dominant_taxa_mag(gtdbtk.coverm.drep.rescaled,"Family")

head(dominant.families.gtdbtk.drep) %>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## Phylum          Family                  min_old     max_old    mean_old      sd_old   n_mags    min_new    max_new    mean_new      sd_new
## --------------  -------------------  ----------  ----------  ----------  ----------  -------  ---------  ---------  ----------  ----------
## Bacillota       Lachnospiraceae       0.1893376   5.8206810   0.7256053   1.3302601       17   1.398169   49.15358    9.030962   11.439526
## Bacteroidota    Muribaculaceae        0.2837276   1.8713269   0.7118574   0.4227091       14   1.644123   23.54055    8.471589    6.340747
## Bacteroidota    Bacteroidaceae        0.2898562   5.7907840   1.8700276   1.8312858       11   2.447732   61.08841   19.529354   16.718545
## Bacillota       Acutalibacteraceae    0.1794458   1.4522654   0.5198187   0.3603678       11   2.040294   13.86416    6.367607    3.852772
## Spirochaetota   Treponemataceae       0.1611320   0.5982314   0.3558038   0.1388215        8   1.607420   10.80004    5.024315    3.173004
## Bacteroidota    UBA932                0.3821162   1.2218101   0.6907716   0.3315286        7   4.257333   19.98893   10.786099    6.913874
```

``` r
dominant.genera.gtdbtk.drep<-
  get_dominant_taxa_mag(gtdbtk.coverm.drep.rescaled,"Genus")

head(dominant.genera.gtdbtk.drep)%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## Phylum          Family            Genus                       min_old     max_old    mean_old      sd_old   n_mags     min_new     max_new    mean_new      sd_new
## --------------  ----------------  -----------------------  ----------  ----------  ----------  ----------  -------  ----------  ----------  ----------  ----------
## Spirochaetota   Treponemataceae   Treponema_D               0.1611320   0.5982314   0.3558038   0.1388215        8    1.607420   10.800039    5.024315    3.173004
## Bacteroidota    UBA932            Cryptobacteroides         0.3821162   1.2218101   0.6907716   0.3315286        7    4.257333   19.988928   10.786099    6.913874
## Bacteroidota    Bacteroidaceae    Prevotella                1.1587447   5.7907840   2.9338328   1.9191905        6   10.680090   61.088411   28.766138   17.127590
## Bacteroidota    Muribaculaceae    CAG-873                   0.2837276   1.8713269   0.6809238   0.6008244        6    2.576908   14.596166    6.757199    4.985001
## Bacteroidota    Muribaculaceae    Duncaniella               0.3507090   0.9019371   0.7708235   0.2359536        5    3.967523   23.540550   10.860677    8.520790
## Bacillota       Aristaeellaceae   Aristaeellaceae Family    0.1887492   0.3399344   0.2942476   0.0711757        4    2.604795    5.450363    3.659432    1.255809
```

``` r
dominant.phyla.gtdbtk.drep.fname <- file.path(mag.tables,"mag-dominant-phyla-gtdbtk-drep.tsv")
dominant.families.gtdbtk.drep.fname <- file.path(mag.tables,"mag-dominant-families-gtdbtk-drep.tsv")
dominant.genera.gtdbtk.drep.fname <- file.path(mag.tables,"mag-dominant-genera-gtdbtk-drep.tsv")
if(!file.exists(dominant.phyla.gtdbtk.drep.fname)){
  write.table(dominant.phyla.gtdbtk.drep,
              file = dominant.phyla.gtdbtk.drep.fname,
              sep = "\t",
              col.names = T,
              row.names = F)
  
}
if(!file.exists(dominant.families.gtdbtk.drep.fname)){
  write.table(dominant.families.gtdbtk.drep,
              file=dominant.families.gtdbtk.drep.fname,
              sep = "\t",
              col.names = T,
              row.names = F)
  
}

if(!file.exists(dominant.genera.gtdbtk.drep.fname)){
  write.table(dominant.genera.gtdbtk.drep,
              file = dominant.genera.gtdbtk.drep.fname,
              sep = "\t",
              col.names = T,
              row.names = F)  
}
```



## Abundance of specific taxa.


``` r
dominant.phyla.gtdbtk.drep%>%
  filter(grepl("Bacteroidota|Spirochaetota|Desulfobacterota",Phylum))%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## Phylum                  min_old     max_old    mean_old      sd_old   n_mags    min_new    max_new    mean_new      sd_new
## -------------------  ----------  ----------  ----------  ----------  -------  ---------  ---------  ----------  ----------
## Bacteroidota          0.2089209   5.7907840   1.1421252   1.2664750       35   1.644123   61.08841   12.518065   11.573996
## Spirochaetota         0.1611320   0.5982314   0.3612197   0.1229951       10   1.607420   10.80004    5.343509    3.062228
## Desulfobacterota_I    0.2877855   0.7019396   0.5107476   0.1847721        4   3.912310   11.55002    8.127407    3.797980
## Desulfobacterota      0.5666781   0.9980773   0.8053690   0.2193448        3   9.630016   13.65021   11.523779    2.020169
```

``` r
dominant.families.gtdbtk.drep%>%
  filter(grepl("Bacteroidaceae",Family))%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## Phylum         Family              min_old    max_old   mean_old     sd_old   n_mags    min_new    max_new   mean_new     sd_new
## -------------  ---------------  ----------  ---------  ---------  ---------  -------  ---------  ---------  ---------  ---------
## Bacteroidota   Bacteroidaceae    0.2898562   5.790784   1.870028   1.831286       11   2.447732   61.08841   19.52935   16.71854
```

``` r
dominant.genera.gtdbtk.drep%>%
  filter(grepl("Treponema",Genus))%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## Phylum          Family            Genus           min_old     max_old    mean_old      sd_old   n_mags   min_new    max_new   mean_new     sd_new
## --------------  ----------------  ------------  ---------  ----------  ----------  ----------  -------  --------  ---------  ---------  ---------
## Spirochaetota   Treponemataceae   Treponema_D    0.161132   0.5982314   0.3558038   0.1388215        8   1.60742   10.80004   5.024315   3.173004
```



## Count eukaryotic contigs.
Count genera (after the second semicolon from the end):


``` r
blastn.coverm.megahit%>%
  filter(grepl("Eukaryota",classification))%>%
  mutate(classification=sub(".*;[^;]*;([^;]*;[^;]*)$", "\\1", classification))%>%
  separate_wider_delim(delim = ";",cols=classification,
                       names=c("Genus","Species"))%>%
  count(Genus,name = "n_contigs")%>%
  arrange(-n_contigs)%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## Genus                n_contigs
## ------------------  ----------
## Tritrichomonas              13
## Balantidium                  4
## Giardia                      4
## Pichia                       3
## Hanseniaspora                2
## Pycnora                      2
## Ipomoea                      1
## Saccharomycetales            1
```

``` r
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
##  [1] microViz_0.12.1 phyloseq_1.50.0 lubridate_1.9.4 forcats_1.0.0  
##  [5] stringr_1.5.1   dplyr_1.1.4     purrr_1.0.4     readr_2.1.5    
##  [9] tidyr_1.3.1     tibble_3.2.1    ggplot2_4.0.0   tidyverse_2.0.0
## 
## loaded via a namespace (and not attached):
##  [1] ade4_1.7-23             tidyselect_1.2.1        farver_2.1.2           
##  [4] bitops_1.0-9            Biostrings_2.74.1       S7_0.2.0               
##  [7] RCurl_1.98-1.17         fastmap_1.2.0           digest_0.6.37          
## [10] timechange_0.3.0        lifecycle_1.0.5         cluster_2.1.8.1        
## [13] survival_3.8-3          magrittr_2.0.3          compiler_4.4.3         
## [16] rlang_1.1.5             sass_0.4.10             tools_4.4.3            
## [19] utf8_1.2.4              igraph_2.1.4            yaml_2.3.12            
## [22] data.table_1.17.8       knitr_1.51              labeling_0.4.3         
## [25] here_1.0.2              xml2_1.3.8              plyr_1.8.9             
## [28] RColorBrewer_1.1-3      withr_3.0.2             BiocGenerics_0.52.0    
## [31] grid_4.4.3              stats4_4.4.3            multtest_2.62.0        
## [34] biomformat_1.34.0       colorspace_2.1-2        Rhdf5lib_1.28.0        
## [37] scales_1.4.0            iterators_1.0.14        MASS_7.3-65            
## [40] cli_3.6.4               rmarkdown_2.30          vegan_2.6-4            
## [43] crayon_1.5.3            generics_0.1.4          otel_0.2.0             
## [46] rstudioapi_0.18.0       httr_1.4.7              reshape2_1.4.4         
## [49] tzdb_0.5.0              ape_5.8-1               cachem_1.1.0           
## [52] rhdf5_2.50.2            zlibbioc_1.52.0         splines_4.4.3          
## [55] parallel_4.4.3          XVector_0.46.0          vctrs_0.6.5            
## [58] Matrix_1.7-4            jsonlite_2.0.0          bookdown_0.46          
## [61] IRanges_2.40.1          hms_1.1.4               S4Vectors_0.44.0       
## [64] foreach_1.5.2           jquerylib_0.1.4         glue_1.8.0             
## [67] codetools_0.2-20        ggtext_0.1.2            stringi_1.8.4          
## [70] gtable_0.3.6            GenomeInfoDb_1.42.3     UCSC.utils_1.2.0       
## [73] pillar_1.11.1           htmltools_0.5.8.1       rhdf5filters_1.18.1    
## [76] GenomeInfoDbData_1.2.13 R6_2.6.1                rprojroot_2.1.1        
## [79] evaluate_1.0.5          lattice_0.22-6          Biobase_2.66.0         
## [82] gridtext_0.1.5          bslib_0.10.0            Rcpp_1.0.14            
## [85] nlme_3.1-167            permute_0.9-8           mgcv_1.9-1             
## [88] xfun_0.56               pkgconfig_2.0.3
```

``` r
rm(list =setdiff(ls(all.names = TRUE), c("markdown.dir")))
gc()
```

```
##            used  (Mb) gc trigger  (Mb) max used  (Mb)
## Ncells  5671466 302.9    8824807 471.3  8824807 471.3
## Vcells 10085565  77.0   48726814 371.8 54368499 414.8
```

