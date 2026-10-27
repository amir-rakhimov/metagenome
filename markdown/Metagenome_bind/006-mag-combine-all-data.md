---
output: 
  bookdown::html_document2:
     toc: true
---


# Metagenome assembly analysis {-}



# Combine all the data from MAG assembly





## Introduction
Here, we combine BLASTN and GTDB-tk results (on all contigs from MEGAHIT) 
with seqkit stats (on contigs and dereplicated bins).




## Load necessary libraries


``` r
# install.packages(c("tidyverse", "data.table"))
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
library(data.table)
```

```
## 
## Attaching package: 'data.table'
## 
## The following objects are masked from 'package:lubridate':
## 
##     hour, isoweek, mday, minute, month, quarter, second, wday, week,
##     yday, year
## 
## The following objects are masked from 'package:dplyr':
## 
##     between, first, last
## 
## The following object is masked from 'package:purrr':
## 
##     transpose
```

``` r
# load(file="./output/rdafiles/mag-combine-all-data.Rdata")
```



## Load data.


``` r
source(here::here("config/R/config.R"))# config file with global variables
source(here::here("config/R/themes.R"))# config file with themes
dir.create(file.path(mag.rdafiles),recursive = TRUE)
```

```
## Warning in dir.create(file.path(mag.rdafiles), recursive = TRUE):
## '\\wsl.localhost\Ubuntu-22.04\home\rakhimov\projects\metagenome\analysis\2-mag\20250303_17_54_26\rdafiles'
## already exists
```

``` r
dir.create(file.path(mag.tables),recursive = TRUE)
```

```
## Warning in dir.create(file.path(mag.tables), recursive = TRUE):
## '\\wsl.localhost\Ubuntu-22.04\home\rakhimov\projects\metagenome\analysis\2-mag\20250303_17_54_26\tables'
## already exists
```

``` r
dir.create(file.path(mag.figures),recursive = TRUE)
```

```
## Warning in dir.create(file.path(mag.figures), recursive = TRUE):
## '\\wsl.localhost\Ubuntu-22.04\home\rakhimov\projects\metagenome\analysis\2-mag\20250303_17_54_26\figures'
## already exists
```



## Read the seqkit stats.
Seqkit on MEGAHIT [VERY LARGE]


``` r
seqkit.megahit.df<-data.table::fread(seqkit.megahit.fname,
                                     header = TRUE,
                              sep="\t")%>%
  as_tibble()%>%
  rename("sample_id" = "sample")
```

Seqkit on MetaBAT2


``` r
seqkit.metabat2.df<-read.table(seqkit.metabat2.fname,header = TRUE,
                              sep="\t",comment.char = "")%>%
  as_tibble()%>%
  rename("sample_id" = "sample")
```

Seqkit on dRep


``` r
seqkit.drep.df<-read.table(seqkit.drep.fname,header = TRUE,
                              sep="\t",comment.char = "")%>%
  as_tibble()%>%
  rename("sample_id" = "sample")
head(seqkit.megahit.df)%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## sample_id   contig_id      contig_length      GC      AT
## ----------  ------------  --------------  ------  ------
## 2D10        k141_0                   358    6.15   93.85
## 2D10        k141_164181              658   39.97   60.03
## 2D10        k141_97019               321   34.89   65.11
## 2D10        k141_37315               446   48.43   51.57
## 2D10        k141_111945              456   59.65   40.35
## 2D10        k141_1                   351   68.38   31.62
```

``` r
head(seqkit.metabat2.df)%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## sample_id    bin_id  contig_id      contig_length      GC      AT
## ----------  -------  ------------  --------------  ------  ------
## 2D10            100  k141_14959              8429   61.92   38.08
## 2D10            100  k141_104486             5478   56.97   43.03
## 2D10            100  k141_37383              3741   59.48   40.52
## 2D10            100  k141_223951             3933   56.88   43.12
## 2D10            100  k141_231423            41915   62.00   38.00
## 2D10            100  k141_37478              5856   61.63   38.37
```

``` r
head(seqkit.drep.df)%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## sample_id    bin_id  contig_id      contig_length      GC      AT
## ----------  -------  ------------  --------------  ------  ------
## 2D10            106  k141_142200            16466   57.49   42.51
## 2D10            106  k141_998               52537   61.47   38.53
## 2D10            106  k141_82854             25032   62.32   37.68
## 2D10            106  k141_120674            28507   62.41   37.59
## 2D10            106  k141_121046            10463   60.49   39.51
## 2D10            106  k141_90798             16624   62.53   37.47
```

3118148 contigs from MEGAHIT


``` r
nrow(seqkit.megahit.df)
```

```
## [1] 3118148
```

156182 contigs from MetaBAT2


``` r
nrow(seqkit.metabat2.df)
```

```
## [1] 156182
```

11456 contigs from dRep


``` r
nrow(seqkit.drep.df)
```

```
## [1] 11456
```

``` r
summary(seqkit.megahit.df)
```

```
##   sample_id          contig_id         contig_length          GC       
##  Length:3118148     Length:3118148     Min.   :   200   Min.   : 0.87  
##  Class :character   Class :character   1st Qu.:   383   1st Qu.:41.45  
##  Mode  :character   Mode  :character   Median :   528   Median :48.75  
##                                        Mean   :  1370   Mean   :48.15  
##                                        3rd Qu.:   868   3rd Qu.:55.59  
##                                        Max.   :821016   Max.   :99.18  
##        AT       
##  Min.   : 0.82  
##  1st Qu.:44.41  
##  Median :51.25  
##  Mean   :51.85  
##  3rd Qu.:58.55  
##  Max.   :99.13
```

``` r
summary(seqkit.metabat2.df)
```

```
##   sample_id             bin_id        contig_id         contig_length   
##  Length:156182      Min.   :  1.00   Length:156182      Min.   :  2500  
##  Class :character   1st Qu.: 35.00   Class :character   1st Qu.:  3479  
##  Mode  :character   Median : 60.00   Mode  :character   Median :  5504  
##                     Mean   : 62.16                      Mean   : 13314  
##                     3rd Qu.: 86.00                      3rd Qu.: 12143  
##                     Max.   :145.00                      Max.   :821016  
##        GC              AT       
##  Min.   :14.45   Min.   :23.68  
##  1st Qu.:41.53   1st Qu.:42.29  
##  Median :49.40   Median :50.60  
##  Mean   :48.90   Mean   :51.10  
##  3rd Qu.:57.71   3rd Qu.:58.47  
##  Max.   :76.32   Max.   :85.55
```

``` r
summary(seqkit.drep.df)
```

```
##   sample_id             bin_id        contig_id         contig_length   
##  Length:11456       Min.   :  1.00   Length:11456       Min.   :  2500  
##  Class :character   1st Qu.: 34.00   Class :character   1st Qu.:  6361  
##  Mode  :character   Median : 63.00   Mode  :character   Median : 13319  
##                     Mean   : 62.65                      Mean   : 26202  
##                     3rd Qu.: 86.00                      3rd Qu.: 29251  
##                     Max.   :142.00                      Max.   :674940  
##        GC              AT       
##  Min.   :20.84   Min.   :26.32  
##  1st Qu.:40.69   1st Qu.:45.76  
##  Median :46.94   Median :53.06  
##  Mean   :47.47   Mean   :52.53  
##  3rd Qu.:54.24   3rd Qu.:59.31  
##  Max.   :73.68   Max.   :79.16
```



## Add drep clusters information.
We want to know which bins were assigned to which clusters, and which
bins weren't because they had low quality.
We need it to get a more precise estimate of mapped reads. 
Read the data with all high-quality MAGs and their clusters


``` r
drep.clusters<-read.csv(drep.clusters.fname)%>%
  as_tibble()%>%
  mutate(genome=gsub(paste0(metabat2.date_time,"_"),"",genome),
         genome=gsub(".fa","",genome))%>%
  separate_wider_delim(genome,delim = "_bin.",names = c("sample_id","bin_id"))%>%
  mutate(bin_id=as.integer(bin_id))
head(drep.clusters)%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## sample_id    bin_id  secondary_cluster    threshold  cluster_method   comparison_algorithm    primary_cluster
## ----------  -------  ------------------  ----------  ---------------  ---------------------  ----------------
## H3               42  1_1                       0.05  average          fastANI                               1
## Y66b             15  1_1                       0.05  average          fastANI                               1
## G18              52  2_0                       0.05  average          fastANI                               2
## H21             110  3_0                       0.05  average          fastANI                               3
## 2D14             98  4_1                       0.05  average          fastANI                               4
## G14              96  4_1                       0.05  average          fastANI                               4
```

Read the data with representative genomes and their clusters


``` r
drep.winning_genomes<-read.csv(drep.winning_genomes.fname)%>%
  as_tibble()%>%
  mutate(genome=gsub(paste0(metabat2.date_time,"_"),"",genome),
         genome=gsub(".fa","",genome))%>%
  separate_wider_delim(genome,delim = "_bin.",names = c("sample_id","bin_id"))%>%
  mutate(bin_id=as.integer(bin_id))
head(drep.winning_genomes)%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## sample_id    bin_id  cluster        score
## ----------  -------  --------  ----------
## Y66b             15  1_1         94.02366
## G18              52  2_0         87.26024
## H21             110  3_0         97.50959
## 2D14             98  4_1        100.19536
## Y51b             13  5_1        101.36102
## H3              108  6_0         90.23830
```

Combine drep MAGs with other checkm2 MAGs


``` r
drep.clusters.joined<-drep.winning_genomes%>%
  select(sample_id,bin_id,cluster)%>%
  rename("secondary_cluster" = "cluster",
         "drep_sample"= "sample_id",
         "drep_bin_id"= "bin_id")%>%
  full_join(drep.clusters[,c("sample_id","bin_id","secondary_cluster")],
            join_by("secondary_cluster"))%>%
  # Add a column that shows if a bin is representative or not.
  # Then, remove drep_sample and drep_bin_id because it's the same information
  mutate(is_representative_bin=ifelse(sample_id==drep_sample&bin_id==drep_bin_id,
                            TRUE,FALSE))%>%
  select(-drep_sample,-drep_bin_id)
head(drep.clusters.joined)%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## secondary_cluster   sample_id    bin_id  is_representative_bin 
## ------------------  ----------  -------  ----------------------
## 1_1                 H3               42  FALSE                 
## 1_1                 Y66b             15  TRUE                  
## 2_0                 G18              52  TRUE                  
## 3_0                 H21             110  TRUE                  
## 4_1                 2D14             98  TRUE                  
## 4_1                 G14              96  FALSE
```

``` r
rm(drep.clusters)
rm(drep.clusters.fname)
rm(drep.winning_genomes)
rm(drep.winning_genomes.fname)
```



### Join cluster information with seqkit stats.


``` r
seqkit.drep.df<-seqkit.drep.df%>%
  left_join(drep.clusters.joined,by = join_by(sample_id, bin_id))

seqkit.metabat2.df <- seqkit.metabat2.df%>%
  left_join(drep.clusters.joined,by = join_by(sample_id, bin_id))%>%
  mutate(secondary_cluster = ifelse(is.na(secondary_cluster), "low_quality_mag",
                                    secondary_cluster),
         is_representative_bin = ifelse(is.na(is_representative_bin), "low_quality_mag",
                                        is_representative_bin)
         )
  
head(seqkit.drep.df)%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## sample_id    bin_id  contig_id      contig_length      GC      AT  secondary_cluster   is_representative_bin 
## ----------  -------  ------------  --------------  ------  ------  ------------------  ----------------------
## 2D10            106  k141_142200            16466   57.49   42.51  47_1                TRUE                  
## 2D10            106  k141_998               52537   61.47   38.53  47_1                TRUE                  
## 2D10            106  k141_82854             25032   62.32   37.68  47_1                TRUE                  
## 2D10            106  k141_120674            28507   62.41   37.59  47_1                TRUE                  
## 2D10            106  k141_121046            10463   60.49   39.51  47_1                TRUE                  
## 2D10            106  k141_90798             16624   62.53   37.47  47_1                TRUE
```

``` r
head(seqkit.metabat2.df)
```

```
## # A tibble: 6 × 8
##   sample_id bin_id contig_id   contig_length    GC    AT secondary_cluster
##   <chr>      <int> <chr>               <int> <dbl> <dbl> <chr>            
## 1 2D10         100 k141_14959           8429  61.9  38.1 low_quality_mag  
## 2 2D10         100 k141_104486          5478  57.0  43.0 low_quality_mag  
## 3 2D10         100 k141_37383           3741  59.5  40.5 low_quality_mag  
## 4 2D10         100 k141_223951          3933  56.9  43.1 low_quality_mag  
## 5 2D10         100 k141_231423         41915  62    38   low_quality_mag  
## 6 2D10         100 k141_37478           5856  61.6  38.4 low_quality_mag  
## # ℹ 1 more variable: is_representative_bin <chr>
```

``` r
rm(drep.clusters.joined)
```



## Combine all seqkit statistics into one table.
First, we create a column with each contig + its bin + sample_id.


``` r
seqkit.metabat2.df <- seqkit.metabat2.df%>%
  mutate(sample_contig=paste(sample_id,contig_id,sep = "_"),
         sample_bin_contig=paste(sample_id,bin_id,contig_id,sep = "_"))
seqkit.drep.df <- seqkit.drep.df%>%
  mutate(sample_contig=paste(sample_id,contig_id,sep = "_"),
         sample_bin_contig=paste(sample_id,bin_id,contig_id,sep = "_"))
seqkit.megahit.df <- seqkit.megahit.df%>%
  mutate(sample_contig=paste(sample_id,contig_id,sep = "_"))
head(seqkit.drep.df)%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## sample_id    bin_id  contig_id      contig_length      GC      AT  secondary_cluster   is_representative_bin   sample_contig      sample_bin_contig    
## ----------  -------  ------------  --------------  ------  ------  ------------------  ----------------------  -----------------  ---------------------
## 2D10            106  k141_142200            16466   57.49   42.51  47_1                TRUE                    2D10_k141_142200   2D10_106_k141_142200 
## 2D10            106  k141_998               52537   61.47   38.53  47_1                TRUE                    2D10_k141_998      2D10_106_k141_998    
## 2D10            106  k141_82854             25032   62.32   37.68  47_1                TRUE                    2D10_k141_82854    2D10_106_k141_82854  
## 2D10            106  k141_120674            28507   62.41   37.59  47_1                TRUE                    2D10_k141_120674   2D10_106_k141_120674 
## 2D10            106  k141_121046            10463   60.49   39.51  47_1                TRUE                    2D10_k141_121046   2D10_106_k141_121046 
## 2D10            106  k141_90798             16624   62.53   37.47  47_1                TRUE                    2D10_k141_90798    2D10_106_k141_90798
```

``` r
head(seqkit.metabat2.df)%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## sample_id    bin_id  contig_id      contig_length      GC      AT  secondary_cluster   is_representative_bin   sample_contig      sample_bin_contig    
## ----------  -------  ------------  --------------  ------  ------  ------------------  ----------------------  -----------------  ---------------------
## 2D10            100  k141_14959              8429   61.92   38.08  low_quality_mag     low_quality_mag         2D10_k141_14959    2D10_100_k141_14959  
## 2D10            100  k141_104486             5478   56.97   43.03  low_quality_mag     low_quality_mag         2D10_k141_104486   2D10_100_k141_104486 
## 2D10            100  k141_37383              3741   59.48   40.52  low_quality_mag     low_quality_mag         2D10_k141_37383    2D10_100_k141_37383  
## 2D10            100  k141_223951             3933   56.88   43.12  low_quality_mag     low_quality_mag         2D10_k141_223951   2D10_100_k141_223951 
## 2D10            100  k141_231423            41915   62.00   38.00  low_quality_mag     low_quality_mag         2D10_k141_231423   2D10_100_k141_231423 
## 2D10            100  k141_37478              5856   61.63   38.37  low_quality_mag     low_quality_mag         2D10_k141_37478    2D10_100_k141_37478
```

``` r
head(seqkit.megahit.df)%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## sample_id   contig_id      contig_length      GC      AT  sample_contig    
## ----------  ------------  --------------  ------  ------  -----------------
## 2D10        k141_0                   358    6.15   93.85  2D10_k141_0      
## 2D10        k141_164181              658   39.97   60.03  2D10_k141_164181 
## 2D10        k141_97019               321   34.89   65.11  2D10_k141_97019  
## 2D10        k141_37315               446   48.43   51.57  2D10_k141_37315  
## 2D10        k141_111945              456   59.65   40.35  2D10_k141_111945 
## 2D10        k141_1                   351   68.38   31.62  2D10_k141_1
```

``` r
nrow(seqkit.drep.df)
```

```
## [1] 11456
```

``` r
nrow(seqkit.metabat2.df)
```

```
## [1] 156182
```

``` r
nrow(seqkit.megahit.df)
```

```
## [1] 3118148
```

Get the common IDs between datasets to know whether a contig
was assembled (metabat2 or metabat2_drep) or not (megahit).
If it's metabat2, it is not representative. If it's metabat2_drep,
it is representative.


``` r
megahit.metabat2.common.ids<-intersect(seqkit.megahit.df$sample_contig,
                                      seqkit.metabat2.df$sample_contig)
metabat2.drep.common.ids<-intersect(seqkit.metabat2.df$sample_bin_contig,
                                   seqkit.drep.df$sample_bin_contig)
megahit.drep.common.ids<-intersect(seqkit.megahit.df$sample_contig,
                                   seqkit.drep.df$sample_contig)
length(megahit.metabat2.common.ids)
```

```
## [1] 156182
```

``` r
length(metabat2.drep.common.ids)
```

```
## [1] 11456
```

``` r
length(megahit.drep.common.ids)
```

```
## [1] 11456
```

All contigs from drep are found in MetaBAT2 (no differences):


``` r
length(setdiff(seqkit.drep.df$sample_bin_contig,
        seqkit.metabat2.df$sample_bin_contig))
```

```
## [1] 0
```

All contigs from MetaBAT2 are found in MEGAHIT (no differences):


``` r
length(setdiff(seqkit.metabat2.df$sample_contig,
        seqkit.megahit.df$sample_contig))
```

```
## [1] 0
```

Finally combine the data:


``` r
seqkit.all.df <- seqkit.megahit.df %>%
  full_join(subset(seqkit.metabat2.df,select = -c(AT,GC,sample_contig,
                                                      contig_length)),
            by = join_by(sample_id, contig_id))%>%
  relocate(sample_id,contig_id,bin_id)%>%
  # If bin_id, secondary_cluster, and is_representative_bin are NA, the contig was not assembled
  mutate(bin_id = ifelse(is.na(bin_id), "unassembled", bin_id),
         secondary_cluster = ifelse(is.na(secondary_cluster), "unassembled", secondary_cluster),
         is_representative_bin = ifelse(is.na(is_representative_bin), "unassembled", is_representative_bin))%>%
  select(-sample_contig,-sample_bin_contig)
head(seqkit.all.df)%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## sample_id   contig_id     bin_id         contig_length      GC      AT  secondary_cluster   is_representative_bin 
## ----------  ------------  ------------  --------------  ------  ------  ------------------  ----------------------
## 2D10        k141_0        unassembled              358    6.15   93.85  unassembled         unassembled           
## 2D10        k141_164181   unassembled              658   39.97   60.03  unassembled         unassembled           
## 2D10        k141_97019    unassembled              321   34.89   65.11  unassembled         unassembled           
## 2D10        k141_37315    unassembled              446   48.43   51.57  unassembled         unassembled           
## 2D10        k141_111945   unassembled              456   59.65   40.35  unassembled         unassembled           
## 2D10        k141_1        unassembled              351   68.38   31.62  unassembled         unassembled
```

``` r
nrow(seqkit.all.df)
```

```
## [1] 3118148
```

``` r
str(seqkit.all.df)
```

```
## tibble [3,118,148 × 8] (S3: tbl_df/tbl/data.frame)
##  $ sample_id            : chr [1:3118148] "2D10" "2D10" "2D10" "2D10" ...
##  $ contig_id            : chr [1:3118148] "k141_0" "k141_164181" "k141_97019" "k141_37315" ...
##  $ bin_id               : chr [1:3118148] "unassembled" "unassembled" "unassembled" "unassembled" ...
##  $ contig_length        : int [1:3118148] 358 658 321 446 456 351 1154 1561 232 385 ...
##  $ GC                   : num [1:3118148] 6.15 39.97 34.89 48.43 59.65 ...
##  $ AT                   : num [1:3118148] 93.8 60 65.1 51.6 40.4 ...
##  $ secondary_cluster    : chr [1:3118148] "unassembled" "unassembled" "unassembled" "unassembled" ...
##  $ is_representative_bin: chr [1:3118148] "unassembled" "unassembled" "unassembled" "unassembled" ...
##  - attr(*, ".internal.selfref")=<externalptr>
```

``` r
summary(seqkit.all.df)
```

```
##   sample_id          contig_id            bin_id          contig_length   
##  Length:3118148     Length:3118148     Length:3118148     Min.   :   200  
##  Class :character   Class :character   Class :character   1st Qu.:   383  
##  Mode  :character   Mode  :character   Mode  :character   Median :   528  
##                                                           Mean   :  1370  
##                                                           3rd Qu.:   868  
##                                                           Max.   :821016  
##        GC              AT        secondary_cluster  is_representative_bin
##  Min.   : 0.87   Min.   : 0.82   Length:3118148     Length:3118148       
##  1st Qu.:41.45   1st Qu.:44.41   Class :character   Class :character     
##  Median :48.75   Median :51.25   Mode  :character   Mode  :character     
##  Mean   :48.15   Mean   :51.85                                           
##  3rd Qu.:55.59   3rd Qu.:58.55                                           
##  Max.   :99.18   Max.   :99.13
```

``` r
rm(megahit.drep.common.ids)
rm(megahit.metabat2.common.ids)
rm(metabat2.drep.common.ids)
```

Barplot to count how many contigs were retained: most contigs weren't assembled


``` r
seqkit.all.df%>%
  mutate(assembly_status = ifelse(! secondary_cluster %in% c("low_quality_mag",
                                                               "unassembled"),
                                    "high_quality_mag", secondary_cluster))%>%
  count(sample_id,assembly_status )%>%
  # pivot_wider(names_from = "assembly_status",
  #             values_from="n")%>%
  ggplot(aes(x=sample_id,fill=assembly_status,y=n))+
  geom_bar(stat="identity")
```

<img src="006-mag-combine-all-data_files/figure-html/unnamed-chunk-26-1.png" alt="" width="672" />

``` r
rm(seqkit.drep.df)
rm(seqkit.megahit.df)
rm(seqkit.metabat2.df)
gc()
```

```
##            used  (Mb) gc trigger  (Mb) max used  (Mb)
## Ncells  2065489 110.4    7786230 415.9  6929383 370.1
## Vcells 29150506 222.5   99543783 759.5 98169460 749.0
```



## Add taxonomy to the seqkit stats data.




### Add GTDB-tk taxonomy (for dereplicated bins).


``` r
gtdbtk.taxonomy<-read.table(gtdbtk.taxonomy.fname,header = T,
                            sep="\t",comment.char = "")%>%
  as_tibble()%>%
  select(user_genome,classification) # remove unnecessary columns
head(gtdbtk.taxonomy)%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## user_genome                      classification                                                                                       
## -------------------------------  -----------------------------------------------------------------------------------------------------
## 20250612_13_37_47_2D10_bin.10    d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Muribaculaceae;g__Duncaniella;s__     
## 20250612_13_37_47_2D10_bin.106   d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Muribaculaceae;g__CAG-873;s__         
## 20250612_13_37_47_2D10_bin.11    d__Bacteria;p__Bacillota;c__Clostridia;o__Oscillospirales;f__Acutalibacteraceae;g__Acutalibacter;s__ 
## 20250612_13_37_47_2D10_bin.26    d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__      
## 20250612_13_37_47_2D10_bin.39    d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Prevotella;s__      
## 20250612_13_37_47_2D10_bin.42    d__Bacteria;p__Bacillota;c__Clostridia;o__Christensenellales;f__Aristaeellaceae;g__;s__
```

``` r
# 128 bins
nrow(gtdbtk.taxonomy)
```

```
## [1] 128
```

Tidy up GTDB-tk taxonomy: substitute s__ and g__ with unclassified
Separate sample_bin column into separate columns


``` r
gtdbtk.taxonomy<-gtdbtk.taxonomy %>%
  mutate(classification=gsub(";g__;s__",";unclassified;unclassified",
                             classification))%>%
  mutate(classification=gsub(";s__$",";unclassified",classification))%>%
  mutate(classification=gsub("d__|p__|c__|o__|f__|g__|s__","",classification)) %>%
  mutate(user_genome=gsub(paste0(metabat2.date_time,"_"),"",user_genome))%>%
  separate_wider_delim(user_genome,delim="_bin.",names=c("sample_id","bin_id"))%>%
  mutate(bin_id=as.integer(bin_id))
```



### Add BLAST taxonomy (for all contigs from MEGAHIT).


``` r
blastn.taxonomy<-read.table(blastn.taxonomy.fname,header = F,
                            sep="\t",comment.char = "")%>%
  as_tibble()
head(blastn.taxonomy)%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## V1     V2            V3                          V4     V5     V6  V7                                                                                                                          
## -----  ------------  ---------------------  -------  -----  -----  ----------------------------------------------------------------------------------------------------------------------------
## 2D10   k141_105068   EF400498.1.1478         95.611   1162   1857  Bacteria;Bacillota;Clostridia;Peptostreptococcales-Tissierellales;Anaerovoracaceae;Family XIII UCG-001;uncultured bacterium 
## 2D10   k141_10535    EU463109.1.1398         98.784   1398   2488  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Prevotellaceae;Prevotellaceae UCG-001;uncultured bacterium                  
## 2D10   k141_106439   KY635417.1.1502         91.935    992   1382  Bacteria;Thermodesulfobacteriota;Desulfovibrionia;Desulfovibrionales;Desulfovibrionaceae;Bilophila;uncultured bacterium     
## 2D10   k141_107540   CP015406.84360.85899    95.040    504    791  Bacteria;Bacillota;Clostridia;Oscillospirales;Oscillospiraceae;Flavonifractor;Flavonifractor plautii                        
## 2D10   k141_108239   JN162690.1.1347         95.085    529    832  Bacteria;Actinomycetota;Coriobacteriia;Coriobacteriales;Atopobiaceae;Incertae Sedis;uncultured bacterium                    
## 2D10   k141_108435   EU622679.1.1499         98.679    984   1746  Bacteria;Bacillota;Clostridia;Lachnospirales;Lachnospiraceae;[Eubacterium] ventriosum group;uncultured bacterium
```

``` r
# 1899 contigs
nrow(blastn.taxonomy)
```

```
## [1] 1899
```

``` r
colnames(blastn.taxonomy)<-c("sample_id","contig_id","silva_ref_id","perc_identity",
                             "alignment_length","bitscore","classification")
```



### Join stats with taxonomy.


``` r
###NOT NEEDED?
# seqkit.with_taxonomy<-seqkit.all.df%>%
#   full_join(gtdbtk.taxonomy,
#             by = join_by(sample_id, bin_id))%>%
#   full_join(blastn.taxonomy[,c("sample_id","contig_id","classification")],
#             by = join_by(sample_id, contig_id))%>%
#   rename("gtdbtk_result"="classification.x",
#          "blastn_result"="classification.y")%>%
#   ### 6.4 Add column to show if a contig was classified, and which pipeline classified it ####
#   mutate(tax.classification.type=ifelse(!is.na(gtdbtk_result)&!is.na(blastn_result),
#                         "gtdbtk_and_blastn","unclassified"),
#          tax.classification.type=ifelse(!is.na(gtdbtk_result)&is.na(blastn_result),
#                         "gtdbtk_only",tax.classification.type),
#          tax.classification.type=ifelse(is.na(gtdbtk_result)&!is.na(blastn_result),
#                         "blastn_only",tax.classification.type))
# nrow(seqkit.with_taxonomy)
# head(seqkit.with_taxonomy)
# str(seqkit.with_taxonomy)
# summary(seqkit.with_taxonomy)
# rm(seqkit.all.df)
```



## Add bin and contig coverages from CoverM.




### CoverM data from dRep.


``` r
coverm.drep.df<-read.table(coverm.drep.fname,header = T,
                           sep = "\t")%>%
  as_tibble()%>%
  mutate(genome=gsub(paste0(metabat2.date_time,"_"),"",genome))
head(coverm.drep.df)%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## genome           relative_abundance        mean   trimmed_mean   covered_fraction   covered_bases    variance    length   read_count   reads_per_base        rpkm         tpm
## --------------  -------------------  ----------  -------------  -----------------  --------------  ----------  --------  -----------  ---------------  ----------  ----------
## unmapped_2D10            87.1793300          NA             NA                 NA              NA          NA        NA           NA               NA          NA          NA
## 2D10_bin.106              1.8713269    82.58569       78.43517          0.9999754         2110166   1539.1360   2110218      1225708        0.5808443    60.88153   147317.61
## 2D10_bin.10               0.9019371    39.80443       36.12676          0.9999700         2401556    938.5855   2401628       685978        0.2856304    29.93852    72443.48
## 2D10_bin.11               0.4851722    21.41170       20.91922          1.0000000         1193490    150.0317   1193490       172989        0.1449438    15.19237    36761.61
## 2D10_bin.26               1.3692596    60.42838       55.94852          0.9999195         3378337   1459.7344   3378609      1469758        0.4350187    45.59673   110332.31
## 2D10_bin.39               3.3015501   145.70453      142.67023          0.9998766         2915798   2382.5125   2916158      2906398        0.9966531   104.46478   252777.78
```

``` r
nrow(coverm.drep.df)
```

```
## [1] 139
```



### CoverM data from MetaBAT2.


``` r
coverm.metabat2.df<-read.table(coverm.metabat2.fname,header = T,
                               sep = "\t")%>%
  as_tibble()%>%
  mutate(genome=gsub(paste0(metabat2.date_time,"_"),"",genome))
head(coverm.metabat2.df)%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## genome           relative_abundance        mean   trimmed_mean   covered_fraction   covered_bases    variance    length   read_count   reads_per_base        rpkm        tpm
## --------------  -------------------  ----------  -------------  -----------------  --------------  ----------  --------  -----------  ---------------  ----------  ---------
## unmapped_2D10            20.8537700          NA             NA                 NA              NA          NA        NA           NA               NA          NA         NA
## 2D10_bin.100              0.4797199   31.469807      30.401953          0.9998561         2251352   388.04904   2251676       481084        0.2136560   3.6276174   6098.809
## 2D10_bin.101              0.4770218   31.292807      31.020714          1.0000000          811662    74.06812    811662       169735        0.2091203   3.5506072   5969.336
## 2D10_bin.102              0.1502874    9.858910       9.658344          0.9998679         1498842    19.62712   1499040        99069        0.0660883   1.1220986   1886.491
## 2D10_bin.103              0.0999692    6.558015       6.316923          0.9998813          606527    13.67829    606599        26712        0.0440357   0.7476721   1256.998
## 2D10_bin.104              0.1075584    7.055874       6.653600          0.9999239          525867    32.74173    525907        25613        0.0487025   0.8269095   1390.213
```

``` r
nrow(coverm.metabat2.df)
```

```
## [1] 1352
```



### CoverM data from MEGAHIT (use data.table package).


``` r
coverm.megahit.df<-data.table::fread(coverm.megahit.fname,header = T,
                                  sep = "\t")%>%
  as_tibble()%>%
  rename(sample_contig=contig)%>%
  separate_wider_delim(sample_contig,delim = "_",
                       names = c("sample_id","contig_id"),
                       too_many="merge" )
head(coverm.megahit.df)%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## sample_id   contig_id             mean   trimmed_mean   covered_fraction   covered_bases       variance   length   read_count   reads_per_base          rpkm            tpm
## ----------  ------------  ------------  -------------  -----------------  --------------  -------------  -------  -----------  ---------------  ------------  -------------
## 2D10        k141_0         2150.101000    2156.069000                  1             358   2.036667e+06      358         5671       15.8407820   220.8570000   1222.1907000
## 2D10        k141_164181       3.305118       3.320961                  1             658   1.506325e+00      658           15        0.0227964     0.3178337      1.7588462
## 2D10        k141_97019        1.736842       1.780645                  1             321   1.950464e-01      321            4        0.0124611     0.1737359      0.9614291
## 2D10        k141_37315        4.581081       4.574627                  1             446   3.044251e+00      446           13        0.0291480     0.4063900      2.2489035
## 2D10        k141_111945       3.866013       3.851449                  1             456   1.742644e+00      456           12        0.0263158     0.3669027      2.0303862
## 2D10        k141_1            6.487562       6.607735                  1             351   5.091095e+00      351           12        0.0341880     0.4766600      2.6377685
```

``` r
nrow(coverm.megahit.df)
```

```
## [1] 3118148
```

How much relative abundance is classified in drep bins:


``` r
coverm.drep.df%>%
  mutate(sample_id=ifelse(grepl("unmapped",genome),
                       gsub("unmapped_","",genome),
                       gsub("_bin.[0-9]*","",genome)))%>%
  filter(!grepl("unmapped",genome))%>%
  group_by(sample_id)%>%
  summarise(sum_relative_ab=sum(relative_abundance))
```

```
## # A tibble: 11 × 2
##    sample_id sum_relative_ab
##    <chr>               <dbl>
##  1 2D10                12.8 
##  2 2D14                 6.24
##  3 G14                  3.82
##  4 G18                 20.3 
##  5 H15                  8.84
##  6 H21                  7.36
##  7 H3                   4.26
##  8 H4                  11.8 
##  9 O15                  5.88
## 10 Y51b                 5.53
## 11 Y66b                 3.79
```

How much relative abundance is classified in MetaBAT2 bins:


``` r
coverm.metabat2.df%>%
  mutate(sample_id=ifelse(grepl("unmapped",genome),
                       gsub("unmapped_","",genome),
                       gsub("_bin.[0-9]*","",genome)))%>%
  filter(!grepl("unmapped",genome))%>%
  group_by(sample_id)%>%
  summarise(sum_relative_ab=sum(relative_abundance))
```

```
## # A tibble: 11 × 2
##    sample_id sum_relative_ab
##    <chr>               <dbl>
##  1 2D10                 79.1
##  2 2D14                 80.0
##  3 G14                  80.2
##  4 G18                  78.7
##  5 H15                  77.0
##  6 H21                  80.3
##  7 H3                   77.7
##  8 H4                   79.4
##  9 O15                  79.4
## 10 Y51b                 75.1
## 11 Y66b                 74.1
```

How much relative abundance is classified in MetaBAT2 bins that belonged to
drep clusters:


``` r
# seqkit.with_taxonomy%>%
seqkit.all.df%>%
  filter(!secondary_cluster %in% c("unassembled", "low_quality_mag"))%>%
  mutate(genome=paste(sample_id,bin_id,sep="_bin."))%>%
  distinct(sample_id,bin_id,genome,secondary_cluster)%>%
  left_join(coverm.metabat2.df,by=join_by("genome"))%>%
  group_by(sample_id)%>%
  summarise(sum_relative_ab=sum(relative_abundance))
```

```
## # A tibble: 11 × 2
##    sample_id sum_relative_ab
##    <chr>               <dbl>
##  1 2D10                21.7 
##  2 2D14                15.0 
##  3 G14                  9.75
##  4 G18                 18.5 
##  5 H15                 23.4 
##  6 H21                 15.8 
##  7 H3                  18.7 
##  8 H4                  24.6 
##  9 O15                  8.53
## 10 Y51b                12.3 
## 11 Y66b                17.0
```



## Join relative abundances with taxonomy.




### GTDB-tk + dRep


``` r
gtdbtk.coverm.drep <- seqkit.all.df%>% # seqkit.with_taxonomy%>%
  filter(!secondary_cluster %in% c("unassembled", "low_quality_mag"))%>%
  select(sample_id,bin_id,secondary_cluster)%>%
  distinct()%>%
  mutate(bin_id = as.integer(bin_id))%>% # conver bin_id into integer because we will join it with the gtdbtk data
  right_join(gtdbtk.taxonomy)%>%
  mutate(genome=paste(sample_id,bin_id,sep="_bin."))%>%
  left_join(coverm.drep.df[,c("genome","relative_abundance")],
            by=join_by("genome"))%>% # here, we add the relative abundances
  select(-genome)
```

```
## Joining with `by = join_by(sample_id, bin_id)`
```

``` r
head(gtdbtk.coverm.drep)%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## sample_id    bin_id  secondary_cluster   classification                                                                                   relative_abundance
## ----------  -------  ------------------  ----------------------------------------------------------------------------------------------  -------------------
## 2D10             42  69_0                Bacteria;Bacillota;Clostridia;Christensenellales;Aristaeellaceae;unclassified;unclassified                0.3339523
## 2D10             26  60_0                Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Bacteroidaceae;Prevotella;unclassified                    1.3692596
## 2D10             62  28_1                Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;UBA932;Cryptobacteroides;unclassified                     0.5458188
## 2D10             39  54_1                Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Bacteroidaceae;Prevotella;unclassified                    3.3015501
## 2D10             71  92_0                Bacteria;Bacillota;Clostridia;Oscillospirales;Acutalibacteraceae;Ruminococcoides;unclassified             0.2615794
## 2D10             76  105_0               Bacteria;Spirochaetota;Spirochaetia;Treponematales;Treponemataceae;Treponema_D;unclassified               0.4004118
```



### GTDB-tk + MetaBAT2


``` r
gtdbtk.coverm.metabat2 <- seqkit.all.df%>% # seqkit.with_taxonomy%>%
  filter(!secondary_cluster %in% c("unassembled", "low_quality_mag"))%>%
  distinct(sample_id,bin_id,secondary_cluster)%>%
  mutate(bin_id = as.integer(bin_id))%>%
  # full_join(gtdbtk.taxonomy,by = join_by(sample_id, bin_id)) %>%
  full_join(gtdbtk.coverm.drep[,c("secondary_cluster","classification")],
            by = join_by(secondary_cluster)) %>% # add the previous dataset because it has secondary_cluster info
  group_by(secondary_cluster) %>%
  mutate(classification = classification[!is.na(classification)][1L])%>% # remove NAs
  ungroup%>%
  mutate(genome=paste(sample_id,bin_id,sep="_bin."))%>%
  left_join(coverm.metabat2.df[,c("genome","relative_abundance")], 
            by=join_by("genome"))%>% # add relative abundances
  select(-genome)
head(gtdbtk.coverm.metabat2) %>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## sample_id    bin_id  secondary_cluster   classification                                                                                     relative_abundance
## ----------  -------  ------------------  ------------------------------------------------------------------------------------------------  -------------------
## 2D10             52  57_1                Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Marinilabiliaceae;W3P20-009;unclassified                    0.1983270
## 2D10             42  69_0                Bacteria;Bacillota;Clostridia;Christensenellales;Aristaeellaceae;unclassified;unclassified                  0.2103324
## 2D10             89  108_1               Bacteria;Fibrobacterota;Fibrobacteria;Fibrobacterales;Fibrobacteraceae;Fibrobacter;unclassified             4.8582480
## 2D10              3  88_1                Bacteria;Bacillota;Clostridia;Lachnospirales;Lachnospiraceae;Anaerostipes;Anaerostipes caccae               0.1971572
## 2D10             45  59_1                Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Bacteroidaceae;Prevotella;unclassified                      0.5177506
## 2D10             32  45_1                Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Muribaculaceae;CAG-873;unclassified                         0.7969444
```



### BLASTN + MEGAHIT


``` r
blastn.coverm.megahit <- coverm.megahit.df%>%
  select(sample_id, contig_id,tpm)%>%
  inner_join(blastn.taxonomy,by = join_by(sample_id, contig_id))
head(blastn.coverm.megahit)%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## sample_id   contig_id            tpm  silva_ref_id             perc_identity   alignment_length   bitscore  classification                                                                                                 
## ----------  ------------  ----------  ----------------------  --------------  -----------------  ---------  ---------------------------------------------------------------------------------------------------------------
## 2D10        k141_82202      3.633326  KM098109.1.1519                 94.706                510        791  Bacteria;Bacillota;Clostridia;Oscillospirales;Ruminococcaceae;Ruthenibacterium;Ruthenibacterium lactatiformans 
## 2D10        k141_119595    18.049273  EF404262.1.1485                 99.811                528        970  Bacteria;Bacillota;Clostridia;Oscillospirales;Ruminococcaceae;Acutalibacter;uncultured bacterium               
## 2D10        k141_231553    13.820761  EU777025.1.1401                 99.453                548        996  Bacteria;Bacteroidota;Bacteroidia;Bacteroidales;Prevotellaceae;Xylanibacter;uncultured bacterium               
## 2D10        k141_119623    38.320625  NFZD01007643.88.1579            94.703               1246       1930  Bacteria;Bacillota;Bacilli;Erysipelotrichales;Erysipelatoclostridiaceae;Asteroleplasma;gut metagenome          
## 2D10        k141_127099     1.975347  AM982589.1.1482                 93.514                663        985  Bacteria;Bacillota;Bacilli;Erysipelotrichales;Erysipelatoclostridiaceae;UCG-004;uncultured bacterium           
## 2D10        k141_97310      5.611254  FBXW01000008.267.1775           99.545                660       1203  Bacteria;Bacillota;Clostridia;Oscillospirales;Oscillospiraceae;Incertae Sedis;Intestinimonas sp. GM5
```

``` r
# Remove CoverM data because it's already join with taxonomy
rm(coverm.drep.df)
rm(coverm.metabat2.df)
rm(coverm.megahit.df)
rm(blastn.taxonomy)
```

Save data


``` r
if (!file.exists(seqkit.all.df.fname.rds) & !file.exists(seqkit.all.df.fname.tsv)){
  saveRDS(seqkit.all.df, file = seqkit.all.df.fname.rds)
  data.table::fwrite(seqkit.all.df,
              file = seqkit.all.df.fname.tsv,
              col.names = TRUE,
              row.names = FALSE,
              quote = FALSE,
              sep = '\t',
              na =NA)
  
}

if (!file.exists(gtdbtk.coverm.drep.fname.rds) & !file.exists(gtdbtk.coverm.drep.fname.tsv)){
  saveRDS(gtdbtk.coverm.drep,file = gtdbtk.coverm.drep.fname.rds)
  write.table(gtdbtk.coverm.drep,
              file = gtdbtk.coverm.drep.fname.tsv,
              col.names = TRUE,
              row.names = FALSE,
              quote = FALSE,
              sep = '\t')
  
}
if (!file.exists(gtdbtk.coverm.metabat2.fname.rds) & !file.exists(gtdbtk.coverm.metabat2.fname.tsv)){
  saveRDS(gtdbtk.coverm.metabat2,file = gtdbtk.coverm.metabat2.fname.rds)
  write.table(gtdbtk.coverm.metabat2,
              file = gtdbtk.coverm.metabat2.fname.tsv,
              col.names = TRUE,
              row.names = FALSE,
              quote = FALSE,
              sep = '\t')
  
}
if (!file.exists(blastn.coverm.megahit.fname.rds) & !file.exists(blastn.coverm.megahit.fname.tsv)){
  saveRDS(blastn.coverm.megahit,file = blastn.coverm.megahit.fname.rds)
  write.table(blastn.coverm.megahit,
              file = blastn.coverm.megahit.fname.tsv,
              col.names = TRUE,
              row.names = FALSE,
              quote = FALSE,
              sep = '\t')
  
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
##  [1] data.table_1.17.8 lubridate_1.9.4   forcats_1.0.0     stringr_1.5.1    
##  [5] dplyr_1.1.4       purrr_1.0.4       readr_2.1.5       tidyr_1.3.1      
##  [9] tibble_3.2.1      ggplot2_4.0.0     tidyverse_2.0.0  
## 
## loaded via a namespace (and not attached):
##  [1] sass_0.4.10        utf8_1.2.4         generics_0.1.4     xml2_1.3.8        
##  [5] stringi_1.8.4      hms_1.1.4          digest_0.6.37      magrittr_2.0.3    
##  [9] evaluate_1.0.5     grid_4.4.3         timechange_0.3.0   RColorBrewer_1.1-3
## [13] bookdown_0.46      fastmap_1.2.0      rprojroot_2.1.1    jsonlite_2.0.0    
## [17] ggtext_0.1.2       scales_1.4.0       jquerylib_0.1.4    cli_3.6.4         
## [21] rlang_1.1.5        withr_3.0.2        cachem_1.1.0       yaml_2.3.12       
## [25] otel_0.2.0         tools_4.4.3        tzdb_0.5.0         here_1.0.2        
## [29] vctrs_0.6.5        R6_2.6.1           lifecycle_1.0.5    pkgconfig_2.0.3   
## [33] pillar_1.11.1      bslib_0.10.0       gtable_0.3.6       glue_1.8.0        
## [37] Rcpp_1.0.14        xfun_0.56          tidyselect_1.2.1   rstudioapi_0.18.0 
## [41] knitr_1.51         farver_2.1.2       htmltools_0.5.8.1  labeling_0.4.3    
## [45] rmarkdown_2.30     compiler_4.4.3     S7_0.2.0           gridtext_0.1.5
```

``` r
rm(list = ls(all=TRUE))
gc()
```

```
##           used (Mb) gc trigger  (Mb)  max used  (Mb)
## Ncells 1696615 90.7   11943592 637.9  12980076 693.3
## Vcells 5024806 38.4   95626032 729.6 119532539 912.0
```

