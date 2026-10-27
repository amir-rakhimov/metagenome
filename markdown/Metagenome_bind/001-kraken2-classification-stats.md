---
output: 
  bookdown::html_document2:
     toc: true
---


# Whole metagenome sequencing (shotgun sequencing) analysis {-}



# Create a table with classification statistics





## Introduction
This script creates a table that shows percentage of unclassified reads 
per sample. The number of unclassified reads was taken from the job report.


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
source(here::here("config/R/config.R"))# config file with global variables
dir.create(community.composition.tables,recursive = TRUE)
```

```
## Warning in dir.create(community.composition.tables, recursive = TRUE):
## '\\wsl.localhost\Ubuntu-22.04\home\rakhimov\projects\metagenome\analysis\1-kraken2\20240409_17_32_40\01-community_composition\tables'
## already exists
```

``` r
dir.create(community.composition.rdafiles,recursive = TRUE)
```

```
## Warning in dir.create(community.composition.rdafiles, recursive = TRUE):
## '\\wsl.localhost\Ubuntu-22.04\home\rakhimov\projects\metagenome\analysis\1-kraken2\20240409_17_32_40\01-community_composition\rdafiles'
## already exists
```

``` r
dir.create(community.composition.figures,recursive = TRUE)
```

```
## Warning in dir.create(community.composition.figures, recursive = TRUE):
## '\\wsl.localhost\Ubuntu-22.04\home\rakhimov\projects\metagenome\analysis\1-kraken2\20240409_17_32_40\01-community_composition\figures'
## already exists
```

``` r
kraken2.results<-read.table(kraken2.unclassified.stats.fname,header = T)
kraken2.results<-kraken2.results%>%
  mutate(ClassifiedRate=Classified/Total*100)%>%
  mutate(meanClassifiedRate=round(mean(ClassifiedRate)))%>%
  mutate(SDClassifiedRate=round(sd(ClassifiedRate),3))
knitr::kable(kraken2.results, format = "simple")%>%
  print()
```

```
## 
## 
## Sample    Unclassified   Classified      Total   ClassifiedRate   meanClassifiedRate   SDClassifiedRate
## -------  -------------  -----------  ---------  ---------------  -------------------  -----------------
## 2D10           2380423     34827316   37207739         93.60234                   93              0.428
## 2D14           3564640     42936697   46501337         92.33433                   93              0.428
## G14            3486934     45293219   48780153         92.85174                   93              0.428
## G18            2524918     35471203   37996121         93.35480                   93              0.428
## H15            3556683     44364801   47921484         92.57810                   93              0.428
## H21            2787507     39540507   42328014         93.41451                   93              0.428
## H3             2896352     40294433   43190785         93.29405                   93              0.428
## H4             2687232     38383546   41070778         93.45707                   93              0.428
## O15            3052380     39265255   42317635         92.78698                   93              0.428
## Y51b           3106801     39484956   42591757         92.70563                   93              0.428
## Y66b           2830735     35846488   38677223         92.68113                   93              0.428
```

``` r
classification.stats.fname <- file.path(community.composition.tables,
                                        "kraken2-classification-stats.tsv")
if (! file.exists(classification.stats.fname)){
  write.table(kraken2.results, classification.stats.fname,
              row.names = F,quote = F,sep = "\t")
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
##  [1] lubridate_1.9.4 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4    
##  [5] purrr_1.0.4     readr_2.1.5     tidyr_1.3.1     tibble_3.2.1   
##  [9] ggplot2_4.0.0   tidyverse_2.0.0
## 
## loaded via a namespace (and not attached):
##  [1] gtable_0.3.6       jsonlite_2.0.0     compiler_4.4.3     tidyselect_1.2.1  
##  [5] jquerylib_0.1.4    scales_1.4.0       yaml_2.3.12        fastmap_1.2.0     
##  [9] here_1.0.2         R6_2.6.1           generics_0.1.4     knitr_1.51        
## [13] bookdown_0.46      rprojroot_2.1.1    bslib_0.10.0       pillar_1.11.1     
## [17] RColorBrewer_1.1-3 tzdb_0.5.0         rlang_1.1.5        stringi_1.8.4     
## [21] cachem_1.1.0       xfun_0.56          sass_0.4.10        S7_0.2.0          
## [25] otel_0.2.0         timechange_0.3.0   cli_3.6.4          withr_3.0.2       
## [29] magrittr_2.0.3     digest_0.6.37      grid_4.4.3         rstudioapi_0.18.0 
## [33] hms_1.1.4          lifecycle_1.0.5    vctrs_0.6.5        evaluate_1.0.5    
## [37] glue_1.8.0         farver_2.1.2       rmarkdown_2.30     tools_4.4.3       
## [41] pkgconfig_2.0.3    htmltools_0.5.8.1
```

``` r
rm(list =setdiff(ls(all.names = TRUE), c("markdown.dir")))
gc()
```

```
##           used (Mb) gc trigger  (Mb) max used  (Mb)
## Ncells 1280931 68.5    2518040 134.5  2122542 113.4
## Vcells 2214288 16.9    8388608  64.0  3628830  27.7
```

