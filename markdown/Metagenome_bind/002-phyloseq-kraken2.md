---
output: 
  bookdown::html_document2:
     toc: true
---






# Convert Kraken2 output into a phyloseq object





## Introduction
In this script, we will import the dataset from Kraken2 using 
phyloseq.





## Load necessary libraries


``` r
# install.packages(c("tidyverse"))
# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
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
library(microViz)
```

```
## microViz version 0.12.1 - Copyright (C) 2021-2024 David Barnett
## ! Website: https://david-barnett.github.io/microViz
## ✔ Useful?  For citation details, run: `citation("microViz")`
## ✖ Silence? `suppressPackageStartupMessages(library(microViz))`
```

``` r
library(phyloseq)

source(here::here("config/R/config.R"))# config file with global variables
```



## Import the combined report from Kraken2.


``` r
combined.report<-read.table(combined.report.filename, 
            header = T, 
            sep = "\t",
            fill=TRUE,
            quote = "", 
            comment.char = "@",
            na.strings = c("","NA"))%>%
  as_tibble()
```

Convert into wide format:


``` r
combined.report<-combined.report%>%
  filter(Abundance!=0)%>%
  pivot_wider(names_from = Sample,
              values_from = Abundance)

dim(combined.report)
```

```
## [1] 620  18
```

``` r
head(combined.report)%>%
  knitr::kable(format = "simple")
```



Kingdom       Phylum         Class           Order               Family               Genus             Species                              2D10      2D14       G14       G18      H15       H21        H3        H4      O15      Y51b      Y66b
------------  -------------  --------------  ------------------  -------------------  ----------------  -------------------------------  --------  --------  --------  --------  -------  --------  --------  --------  -------  --------  --------
k__Bacteria   p__Bacillota   c__Clostridia   o__Lachnospirales   f__Lachnospiraceae   g__Anaerostipes   s__Anaerostipes_caccae            2380645   5980486   2272260   1279006    47290   2314957   2442271   2032729   552419   1929362   3049341
k__Bacteria   p__Bacillota   c__Clostridia   o__Lachnospirales   f__Lachnospiraceae   g__Anaerostipes   s__Anaerostipes_rhamnosivorans      10750     15291      8617      4430      561      8448      8456      8231     2400      6845      9031
k__Bacteria   p__Bacillota   c__Clostridia   o__Lachnospirales   f__Lachnospiraceae   g__Anaerostipes   s__Anaerostipes_hadrus              23445      1854    653563      6801    17232        NA      1592   1732918    12775     44529      4063
k__Bacteria   p__Bacillota   c__Clostridia   o__Lachnospirales   f__Lachnospiraceae   g__Blautia        s__Blautia_parvula                 608654   1209122    759221    454882   314325    140405    127138    342366   197539    284080    336088
k__Bacteria   p__Bacillota   c__Clostridia   o__Lachnospirales   f__Lachnospiraceae   g__Blautia        s__Blautia_pseudococcoides          42287     17924      5985     35452     6177      1874      1689      2048     1582      5933      4083
k__Bacteria   p__Bacillota   c__Clostridia   o__Lachnospirales   f__Lachnospiraceae   g__Blautia        s__Blautia_sp._SC05B48               5524      5646      4047     15088     9221      5171      3298     15408     3273      5284      4974

``` r
# all.ranks<-c("Kingdom", "Phylum", "Class", "Order", "Family","Genus","Species")
# combined.report[,all.ranks][is.na(combined.report[,all.ranks])]<-"NA"
# combined.report[,-which(names(combined.report)%in%all.ranks)][is.na(combined.report[,-which(names(combined.report)%in%all.ranks)])]<-0
```

Extract table with Species abundances where rows are Species and columns 
are samples. Here, we combine all the columns into one and call it "OTU".


``` r
df.otus<-combined.report%>%
  unite("OTU",Kingdom:Species,sep = "|")%>%
  column_to_rownames("OTU")

df.otus[is.na(df.otus)]<-0
```

Extract table with taxonomy where rows are OTUS and columns are taxonomy ranks:
Kingdom, ..., Species


``` r
df.taxa<-combined.report%>%
  select(Kingdom:Species)%>%
  unite("OTU",Kingdom:Species,sep = "|",remove = FALSE)%>%
  column_to_rownames("OTU")

df.taxa$Kingdom<-
  gsub("k__","",df.taxa$Kingdom)
df.taxa$Kingdom<-
  gsub("d__","",df.taxa$Kingdom)
df.taxa$Phylum<-
  gsub("d__|p__","",df.taxa$Phylum)
df.taxa$Class<-
  gsub("d__|p__|c__","",df.taxa$Class)
df.taxa$Order<-
  gsub("d__|p__|c__|o__","",df.taxa$Order)
df.taxa$Family<-
  gsub("d__|p__|c__|o__|f__","",df.taxa$Family)
df.taxa$Genus<-
  gsub("d__|p__|c__|o__|f__|g__","",df.taxa$Genus)
df.taxa$Species<-
  gsub("d__|p__|c__|o__|f__|g__|s__","",df.taxa$Species)
head(df.taxa)%>%
  knitr::kable(format = "simple")
```

                                                                                                                             Kingdom    Phylum      Class        Order            Family            Genus          Species                     
---------------------------------------------------------------------------------------------------------------------------  ---------  ----------  -----------  ---------------  ----------------  -------------  ----------------------------
k__Bacteria|p__Bacillota|c__Clostridia|o__Lachnospirales|f__Lachnospiraceae|g__Anaerostipes|s__Anaerostipes_caccae           Bacteria   Bacillota   Clostridia   Lachnospirales   Lachnospiraceae   Anaerostipes   Anaerostipes_caccae         
k__Bacteria|p__Bacillota|c__Clostridia|o__Lachnospirales|f__Lachnospiraceae|g__Anaerostipes|s__Anaerostipes_rhamnosivorans   Bacteria   Bacillota   Clostridia   Lachnospirales   Lachnospiraceae   Anaerostipes   Anaerostipes_rhamnosivorans 
k__Bacteria|p__Bacillota|c__Clostridia|o__Lachnospirales|f__Lachnospiraceae|g__Anaerostipes|s__Anaerostipes_hadrus           Bacteria   Bacillota   Clostridia   Lachnospirales   Lachnospiraceae   Anaerostipes   Anaerostipes_hadrus         
k__Bacteria|p__Bacillota|c__Clostridia|o__Lachnospirales|f__Lachnospiraceae|g__Blautia|s__Blautia_parvula                    Bacteria   Bacillota   Clostridia   Lachnospirales   Lachnospiraceae   Blautia        Blautia_parvula             
k__Bacteria|p__Bacillota|c__Clostridia|o__Lachnospirales|f__Lachnospiraceae|g__Blautia|s__Blautia_pseudococcoides            Bacteria   Bacillota   Clostridia   Lachnospirales   Lachnospiraceae   Blautia        Blautia_pseudococcoides     
k__Bacteria|p__Bacillota|c__Clostridia|o__Lachnospirales|f__Lachnospiraceae|g__Blautia|s__Blautia_sp._SC05B48                Bacteria   Bacillota   Clostridia   Lachnospirales   Lachnospiraceae   Blautia        Blautia_sp._SC05B48         



## Convert into phyloseq taxonomyTable object.


``` r
df.taxa<-tax_table(as.matrix(df.taxa))

df.taxa.pretty<-tax_fix(df.taxa)
df.taxa<-df.taxa.pretty
rm(df.taxa.pretty)
```

Add metadata


``` r
custom.md <- readRDS(custom.md.path)%>%
  filter(sequencing_type =="Naked mole-rat whole metagenome sequencing")%>%
  remove_rownames()%>%
  filter(Sample%in%colnames(combined.report))%>%
  column_to_rownames("Sample")
```



## Create a phyloseq object. 


``` r
ps.q<-phyloseq(otu_table(df.otus,taxa_are_rows = TRUE),
               tax_table(df.taxa),
               sample_data(custom.md))
```

Remove human data 


``` r
ps.q <- subset_taxa(ps.q,Phylum !="Chordata")
```

Number of features in the unfiltered dataset:


``` r
ntaxa(ps.q)
```

```
## [1] 619
```

Total frequency in the unfiltered dataset:


``` r
sum(sample_sums(ps.q))
```

```
## [1] 430634351
```

Summary statistics (min, median, max, quartiles) of the unfiltered dataset:


``` r
summary(sample_sums(ps.q))
```

```
##     Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
## 34359903 36618876 39032318 39148577 41260087 44778833
```



## 5. Calculate relative abundance: Abundance divided by total abundance in a sample. ####


``` r
ps.q.rel <- transform_sample_counts(ps.q, function(x) 100*x/sum(x)) 
```

Save phyloseq objects and the contents (OTU table, taxonomy, tree).


``` r
if(!file.exists(ps.q.raw.fname)){
  saveRDS(ps.q,
          file = ps.q.raw.fname)
}

if(!file.exists(ps.q.rel.raw.fname)){
  saveRDS(ps.q.rel,
          file = ps.q.rel.raw.fname)
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
##  [1] phyloseq_1.50.0 microViz_0.12.1 lubridate_1.9.4 forcats_1.0.0  
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
## [19] igraph_2.1.4            yaml_2.3.12             data.table_1.17.8      
## [22] knitr_1.51              here_1.0.2              plyr_1.8.9             
## [25] RColorBrewer_1.1-3      withr_3.0.2             BiocGenerics_0.52.0    
## [28] grid_4.4.3              stats4_4.4.3            multtest_2.62.0        
## [31] biomformat_1.34.0       colorspace_2.1-2        Rhdf5lib_1.28.0        
## [34] scales_1.4.0            iterators_1.0.14        MASS_7.3-65            
## [37] cli_3.6.4               rmarkdown_2.30          vegan_2.6-4            
## [40] crayon_1.5.3            generics_0.1.4          otel_0.2.0             
## [43] rstudioapi_0.18.0       httr_1.4.7              reshape2_1.4.4         
## [46] tzdb_0.5.0              ape_5.8-1               cachem_1.1.0           
## [49] rhdf5_2.50.2            zlibbioc_1.52.0         splines_4.4.3          
## [52] parallel_4.4.3          XVector_0.46.0          vctrs_0.6.5            
## [55] Matrix_1.7-4            jsonlite_2.0.0          bookdown_0.46          
## [58] IRanges_2.40.1          hms_1.1.4               S4Vectors_0.44.0       
## [61] foreach_1.5.2           jquerylib_0.1.4         glue_1.8.0             
## [64] codetools_0.2-20        stringi_1.8.4           gtable_0.3.6           
## [67] GenomeInfoDb_1.42.3     UCSC.utils_1.2.0        pillar_1.11.1          
## [70] htmltools_0.5.8.1       rhdf5filters_1.18.1     GenomeInfoDbData_1.2.13
## [73] R6_2.6.1                rprojroot_2.1.1         evaluate_1.0.5         
## [76] lattice_0.22-6          Biobase_2.66.0          bslib_0.10.0           
## [79] Rcpp_1.0.14             nlme_3.1-167            permute_0.9-8          
## [82] mgcv_1.9-1              xfun_0.56               pkgconfig_2.0.3
```

``` r
rm(list =setdiff(ls(all.names = TRUE), c("markdown.dir")))
gc()
```

```
##           used  (Mb) gc trigger  (Mb) max used  (Mb)
## Ncells 5408896 288.9    8818158 471.0  7907204 422.3
## Vcells 9303353  71.0   16139409 123.2 13379673 102.1
```

