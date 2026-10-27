---
output: 
  bookdown::html_document2:
     toc: true
---




# Gene annotation analysis.





## Introduction
Import the data from gene annotation and convert into tibble format.





## Load necessary libraries.


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



## Import datasets.


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

Tidy taxonomy table


``` r
gtdbtk.taxonomy <- readRDS(gtdbtk.taxonomy.clean.fname.rds)
```



## Analyse PROKKA mapping between samples, contigs, and genes.


``` r
prokka.contig.gene.map.df <- data.table::fread(prokka.contig.gene.map.fname,
                                      header = FALSE,sep="\t")%>%
  as_tibble()
colnames(prokka.contig.gene.map.df)<-c("sample_id","contig_id",
                                       "locus_tag","gene_length","gene_type")
```

3,847,253 genes (including redundancies):


``` r
nrow(prokka.contig.gene.map.df)
```

```
## [1] 3847253
```

``` r
head(prokka.contig.gene.map.df) %>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## sample_id   contig_id     locus_tag         gene_length  gene_type 
## ----------  ------------  ---------------  ------------  ----------
## 2D10        k141_0        GPCOLHEP_00001            141  CDS       
## 2D10        k141_141795   GPCOLHEP_00002            444  CDS       
## 2D10        k141_171643   GPCOLHEP_00003            789  CDS       
## 2D10        k141_111946   GPCOLHEP_00004            315  CDS       
## 2D10        k141_14926    GPCOLHEP_00005            414  CDS       
## 2D10        k141_67167    GPCOLHEP_00006           1017  CDS
```

1,633,598 contigs had some genes from PROKKA (no matter how small or 
whether they're not CDS):


``` r
prokka.contig.gene.map.df%>%
  distinct(sample_id,contig_id)%>%
  nrow
```

```
## [1] 1633598
```

3771950 CDS in total:


``` r
prokka.contig.gene.map.df<-prokka.contig.gene.map.df%>%
  filter(gene_type=="CDS")
prokka.contig.gene.map.df%>%
  distinct(locus_tag)%>%
  nrow
```

```
## [1] 3771950
```

3764752 CDS >=100 bp:


``` r
prokka.contig.gene.map.df<-prokka.contig.gene.map.df%>%
  filter(gene_type=="CDS",gene_length>=100)
prokka.contig.gene.map.df%>%
  distinct(locus_tag)%>%
  nrow
```

```
## [1] 3764752
```

1,618,027 contigs had CDS that were >=100 bp long


``` r
prokka.contig.gene.map.df%>%
  distinct(sample_id,contig_id)%>%
  nrow
```

```
## [1] 1618027
```



## Add gene count data from HTseq. Then, normalise to TPM.
Sources: https://metagenomics-workshop.readthedocs.io/en/latest/annotation/quantification.html  
https://github.com/EnvGen/metagenomics-workshop/blob/master/in-house/tpm_table.py


``` r
htseq.gene_counts.df<-data.table::fread(htseq.gene_counts.fname,
                                 header = FALSE,sep="\t")%>%
  as_tibble()
colnames(htseq.gene_counts.df)<-c("locus_tag","raw_count")
```

3,847,253 rows, like in prokka.contig.gene.map.df


``` r
nrow(htseq.gene_counts.df)
```

```
## [1] 3847253
```

``` r
head(htseq.gene_counts.df) %>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## locus_tag         raw_count
## ---------------  ----------
## GPCOLHEP_00001            1
## GPCOLHEP_00002            5
## GPCOLHEP_00003           11
## GPCOLHEP_00004            0
## GPCOLHEP_00005            7
## GPCOLHEP_00006           10
```



## MMseqs2 data has mapping between all locus_tag and representative clusters.


``` r
mmseqs.clusters.df <- data.table::fread(mmseqs.clusters.fname,
                               header = FALSE,sep = "\t")%>%
  as_tibble()
colnames(mmseqs.clusters.df)<-c("locus_tag_cluster","locus_tag")
head(mmseqs.clusters.df) %>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## locus_tag_cluster   locus_tag       
## ------------------  ----------------
## GPCOLHEP_142781     GPCOLHEP_142781 
## GPCOLHEP_142781     BPLGMAPI_87195  
## GPCOLHEP_142781     LBKEGPMD_376473 
## GPCOLHEP_142781     HGNGMENJ_233808 
## GPCOLHEP_142781     HFIMOJPH_168376 
## GPCOLHEP_142781     NNLNMLNN_159525
```

3,764,752 rows


``` r
nrow(mmseqs.clusters.df)
```

```
## [1] 3764752
```

1,922,857 unique clusters


``` r
mmseqs.clusters.df%>%
  distinct(locus_tag_cluster)
```

```
## # A tibble: 1,922,857 × 1
##    locus_tag_cluster
##    <chr>            
##  1 GPCOLHEP_142781  
##  2 GPCOLHEP_142782  
##  3 GPCOLHEP_142783  
##  4 GPCOLHEP_142784  
##  5 GPCOLHEP_142785  
##  6 GPCOLHEP_142787  
##  7 GPCOLHEP_142788  
##  8 GPCOLHEP_142790  
##  9 GPCOLHEP_142791  
## 10 GPCOLHEP_142793  
## # ℹ 1,922,847 more rows
```

Map locus_tag_clusters to samples: 3847253 entries (all mappings, 
including redundancies).  
Add gene counts: raw counts, so we need to normalise.


``` r
gene.annotation.df <- mmseqs.clusters.df%>%
  full_join(prokka.contig.gene.map.df,by = join_by(locus_tag))
```

Most genes are short:


``` r
gene.annotation.df%>%
  filter(gene_type=="CDS",!is.na(locus_tag_cluster))%>%
  distinct(locus_tag_cluster,.keep_all = T)%>%
  summary()
```

```
##  locus_tag_cluster   locus_tag          sample_id          contig_id        
##  Length:1922857     Length:1922857     Length:1922857     Length:1922857    
##  Class :character   Class :character   Class :character   Class :character  
##  Mode  :character   Mode  :character   Mode  :character   Mode  :character  
##                                                                             
##                                                                             
##                                                                             
##   gene_length       gene_type        
##  Min.   :  102.0   Length:1922857    
##  1st Qu.:  252.0   Class :character  
##  Median :  405.0   Mode  :character  
##  Mean   :  582.5                     
##  3rd Qu.:  714.0                     
##  Max.   :74937.0
```



## TPM normalisation.
Combine locus_tag with their counts. Since we have locus_tag_cluster info,
we know the abundance for clusters is the sum of all locus_tag 
counts in the cluster. We're working with CDS >=100 bp long.


``` r
gene.annotation.df <- gene.annotation.df%>%
  full_join(htseq.gene_counts.df,by = join_by(locus_tag))%>%
  filter(!is.na(locus_tag_cluster))
head(gene.annotation.df) %>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## locus_tag_cluster   locus_tag         sample_id   contig_id      gene_length  gene_type    raw_count
## ------------------  ----------------  ----------  ------------  ------------  ----------  ----------
## GPCOLHEP_142781     GPCOLHEP_142781   2D10        k141_48106            1260  CDS                108
## GPCOLHEP_142781     BPLGMAPI_87195    2D14        k141_181676           1260  CDS                 39
## GPCOLHEP_142781     LBKEGPMD_376473   G18         k141_157071           1260  CDS                 42
## GPCOLHEP_142781     HGNGMENJ_233808   H15         k141_153424           1260  CDS                 37
## GPCOLHEP_142781     HFIMOJPH_168376   H3          k141_143288           1260  CDS                 51
## GPCOLHEP_142781     NNLNMLNN_159525   H4          k141_216927           1260  CDS                 41
```

Find representative gene lengths (clustered or unclustered gene length):


``` r
representative.gene_lengths <- gene.annotation.df%>%
  group_by(locus_tag_cluster)%>%
  select(locus_tag_cluster,locus_tag,gene_length)%>%
  summarise(rep_gene_len=gene_length[locus_tag == locus_tag_cluster][1])
head(representative.gene_lengths) %>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## locus_tag_cluster    rep_gene_len
## ------------------  -------------
## BPLGMAPI_00001                204
## BPLGMAPI_00002                249
## BPLGMAPI_00006                171
## BPLGMAPI_00007                192
## BPLGMAPI_00008                219
## BPLGMAPI_00009                411
```

Sum counts within each sample:


``` r
tpm.df <- gene.annotation.df%>%
  group_by(sample_id,locus_tag_cluster)%>%
  summarise(locus_tag_cluster_sum_count=sum(raw_count))%>%
  ungroup%>%
  # Add gene_count*read_length/representative_gene_length
  left_join(representative.gene_lengths)%>%
  mutate(rgxrl_over_gene_len=(locus_tag_cluster_sum_count*150)/rep_gene_len)%>%
  # Get the final TPM
  group_by(sample_id)%>%
  mutate(T_denom=sum(rgxrl_over_gene_len),
         tpm=(locus_tag_cluster_sum_count*150*10^6) / (rep_gene_len * T_denom))%>%
  ungroup
```

```
## `summarise()` has grouped output by 'sample_id'. You can override using the
## `.groups` argument.
## Joining with `by = join_by(locus_tag_cluster)`
```

``` r
head(tpm.df) %>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## sample_id   locus_tag_cluster    locus_tag_cluster_sum_count   rep_gene_len   rgxrl_over_gene_len   T_denom         tpm
## ----------  ------------------  ----------------------------  -------------  --------------------  --------  ----------
## 2D10        BPLGMAPI_00010                                 3            126             3.5714286   1757091   2.0325808
## 2D10        BPLGMAPI_00054                                 4            843             0.7117438   1757091   0.4050695
## 2D10        BPLGMAPI_00103                                 6            552             1.6304348   1757091   0.9279173
## 2D10        BPLGMAPI_00105                                 6            144             6.2500000   1757091   3.5570164
## 2D10        BPLGMAPI_00115                                12            429             4.1958042   1757091   2.3879271
## 2D10        BPLGMAPI_00119                                23            597             5.7788945   1757091   3.2888995
```

``` r
rm(prokka.contig.gene.map.df)
rm(mmseqs.clusters.df)
rm(htseq.gene_counts.df)
gc()
```

```
##            used  (Mb) gc trigger   (Mb)  max used   (Mb)
## Ncells  5584245 298.3   19736153 1054.1  24670191 1317.6
## Vcells 83968288 640.7  189170612 1443.3 189077322 1442.6
```



## Add KofamScan annotation.


``` r
kofamscan.df <- read.table(kofamscan.fname,header = TRUE,sep = "\t")%>%
  as_tibble()
colnames(kofamscan.df)[1] <- "locus_tag_cluster"
head(kofamscan.df) %>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## locus_tag_cluster   ko        threshold   score   evalue  ko_definition                                             
## ------------------  -------  ----------  ------  -------  ----------------------------------------------------------
## BPLGMAPI_00013      K00765        40.77   238.4        0  ATP phosphoribosyltransferase [EC:2.4.2.17]               
## BPLGMAPI_00027      K06960        47.93    97.5        0  uncharacterized protein                                   
## BPLGMAPI_00029      K11145        36.47   156.1        0  mini-ribonuclease III [EC:3.1.26.-]                       
## BPLGMAPI_00044      K14742        84.77   129.6        0  tRNA threonylcarbamoyladenosine biosynthesis protein TsaB 
## BPLGMAPI_00053      K06907       143.97   297.3        0  Bacteriophage tail sheath protein                         
## BPLGMAPI_00054      K00858        36.23   194.0        0  NAD+ kinase [EC:2.7.1.23]
```

Not all mmseqs clusters are in kofamscan.df, but all kofamscan.df are 
in clusters (obviously)


``` r
table(unique(gene.annotation.df$locus_tag_cluster)%in%kofamscan.df$locus_tag_cluster)
```

```
## 
##   FALSE    TRUE 
## 1546944  375913
```

``` r
table(kofamscan.df$locus_tag_cluster%in%unique(gene.annotation.df$locus_tag_cluster))
```

```
## 
##   TRUE 
## 375913
```



## Add dbcan annotation.


``` r
dbcan.df <- read.table(dbcan.fname,header=T,sep = "\t",
                     stringsAsFactors = FALSE,comment.char = "")%>%
  as_tibble()
colnames(dbcan.df)<-c("locus_tag_cluster","ec","hmmer","dbcan_sub","diamond","n_of_tools")
head(dbcan.df) %>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## locus_tag_cluster   ec           hmmer            dbcan_sub   diamond    n_of_tools
## ------------------  -----------  ---------------  ----------  --------  -----------
## BPLGMAPI_00136      -            GT4(292-431)     GT4_e1362   -                   2
## BPLGMAPI_00138      -            GT2(6-133)       GT2         -                   2
## BPLGMAPI_00168      -            GH57(20-326)     GH57_e10    GH57                3
## BPLGMAPI_00169      -            GT4(240-384)     GT4_e602    GT4                 3
## BPLGMAPI_00170      -            GH133(285-636)   GH133_e1    GH133               3
## BPLGMAPI_01197      3.2.1.23:1   GH2(21-587)      GH2_e112    GH2                 3
```

Not all mmseqs clusters are in dbcan, but all dbcan are in clusters (obviously)


``` r
table(unique(gene.annotation.df$locus_tag_cluster)%in%dbcan.df$locus_tag_cluster)
```

```
## 
##   FALSE    TRUE 
## 1878209   44648
```

``` r
table(dbcan.df$locus_tag_cluster%in%unique(gene.annotation.df$locus_tag_cluster))
```

```
## 
##  TRUE 
## 44648
```

16752 CAZymes were identified by 1 tool


``` r
dbcan.df%>%
  count(n_of_tools)
```

```
## # A tibble: 3 × 2
##   n_of_tools     n
##        <int> <int>
## 1          1 16752
## 2          2  9769
## 3          3 18127
```

Filter for confident CAZyme hits (≥2 tools agree): count how many tools
gave a hit per gene:


``` r
dbcan.df <- dbcan.df %>%
  # By default, use hmmer. But if hmmer is unclassified, use dbcan_sub
  mutate(caz_with_aa_coord = ifelse(hmmer != "-", hmmer, dbcan_sub),
         caz_with_aa_coord = ifelse(hmmer=="-" & dbcan_sub == "-", diamond, caz_with_aa_coord))%>%
  # Remove empty entries ('-') for tools that didn't identify a locus_tag_cluster
  select(-ec,-n_of_tools,-hmmer,-dbcan_sub,-diamond)%>%
  distinct(locus_tag_cluster, caz_with_aa_coord, .keep_all = T)  %>% # Remove duplicates
  # Remove text in brackets and CBM (e.g. "GH5(1-100)+CBM50(110-180)")
  mutate(caz_subclass_combination = gsub("_[0-9]*|_e[0-9]*","", caz_with_aa_coord),
         caz_subclass_combination = gsub("\\([0-9]+\\-[0-9]+\\)","", caz_subclass_combination),
         caz_subclass_combination = gsub("in[A-Z]+[0-9]+","", caz_subclass_combination)
         )%>%
  # Split multi-domain strings (e.g. "GH5+CBM50")
  mutate(caz_subclass = caz_subclass_combination)%>%
  separate_longer_delim(caz_subclass, delim = "+") %>%
  mutate(caz_class=ifelse(grepl("[A-Za-z]+[0-9]*$",caz_subclass),
                          gsub("[0-9]*$","",caz_subclass),caz_subclass))
```

Extract caz_subclass name: at the start of the string, match any character
except ( for one or more times (that's why we add + sign).
Matches a sequence of one or more characters until the first () from the
beginning of the string.  
```mutate(caz_subclass = str_extract(caz, "^[^\\(]+"),  
       caz_subclass=gsub("_[0-9]*|_e[0-9]*","",caz_subclass))```
CAZ class is caz subclass without numbers at the end. For example, GH1 becomes GH.


``` r
head(dbcan.df)%>%
  knitr::kable(format = "simple")%>%
  print()
```

```
## 
## 
## locus_tag_cluster   caz_with_aa_coord   caz_subclass_combination   caz_subclass   caz_class 
## ------------------  ------------------  -------------------------  -------------  ----------
## BPLGMAPI_00136      GT4(292-431)        GT4                        GT4            GT        
## BPLGMAPI_00138      GT2(6-133)          GT2                        GT2            GT        
## BPLGMAPI_00168      GH57(20-326)        GH57                       GH57           GH        
## BPLGMAPI_00169      GT4(240-384)        GT4                        GT4            GT        
## BPLGMAPI_00170      GH133(285-636)      GH133                      GH133          GH        
## BPLGMAPI_01197      GH2(21-587)         GH2                        GH2            GH
```

``` r
dbcan.df%>%
  nrow()
```

```
## [1] 46238
```

24063 CAZymes


``` r
dbcan.df%>%
  summarise(n_caz=n_distinct(caz_with_aa_coord),
            n_caz_subclass_combination=n_distinct(caz_subclass_combination),
            n_caz_subclass=n_distinct(caz_subclass),
            n_caz_class=n_distinct(caz_class))
```

```
## # A tibble: 1 × 4
##   n_caz n_caz_subclass_combination n_caz_subclass n_caz_class
##   <int>                      <int>          <int>       <int>
## 1 24063                        553            294          20
```



## Map dbCAN and kofamscan locus_tag_clusters with MMseqs locus_tag_clusters and samples.


``` r
# Then, add columns that show whether a locus_tag_cluster was classified by dbcan
# or by kofamscan, or both
# gene.annotation.df <- gene.annotation.df %>%
#   full_join(dbcan.df,by = join_by(locus_tag_cluster),
#             relationship = "many-to-many")%>%
#   full_join(kofamscan.df[,c("locus_tag_cluster","ko")],
#             by = join_by(locus_tag_cluster))%>%
#   mutate(annotation_type=ifelse(is.na(caz)&is.na(ko),"None","Annotated"),
#          annotation_type=ifelse(!is.na(caz)&is.na(ko),"dbCAN",annotation_type),
#          annotation_type=ifelse(is.na(caz)&!is.na(ko),"KofamScan",annotation_type),
#          annotation_type=ifelse(!is.na(caz)&!is.na(ko),"dbCAN+KofamScan",annotation_type))
# head(gene.annotation.df)
# rm(dbcan.df)
# rm(kofamscan.df)
gc()
```

```
##            used  (Mb) gc trigger   (Mb)  max used   (Mb)
## Ncells  5632768 300.9   15788923  843.3  24670191 1317.6
## Vcells 86577526 660.6  189170612 1443.3 189077322 1442.6
```

Save data:


``` r
if(!file.exists(gene.annotation.df.fname.rds) & !file.exists(gene.annotation.df.fname.tsv)){
  saveRDS(gene.annotation.df, file = gene.annotation.df.fname.rds)
  data.table::fwrite(gene.annotation.df, file = gene.annotation.df.fname.tsv,
                     quote = TRUE,sep = '\t',row.names = FALSE, col.names = TRUE)
}

if(!file.exists(dbcan.df.fname.rds) & !file.exists(dbcan.df.fname.tsv)){
  saveRDS(dbcan.df, file = dbcan.df.fname.rds)
  data.table::fwrite(dbcan.df, file = dbcan.df.fname.tsv,
                     quote = TRUE,sep = '\t',row.names =FALSE,col.names = TRUE)
}

if(!file.exists(kofamscan.df.fname.rds) & !file.exists(kofamscan.df.fname.tsv)){
  saveRDS(kofamscan.df, file = kofamscan.df.fname.rds)
  data.table::fwrite(kofamscan.df, file = kofamscan.df.fname.tsv,
                     quote = TRUE,sep = '\t',row.names = FALSE,col.names = TRUE)
}

if(!file.exists(representative.gene_lengths.fname.rds) & !file.exists(representative.gene_lengths.fname.tsv)){
  saveRDS(representative.gene_lengths, file = representative.gene_lengths.fname.rds)
  data.table::fwrite(representative.gene_lengths, file = representative.gene_lengths.fname.tsv,
                     quote = TRUE,sep = '\t',row.names = FALSE,col.names = TRUE)
}
if(!file.exists(tpm.df.fname.rds) & !file.exists(tpm.df.fname.tsv)){
  saveRDS(tpm.df, file = tpm.df.fname.rds)
  data.table::fwrite(tpm.df, file = tpm.df.fname.tsv,
                     quote = TRUE,sep = '\t',row.names = FALSE,col.names = TRUE)
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
## [41] knitr_1.51         farver_2.1.2       htmltools_0.5.8.1  rmarkdown_2.30    
## [45] compiler_4.4.3     S7_0.2.0           gridtext_0.1.5
```

``` r
rm(list =setdiff(ls(all.names = TRUE), c("markdown.dir")))
gc()
```

```
##           used (Mb) gc trigger   (Mb)  max used   (Mb)
## Ncells 1463421 78.2   12631139  674.6  24670191 1317.6
## Vcells 4607883 35.2  151336490 1154.7 189077322 1442.6
```

