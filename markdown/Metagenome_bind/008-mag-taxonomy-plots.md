---
output: 
  bookdown::html_document2:
     toc: true
---




# MAG taxonomy barplots.





## Introduction
Plot taxonomy barplots from GTDB-tk and BLAST.





## Load necessary libraries.


``` r
# install.packages(c("tidyverse","Polychrome","ggtext"))
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
library(Polychrome)
library(ggtext)
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

Add CheckM2 quality information for all MAGs:


``` r
checkm2.high_quality_mags <- read.table(checkm2.high_quality_mags.fname,
                                      sep="\t",header = T)%>%
  as_tibble()%>%
  rename("genome"="Name")%>%
  mutate(genome=gsub(paste0(metabat2.date_time,"_"),"",genome))%>%
  separate_wider_delim(genome,delim = "_bin.",names = c("sample_id","bin_id"))%>%
  mutate(bin_id=as.integer(bin_id))
```



## Taxonomy barplots (GTDB-tk) (**Figure 5**).
Prepare the dRep classification dataset for plotting: hide rare taxa,
italicise names, highlight non-bacterial species in the legend.


``` r
avg.taxon.len <- round(mean(nchar(gtdbtk.coverm.drep$classification)))
gtdbtk.coverm.drep.for_plot<-gtdbtk.coverm.drep%>%
  group_by(classification)%>%
  mutate(mean_relab=mean(relative_abundance))%>%
  ungroup()%>%
  separate_wider_delim(cols = classification,
                       delim = ";",
                       names = c("Kingdom","Phylum","Class","Order",
                                 "Family","Genus","Species"),
                       cols_remove = FALSE)%>%
  mutate(taxon = case_when(mean_relab < 0.7 ~ "Remainder (mean relative abundance < 0.7%)",
                           Species == "Fimivicinus sp900544375"~ paste0("<i>",Genus,"</i>",  " (", "<i>",Family,"</i>",")"),
                           Genus == "Phil12" ~ paste0(Genus,  " (", "<i>",Order,"</i>",")"),
                           grepl("UBA",Genus)~ paste0(Genus,  " (", "<i>",Family,"</i>",")"),
                           grepl("MGBC",Genus)~ paste0(Genus,  " (", "<i>",Family,"</i>",")"),
                           Species =="unclassified" & Genus!="unclassified" ~ paste0("<i>",Genus,"</i>",  " (", "<i>",Family,"</i>",")"),
                           Species != "unclassified" & Genus!="unclassified" ~ Species,
                           Species =="unclassified"  & Genus=="unclassified" ~ paste0("<i>",Family,"</i>")
  ),
  taxon = ifelse(! grepl("\\(", taxon) & !grepl("Remainder",taxon)& grepl(" ",taxon), 
                 paste0("<i>",taxon,"</i>"), taxon),
  taxon = case_when(Kingdom =="Archaea" & !grepl("Remainder",taxon) ~ paste0("<span style='color: violetred4'><b>",taxon,"</b></span>"),
                    Kingdom =="Bacteria" & !grepl("Remainder",taxon) ~ paste0("<span style='color: orange'><b>",taxon,"</b></span>"),
                    grepl("Remainder", taxon) ~ taxon
  ))%>%
  arrange(Family)%>%
  mutate(taxon=factor(taxon,levels=unique(taxon)))%>%
  mutate(taxon=fct_relevel(taxon,
                           "Remainder (mean relative abundance < 0.7%)",
                           after = 0))
```

Generate a color palette.


``` r
set.seed(1)
plot.cols<-createPalette(length(levels(unique(gtdbtk.coverm.drep.for_plot$taxon))),
                         seedcolors =rainbow(7))# input: number of rows

plot.cols<-c("#C1CDCD",plot.cols[1:length(plot.cols)-1])
col.vec<-setNames(plot.cols,levels(gtdbtk.coverm.drep.for_plot$taxon))
```

The barplot:


``` r
gtdbtk.coverm.drep.plot <- gtdbtk.coverm.drep.for_plot%>%
  ggplot(aes(x=sample_id,y=relative_abundance,fill=taxon))+
  geom_bar(stat="identity")+
  scale_fill_manual(values = col.vec)+ # custom fill that is based on our 
  guides(fill=guide_legend(ncol=1))+
  theme(legend.position = "right")+
  theme_bw()+
  coord_cartesian(expand=FALSE) +
  labs(x="",
       y="Relative Abundance (%)",
       fill="Classification"
       )
```

Save the plot as RDS for later.


``` r
saveRDS(gtdbtk.coverm.drep.plot, file = file.path(mag.manuscript.figures, 
                                                  "gtdbtk.coverm.drep.plot.rds"))

gtdbtk.coverm.drep.plot <- gtdbtk.coverm.drep.plot +
  mag.taxonomy.barplot.theme +
  theme(plot.title = element_text(size=15), 
        panel.spacing = unit(0.8, "cm"), # increase distance between facets
        legend.text = element_markdown(size = 13), # size of legend text
  )
```

``` r
print(gtdbtk.coverm.drep.plot + 
        ggtitle("Relative abundances of representative MAGs"))
```

<img src="008-mag-taxonomy-plots_files/figure-html/unnamed-chunk-18-1.png" alt="" width="1056" />

MetaBAT2 classification plot wasn't included in the paper, so we leave it
as it is


``` r
gtdbtk.coverm.metabat2.for_plot <- gtdbtk.coverm.metabat2%>%
  group_by(classification)%>%
  mutate(mean_relab=mean(relative_abundance))%>%
  ungroup()%>%
  # mutate(classification=ifelse(mean_relab<1,"Remainder",classification))%>%
  mutate(classification=ifelse(mean_relab<0.7,
                               "Remainder (mean relative abundance < 0.7%)",
                               classification))%>%
  mutate(classification=gsub(";","; ",classification),
         classification=str_wrap(classification,width=avg.taxon.len),
         classification=gsub("; ",";",classification))%>%
  arrange(classification)%>%
  mutate(classification=factor(classification,levels=unique(classification)))%>%
  mutate(classification=fct_relevel(classification,
                                    "Remainder (mean relative abundance < 0.7%)",
                                    after = 0))

set.seed(1)
plot.cols<-createPalette(length(levels(unique(gtdbtk.coverm.metabat2.for_plot$classification))),
                         seedcolors =rainbow(7))# input: number of rows

plot.cols<-c("#C1CDCD",plot.cols[1:length(plot.cols)-1])
col.vec<-setNames(plot.cols,levels(gtdbtk.coverm.metabat2.for_plot$classification))

gtdbtk.coverm.metabat2.plot <- gtdbtk.coverm.metabat2.for_plot%>%
  ggplot(aes(x=sample_id,y=relative_abundance,fill=classification))+
  geom_bar(stat="identity")+
  scale_fill_manual(values = col.vec)+ # custom fill that is based on our 
  guides(fill=guide_legend(ncol=1))+
  theme_bw()+
  theme(legend.position = "right")+
  # labs(title = Sys.time())+
  coord_cartesian(expand=FALSE) +
  labs(x="",
       y="Relative Abundance (%)",
       fill="Classification")+
  mag.taxonomy.barplot.theme +
  theme(plot.title = element_text(size=15), 
        panel.spacing = unit(0.8, "cm"), # increase distance between facets
        legend.text = element_text(size = 8)  # size of legend text
) 
```

``` r
print(gtdbtk.coverm.metabat2.plot +
        ggtitle("Relative abundances of all high-quality MAGs (representative + redundant)"))
```

<img src="008-mag-taxonomy-plots_files/figure-html/unnamed-chunk-20-1.png" alt="" width="1056" />



## Taxonomy barplots (BLASTN) (**Figure 5**).
Prepare the dataset in a similar manner to previous plots:


``` r
blastn.coverm.megahit.for_plot <- blastn.coverm.megahit%>%
  filter(grepl("Eukaryota",classification))%>%
  mutate(classification = str_extract(classification, "[^;]+;[^;]+;[^;]+$"))%>%
  separate_wider_delim(classification,
                       delim = ";",
                       names = c("Family","Genus", "Species"),
                       cols_remove = F)%>%
  mutate(Species = gsub(" \\(sweet potato\\)","",Species))%>%
  mutate(taxon = ifelse(grepl("sp\\.", Species),
                        paste0("<i>",Genus,"</i>",  " (", "<i>",Family,"</i>",")"),
                        paste0("<i>",Species,"</i>", " (", "<i>",Genus,"</i>",")")),
         taxon = str_wrap(taxon,width=avg.taxon.len))%>%
  mutate(is_unclassified = ifelse(grepl("sp\\.", Species),
                                  TRUE,
                                  FALSE))%>%
  arrange(is_unclassified,Species)%>%
  mutate(taxon = factor(taxon, levels = unique(taxon)))%>%
  select(-is_unclassified)

set.seed(1)
blastn.plot.cols<-createPalette(length(levels(blastn.coverm.megahit.for_plot$taxon)),
                                seedcolors =rainbow(7))# input: number of rows

# blastn.plot.cols<-c("#C1CDCD",blastn.plot.cols[1:length(plot.cols)-1])
blastn.col.vec<-setNames(blastn.plot.cols,levels(blastn.coverm.megahit.for_plot$taxon))

blastn.coverm.megahit.plot <- blastn.coverm.megahit.for_plot%>%
  ggplot(aes(x=sample_id,y=tpm,fill=taxon))+
  geom_bar(stat="identity")+
  coord_cartesian(expand = F)+
  labs(x="Sample",
       y="TPM",
       # title="TPM of eukaryotic contigs",
       fill="Classification")+
  theme_bw()+
  scale_fill_manual(values = blastn.col.vec) # custom fill that is based on col.vec 

saveRDS(blastn.coverm.megahit.plot, 
        file = file.path(mag.manuscript.figures, "blastn.coverm.megahit.plot.rds"))

blastn.coverm.megahit.plot <- blastn.coverm.megahit.plot +
  mag.taxonomy.barplot.theme +
  theme(plot.title = element_text(size = 25), # size of plot title
        legend.text = element_markdown(size = 13), # size of legend text
        legend.title = element_text(size = 15), # size of legend title
        ) # legend on the right
```

``` r
print(blastn.coverm.megahit.plot + 
        ggtitle("Classification of eukaryotic contigs with BLASTN"))
```

<img src="008-mag-taxonomy-plots_files/figure-html/unnamed-chunk-23-1.png" alt="" width="1056" />

``` r
for(image.format in image.formats){
  ggsave(filename = paste("gtdbtk-coverm-drep-barplot",
                      image.format, sep = "."),
         plot = gtdbtk.coverm.drep.plot,
         path = mag.figures,
         width=11, height=5.5,
         units="in",
         dpi=300,device = image.format)

  ggsave(filename = paste("gtdbtk-coverm-metabat2-barplot", 
                          image.format, sep = "."),
         plot=gtdbtk.coverm.metabat2.plot,
         path = mag.figures,
         width=11, height=7,units="in",
         dpi=300,device = image.format)
  ggsave(filename = paste("blastn-coverm-megahit-barplot", 
                          image.format, sep = "."),
         plot = blastn.coverm.megahit.plot,
         path = mag.figures,
         width=11, height=5,units="in",
         dpi=300,device = image.format)
}
```



## Plot MAG genome sizes. 


``` r
mag.genome.sizes.plot <- seqkit.all.df %>%
  filter(is_representative_bin == TRUE) %>%
  filter(!secondary_cluster %in% c("unassembled"))%>%
  select(sample_id,contig_id,bin_id,contig_length,secondary_cluster,GC)%>%
  group_by(secondary_cluster)%>%
  mutate(genome_size_mbp=sum(contig_length)/10^6,
         mag_gc=mean(GC),
         n_contig=n_distinct(contig_id))%>%
  ungroup() %>%
  distinct(secondary_cluster,.keep_all = T)%>%
  mutate(bin_id = as.integer(bin_id))%>%
  left_join(gtdbtk.taxonomy, by = join_by(sample_id, bin_id, secondary_cluster))%>%
  group_by(Phylum)%>%
  mutate(median_genome_size_mbp=median(genome_size_mbp))%>%
  ungroup%>%
  arrange(-median_genome_size_mbp)%>%
  select(-median_genome_size_mbp)%>%
  mutate(Phylum=factor(Phylum,levels=unique(Phylum)))%>%
  ggplot(aes(x=Phylum,y=genome_size_mbp))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()+
  labs(y="Genome size (Mbp)"#,
       # title="MAG genome sizes"
       )+
  mag.boxplot.theme +
  theme(axis.text.y = element_text(size=15), # size of y axis ticks
        plot.caption = element_text(size=13), # size of plot caption
        legend.title = element_text(size = 15) # size of legend title
  )
```

``` r
print(mag.genome.sizes.plot + 
        ggtitle ("MAG genome sizes"))
```

<img src="008-mag-taxonomy-plots_files/figure-html/unnamed-chunk-26-1.png" alt="" width="960" />

``` r
for(image.format in image.formats){
  ggsave(filename = paste("gtdbtk-mag-genome-sizes", image.format, sep = "."),
         path = mag.figures,
         plot = mag.genome.sizes.plot,
         dpi =300,
         width=10, height=7,units="in",
         # width = 4000, height = 2000,units="px",
         device = image.format)
}
```



## Plot MAG GC%. 


``` r
mag.gc.plot<- seqkit.all.df%>%
  filter(is_representative_bin ==TRUE)%>%
  select(sample_id,contig_id,bin_id,contig_length,secondary_cluster,GC)%>%
  filter(!is.na(secondary_cluster))%>%
  group_by(secondary_cluster)%>%
  mutate(genome_size_mbp=sum(contig_length)/10^6,
         mag_gc=mean(GC),
         n_contig=n_distinct(contig_id))%>%
  ungroup%>%
  distinct(secondary_cluster,.keep_all = T)%>%
  mutate(bin_id = as.integer(bin_id))%>%
  left_join(gtdbtk.taxonomy)%>%
  group_by(Phylum)%>%
  mutate(median_gc=median(mag_gc))%>%
  ungroup%>%
  arrange(-median_gc)%>%
  select(-median_gc)%>%
  mutate(Phylum=factor(Phylum,levels=unique(Phylum)))%>%
  ggplot(aes(x=Phylum,y=mag_gc))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()+
  labs(y="GC%")+
  mag.boxplot.theme + 
  theme(
    plot.margin = unit(c(0.2,0.2,0.2,0.5), "cm"),
    axis.text.y = element_text(size=13), # size of y axis ticks
    plot.caption = element_text(size=13), # size of plot caption
    legend.title = element_text(size = 15) # size of legend title
  )
```

```
## Joining with `by = join_by(sample_id, bin_id, secondary_cluster)`
```

``` r
print(mag.gc.plot + 
        ggtitle("MAG GC%"))
```

<img src="008-mag-taxonomy-plots_files/figure-html/unnamed-chunk-29-1.png" alt="" width="768" />

``` r
for(image.format in image.formats){
  ggsave(filename = paste("gtdbtk-mag-gc",image.format, sep = "."),
         plot = mag.gc.plot,
         path = mag.figures,
         dpi =300,
         width=8, height=7,units="in",
         device = image.format)
}
```



### GC% in detail.


``` r
seqkit.all.df %>%
  filter(is_representative_bin == TRUE)%>%
  select(sample_id,contig_id,bin_id,contig_length,secondary_cluster,GC)%>%
  mutate(bin_id = as.integer(bin_id))%>%
  left_join(gtdbtk.taxonomy)%>%
  filter(!is.na(secondary_cluster))%>%
  group_by(secondary_cluster)%>%
  mutate(mag_gc=mean(GC))%>%
  ungroup%>%
  distinct(secondary_cluster,.keep_all = T)%>%
  group_by(Phylum)%>%
  mutate(min_gc=min(mag_gc),
         max_gc=max(mag_gc),
         mean_gc=mean(mag_gc),
         sd_gc=sd(mag_gc),
         n=n())%>%
  ungroup%>%
  distinct(Phylum,min_gc,max_gc,mean_gc,sd_gc,n)%>%
  arrange(-mean_gc)
```

```
## Joining with `by = join_by(sample_id, bin_id, secondary_cluster)`
```

```
## # A tibble: 11 × 6
##    Phylum             min_gc max_gc mean_gc   sd_gc     n
##    <chr>               <dbl>  <dbl>   <dbl>   <dbl> <int>
##  1 Desulfobacterota     65.1   66.0    65.7  0.500      3
##  2 Pseudomonadota       47.9   65.9    54.7  9.76       3
##  3 Bacteroidota         41.3   61.8    52.8  6.76      35
##  4 Elusimicrobiota      51.6   54.0    52.8  1.72       2
##  5 Actinomycetota       35.7   63.9    45.8  9.51       8
##  6 Bacillota            27.7   69.0    44.4  9.25      57
##  7 Desulfobacterota_I   40.8   44.0    41.7  1.49       4
##  8 Spirochaetota        32.7   57.7    41.6  8.59      10
##  9 Cyanobacteriota      36.2   40.8    38.0  2.47       3
## 10 Fibrobacterota       37.0   37.0    37.0 NA          1
## 11 Methanobacteriota    28.8   29.0    28.9  0.0905     2
```

GC% in BLASTN.


``` r
blastn.coverm.megahit.eukaryotes.gc.plot <- blastn.coverm.megahit%>%
  filter(grepl("Eukaryota",classification))%>%
  left_join(seqkit.all.df[,c("sample_id","contig_id","GC")])%>%
  # filter(Kingdom=="Eukaryota")%>%
  mutate(classification=sub(".*;[^;]*;([^;]*;[^;]*)$", "\\1", classification))%>%
  separate_wider_delim(delim = ";",cols=classification,
                       names=c("Genus","Species"))%>%
  select(sample_id,contig_id,GC,Genus)%>%
  arrange(-GC)%>%
  mutate(Genus=factor(Genus,levels=unique(Genus)))%>%
  ggplot(aes(x=Genus,y=GC))+
  geom_boxplot()+
  geom_jitter()+
  theme_bw()+
  labs(y="Contig GC %")+
  mag.boxplot.theme +
  theme(
    axis.text.y = element_text(size=12), # size of y axis ticks
    plot.caption = element_text(size=10), # size of plot caption
    legend.text = element_text(size = 10), # size of legend text
    legend.title = element_text(size = 10), # size of legend title
    legend.position = "right") # legend on the right
```

```
## Joining with `by = join_by(sample_id, contig_id)`
```

``` r
print(blastn.coverm.megahit.eukaryotes.gc.plot + 
        ggtitle("GC% of eukaryotic contigs identified by BLASTN"))
```

<img src="008-mag-taxonomy-plots_files/figure-html/unnamed-chunk-33-1.png" alt="" width="768" />

``` r
for(image.format in image.formats){
  ggsave(filename = paste("blastn-coverm-megahit-eukaryotes-gc-boxplot",image.format, sep = "."),
         path = mag.figures,
         plot = blastn.coverm.megahit.eukaryotes.gc.plot,
         width = 8, height = 6,units="in",
         dpi=300,device = image.format)
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
##  [1] ggtext_0.1.2     Polychrome_1.5.4 lubridate_1.9.4  forcats_1.0.0   
##  [5] stringr_1.5.1    dplyr_1.1.4      purrr_1.0.4      readr_2.1.5     
##  [9] tidyr_1.3.1      tibble_3.2.1     ggplot2_4.0.0    tidyverse_2.0.0 
## 
## loaded via a namespace (and not attached):
##  [1] utf8_1.2.4           sass_0.4.10          generics_0.1.4      
##  [4] xml2_1.3.8           stringi_1.8.4        hms_1.1.4           
##  [7] digest_0.6.37        magrittr_2.0.3       evaluate_1.0.5      
## [10] grid_4.4.3           timechange_0.3.0     RColorBrewer_1.1-3  
## [13] bookdown_0.46        fastmap_1.2.0        rprojroot_2.1.1     
## [16] jsonlite_2.0.0       scales_1.4.0         textshaping_1.0.3   
## [19] jquerylib_0.1.4      cli_3.6.4            rlang_1.1.5         
## [22] scatterplot3d_0.3-44 litedown_0.9         commonmark_2.0.0    
## [25] withr_3.0.2          cachem_1.1.0         yaml_2.3.12         
## [28] otel_0.2.0           tools_4.4.3          tzdb_0.5.0          
## [31] colorspace_2.1-2     here_1.0.2           vctrs_0.6.5         
## [34] R6_2.6.1             lifecycle_1.0.5      ragg_1.5.0          
## [37] pkgconfig_2.0.3      pillar_1.11.1        bslib_0.10.0        
## [40] gtable_0.3.6         glue_1.8.0           Rcpp_1.0.14         
## [43] systemfonts_1.2.3    xfun_0.56            tidyselect_1.2.1    
## [46] rstudioapi_0.18.0    knitr_1.51           farver_2.1.2        
## [49] htmltools_0.5.8.1    labeling_0.4.3       rmarkdown_2.30      
## [52] compiler_4.4.3       S7_0.2.0             markdown_2.0        
## [55] gridtext_0.1.5
```

``` r
rm(list =setdiff(ls(all.names = TRUE), c("markdown.dir")))
gc()
```

```
##           used (Mb) gc trigger  (Mb) max used  (Mb)
## Ncells 1822008 97.4    3322796 177.5  3322796 177.5
## Vcells 3420301 26.1   44816464 342.0 52648669 401.7
```

