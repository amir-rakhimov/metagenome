---
output: 
  bookdown::html_document2:
     toc: true
---




# Barplots





## Introduction
After importing the kraken2 output files into R and processing them 
with phyloseq, it's time to explore the taxonomic composition of our data.
We will use the Polychrome package to create a custom palette for the 
barplots.





## Load necessary libraries.


``` r
# install.packages(c("tidyverse","ggtext","Polychrome"))
library(phyloseq)
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



## Specifying parameters and directory/file names. 


``` r
source(here::here("config/R/config.R"))# config file with global variables
source(here::here("config/R/themes.R"))# config file with themes
```

The taxonomic rank that was used for agglomeration:


``` r
agglom.rank <- "Species" 
```

Import phyloseq object as an rds file (NOT rarefied): 


``` r
ps.q <- readRDS(file= ps.q.raw.fname)
ps.q.rel <-readRDS(file= ps.q.rel.raw.fname)
custom.md <- readRDS(custom.md.path)%>%
  filter(sequencing_type =="Naked mole-rat whole metagenome sequencing")

if(agglom.rank=="OTU" | agglom.rank=="Species"){
  ps.q.agg.rel<-ps.q.rel%>%
    psmelt()
  ps.q.agg <- ps.q%>%
    psmelt()
}else{
  ps.q.agg.rel<-ps.q.rel%>%
    tax_glom(taxrank = agglom.rank)%>%
    psmelt()%>%
    select(-OTU)
  ps.q.agg<-ps.q%>%
    tax_glom(taxrank = agglom.rank)%>%
    psmelt()%>%
    select(-OTU)
}
# Join the dataframe of relative abundances with absolute abundances
ps.q.agg <-  ps.q.agg.rel %>%
  rename("RelativeAbundance" = "Abundance")%>%
  left_join(ps.q.agg[,c("Sample",agglom.rank,"Abundance")])
```

```
## Joining with `by = join_by(Sample, Species)`
```

``` r
mean_sd_relab.nmr_ages<- read.table(file.path(community.composition.tables,
                                               paste0("mean-max-min-sd-relab-nmr_ages-",agglom.rank,".tsv")),
                                    sep = "\t", header = T)
```

Check the total number of taxa (not filtered by mean relative abundance)


``` r
ps.q.agg%>%
  ungroup()%>%
  select(matches(paste0("^",agglom.rank,"$")))%>%
  pull(.)%>%
  unique()%>%
  sort()%>%
  length()
```

```
## [1] 619
```

select(matches(paste0("^",agglom.rank,"$"))) will select variables that 
match a pattern (regex).
So, if we are agglomerating by Phylum, the command will select the Phylum column




## Set up the barplot. 
Find the column number that precedes Species (Genus column).


``` r
if(agglom.rank=="Species"){
  agglom.rank.col<-which(colnames(ps.q.agg) =="Species")
  preceding.rank.col<-which(colnames(ps.q.agg) ==agglom.rank)-1
  preceding.rank<-colnames(ps.q.agg)[preceding.rank.col]
}else{
  agglom.rank.col<-which(colnames(ps.q.agg) ==agglom.rank)
  preceding.rank.col<-which(colnames(ps.q.agg) ==agglom.rank)-1
  preceding.rank<-colnames(ps.q.agg)[preceding.rank.col]
}
custom_order <- c("Kingdom", "Phylum", "Class", "Order", "Family","Genus","Species")
```

If we agglomerate by higher level (Order,Class, etc), need to adjust the rank


``` r
if(agglom.rank%in%custom_order){
  agglom.rank.index<-match(agglom.rank,custom_order)
  custom_order<-custom_order[1:agglom.rank.index-1]
}
```

Find the average number of characters in the taxa to wrap long names:


``` r
avg.taxon.len <- ps.q.agg%>%
  select(all_of(agglom.rank))%>%
  pull%>%
  nchar%>%
  mean%>%
  round
```

Hide rare taxa, wrap long names, highlight non-bacterial taxa.


``` r
new.df<-ps.q.agg%>%
  unite("taxon",Kingdom:all_of(agglom.rank),sep = ";",remove = FALSE)%>%
  mutate(is_unclassified = grepl
         ("Kingdom|Phylum|Class|Order|Family|Genus|Species",taxon))%>% # add column to show that a taxon was unclassified
  left_join(mean_sd_relab.nmr_ages) %>%
  mutate(MeanRelativeAbundance = ifelse(is.na(MeanRelativeAbundance), 0, MeanRelativeAbundance))%>%
  mutate(is_unclassified=ifelse(is_unclassified==TRUE &!grepl
                                ("Kingdom|Phylum|Class|Order|Family|Genus|Species",
                                  get(agglom.rank)),
                                FALSE, is_unclassified))%>%
  mutate(
    taxon=ifelse(is_unclassified,
                      get(agglom.rank),
                      paste0(get(agglom.rank)#, " (",get(preceding.rank),
                             #")"
                             )),
    taxon = ifelse(mapply(grepl,get(preceding.rank),get(agglom.rank)),
                   taxon,
                   paste0(get(agglom.rank), " (", get(preceding.rank),")")),
    taxon=ifelse(MeanRelativeAbundance<1,
                 "Remainder (Mean relative abundance < 1%)",
                 taxon))%>% # if a taxon is classified, we change it to show the preceding rank, e.g Family (Genus)
  mutate(taxon=gsub("_"," ",taxon),
         taxon=str_wrap(taxon,width=avg.taxon.len))%>%
  mutate(taxon=ifelse(Kingdom=="Bacteria",
                      paste0("<span style='color: orange'><b><i>",taxon,"</i></b></span>"),
                      taxon),
         taxon=ifelse(Kingdom=="Eukaryota",
                      paste0("<span style='color: palegreen4'><b><i>",taxon,"</i></b></span>"),
                      taxon),
         taxon=ifelse(Kingdom=="Archaea",
                      paste0("<span style='color: violetred4'><b><i>",taxon,"</i></b></span>"),
                      taxon),
         taxon=ifelse(grepl("Remainder",taxon),
                      "Remainder (Mean relative abundance < 1%)",taxon))%>%
  group_by_at(c("is_unclassified",preceding.rank))%>%
  arrange(desc(is_unclassified),
          get(preceding.rank),
          taxon,
          .by_group = TRUE)%>% # sort within unclassified and classified taxa (by group)
  mutate(taxon=factor(taxon,levels=unique(taxon)),
         taxon=taxon%>%
           fct_relevel(c("Remainder (Mean relative abundance < 1%)",
                         unique(unlist(lapply(custom_order,
                                              function(pat){levels(.)[str_detect(levels(.),
                                                                                 fixed(pat))]
                                              })))
           ),
           after=0))%>% # change taxon column to factor. Then, we set "Remainder" as the first level, unclassified taxa next, then classified taxa. we use the str_detect to find unclassified taxa because these have Kingdom or Phylum in the name
  ungroup()
```

```
## Joining with `by = join_by(OTU, agegroup, Kingdom, Phylum, Family, Genus,
## Species)`
```

Create a color palette.  
We need to choose colors for the taxa in our barplot. They should be
distinguishable, so we can't choose similar colors. Or at least we shouldn't
put them next to each other.

We will use `createPalette` function from the `Polychrome` package.
The function takes the number of colors for the palette and seed colors that
will be used for palette generation. In this case, we use the number of
rows in our legend (taxa) and rainbow colors (to avoid having similar colors
next to each other).
We also need to set the random seed because the output
is a bit random. The output is a vector of colors.


``` r
set.seed(1)
plot.cols<-createPalette(length(levels(unique(new.df$taxon))),
                         seedcolors =rainbow(7))# input: number of rows
```

Grey color is added ahead of others because it corresponds to Remainder 
taxa. 


``` r
plot.cols<-c("#C1CDCD",plot.cols[1:length(plot.cols)-1])
col.vec<-setNames(plot.cols,levels(new.df$taxon))
```



## Plot the barplot. 


``` r
mainplot <- new.df %>%
  left_join(custom.md[,c("Sample","age","host_birthday")])%>%
  arrange(age,host_birthday)%>%
  mutate(age = paste0("(",age,")"),
         age = factor(age),
         Sample = paste(Sample,age),
         Sample = factor(Sample,levels=unique(Sample)))%>%

  ggplot(aes(x=Sample,y=RelativeAbundance,fill=taxon))+
  geom_bar(stat="identity")+
  scale_fill_manual(values = col.vec) + # custom fill that is based on our custom palette
  coord_cartesian(expand=FALSE) +
  xlab("") +
  ylab("Relative Abundance (%)")+
  labs(fill="Taxon")+
  theme_bw() +
  guides(fill=guide_legend(ncol=1)) # legend as one column
```

```
## Joining with `by = join_by(Sample, host_birthday, age)`
```

Save the plot as an RDS object.


``` r
saveRDS(mainplot,
        file = file.path(kraken2.manuscript.figures,
                         "kraken2.mainplot.rds"))

mainplot <- mainplot +
  theme(plot.margin=unit(c(0.4,0.4,0.4,1), 'cm'),
        strip.text.x = ggtext::element_markdown(size = 10),# the name of
        # each facet will be recognised as a markdown object, so we can
        # add line breaks (cause host names are too long)
        panel.spacing = unit(0.8, "cm"), # increase distance between facets
        axis.text.x = element_text(angle=45,size=13,hjust=1),# rotate
        # the x-axis labels by 45 degrees and shift to the right
        axis.text.y = element_text(size=13), # size of y axis ticks
        axis.title = element_text(size = 15), # size of axis names
        plot.title = element_text(size = 25), # size of plot title
        plot.caption = element_text(size=23), # size of plot caption
        legend.text = element_markdown(size = 13), # size of legend text
        legend.title = element_text(size = 17), # size of legend title
        legend.position = "right") # legend on the right
```

``` r
print(mainplot+
        ggtitle(paste(agglom.rank,"level gut microbiota profiles of fecal samples from naked mole-rats"))+
        theme(plot.title = element_text(size=15))
)
```

<img src="004-barplots-kraken2_files/figure-html/unnamed-chunk-20-1.png" alt="" width="1056" />

``` r
for(image.format in image.formats){
  ggsave(filename = paste0("kraken2-barplot-NMR.", image.format),
         path = community.composition.figures,
         plot=mainplot,
         width=11, height=7,units="in",
         dpi=300,device = image.format)
}
```



## Relative abundance per kingdom (including bacteria).
Read table with unclassified read numbers.


``` r
unclassified_reads<-read.table(kraken2.unclassified.stats.fname,
                               header = T)
colnames(unclassified_reads)[which(colnames(unclassified_reads)=="Unclassified")]<-"Abundance"
unclassified_reads$Kingdom<-"Unclassified"

kingdom.plot.with_bact<-ps.q.agg%>%
  group_by(Sample,Kingdom)%>%
  summarise(Abundance=sum(Abundance))%>%
  bind_rows(unclassified_reads)%>%
  group_by(Kingdom)%>%
  summarise(TotalAbundance=sum(Abundance))%>%
  mutate(TotalAbundancePerTax=TotalAbundance/sum(TotalAbundance)*100)%>%
  arrange(-TotalAbundancePerTax)%>%
  ggplot(aes(x=reorder(Kingdom,-TotalAbundancePerTax),
             y=TotalAbundancePerTax,
             fill=Kingdom))+
  geom_bar(stat = "identity")+
  scale_fill_viridis_d(option = "C")+
  ylab("Relative abundance from all samples (%)")+
  xlab("")+
  theme_bw()
```

```
## `summarise()` has grouped output by 'Sample'. You can override using the
## `.groups` argument.
```

Save the plot as an RDS object.


``` r
saveRDS(kingdom.plot.with_bact,
        file = file.path(kraken2.manuscript.figures,
                         "kraken2.kingdom.plot.with_bact.rds"))

kingdom.plot.with_bact <- kingdom.plot.with_bact +
  coord_cartesian(expand = c("bottom" = F,
                             "left" = F,
                             "right" = F)) +
  geom_text(aes(x=Kingdom, 
                y=TotalAbundancePerTax, 
                label=round(TotalAbundancePerTax,digits = 3)),
            vjust=-0.5,size = 5)+ # add text with percentage above the bar
  theme(axis.line = element_blank(), 
        axis.text.x = element_text(size=13),
        axis.text.y = element_text(size=15), # size of y axis ticks
        axis.title = element_text(size = 15), # size of axis names
        plot.title = ggtext::element_markdown(size = 25), # the plot 
        # title will be recognised as a markdown object, so we can
        # add line breaks (cause host names are too long)
        plot.caption = element_text(size=23),# size of plot caption
        legend.text = element_text(size = 13),# size of legend text
        legend.title = element_text(size = 17), # size of legend title
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) # no legend 
```

``` r
print(kingdom.plot.with_bact+
        ggtitle("Whole metagenome sequencing profile of naked mole-rat gut microbiota
               (Kingdom level)")+
        theme(plot.title = ggtext::element_markdown(size=15))
)
```

<img src="004-barplots-kraken2_files/figure-html/unnamed-chunk-24-1.png" alt="" width="1056" />

``` r
for(image.format in image.formats){
  ggsave(filename = paste0("NMR-kingdom.plot.with_unclas.", image.format),
         path = community.composition.figures,
         plot = kingdom.plot.with_bact,
         width=11, height=6,units="in",
         dpi=300,device = image.format)
}
```

Check library size per sample (not included)


``` r
ps.q.agg%>%
  group_by(Sample,Kingdom)%>%
  summarise(Abundance=sum(Abundance))%>%
  bind_rows(unclassified_reads)%>%
  group_by(Sample)%>%
  summarise(TotalSample=sum(Abundance))%>%
  ggplot(aes(x=Sample,y=TotalSample,fill=TotalSample))+
  geom_bar(stat = "identity")
```

```
## `summarise()` has grouped output by 'Sample'. You can override using the
## `.groups` argument.
```

<img src="004-barplots-kraken2_files/figure-html/unnamed-chunk-25-1.png" alt="" width="672" />

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
##  [1] ggtext_0.1.2     Polychrome_1.5.4 lubridate_1.9.4  forcats_1.0.0   
##  [5] stringr_1.5.1    dplyr_1.1.4      purrr_1.0.4      readr_2.1.5     
##  [9] tidyr_1.3.1      tibble_3.2.1     ggplot2_4.0.0    tidyverse_2.0.0 
## [13] phyloseq_1.50.0 
## 
## loaded via a namespace (and not attached):
##  [1] ade4_1.7-23             tidyselect_1.2.1        viridisLite_0.4.3      
##  [4] farver_2.1.2            Biostrings_2.74.1       S7_0.2.0               
##  [7] fastmap_1.2.0           digest_0.6.37           timechange_0.3.0       
## [10] lifecycle_1.0.5         cluster_2.1.8.1         survival_3.8-3         
## [13] magrittr_2.0.3          compiler_4.4.3          rlang_1.1.5            
## [16] sass_0.4.10             tools_4.4.3             igraph_2.1.4           
## [19] yaml_2.3.12             data.table_1.17.8       knitr_1.51             
## [22] labeling_0.4.3          scatterplot3d_0.3-44    here_1.0.2             
## [25] xml2_1.3.8              plyr_1.8.9              RColorBrewer_1.1-3     
## [28] withr_3.0.2             BiocGenerics_0.52.0     grid_4.4.3             
## [31] stats4_4.4.3            multtest_2.62.0         biomformat_1.34.0      
## [34] colorspace_2.1-2        Rhdf5lib_1.28.0         scales_1.4.0           
## [37] iterators_1.0.14        MASS_7.3-65             cli_3.6.4              
## [40] rmarkdown_2.30          vegan_2.6-4             crayon_1.5.3           
## [43] ragg_1.5.0              generics_0.1.4          otel_0.2.0             
## [46] rstudioapi_0.18.0       httr_1.4.7              reshape2_1.4.4         
## [49] tzdb_0.5.0              commonmark_2.0.0        ape_5.8-1              
## [52] cachem_1.1.0            rhdf5_2.50.2            zlibbioc_1.52.0        
## [55] splines_4.4.3           parallel_4.4.3          XVector_0.46.0         
## [58] vctrs_0.6.5             Matrix_1.7-4            jsonlite_2.0.0         
## [61] litedown_0.9            bookdown_0.46           IRanges_2.40.1         
## [64] hms_1.1.4               S4Vectors_0.44.0        systemfonts_1.2.3      
## [67] foreach_1.5.2           jquerylib_0.1.4         glue_1.8.0             
## [70] codetools_0.2-20        stringi_1.8.4           gtable_0.3.6           
## [73] GenomeInfoDb_1.42.3     UCSC.utils_1.2.0        pillar_1.11.1          
## [76] htmltools_0.5.8.1       rhdf5filters_1.18.1     GenomeInfoDbData_1.2.13
## [79] R6_2.6.1                textshaping_1.0.3       rprojroot_2.1.1        
## [82] evaluate_1.0.5          lattice_0.22-6          Biobase_2.66.0         
## [85] markdown_2.0            gridtext_0.1.5          bslib_0.10.0           
## [88] Rcpp_1.0.14             nlme_3.1-167            permute_0.9-8          
## [91] mgcv_1.9-1              xfun_0.56               pkgconfig_2.0.3
```

``` r
rm(list =setdiff(ls(all.names = TRUE), c("markdown.dir")))
gc()
```

```
##           used  (Mb) gc trigger  (Mb) max used  (Mb)
## Ncells 5547312 296.3    8349707 446.0  8349707 446.0
## Vcells 9572532  73.1   18551710 141.6 15393092 117.5
```

