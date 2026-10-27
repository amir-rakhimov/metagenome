---
output: 
  bookdown::html_document2:
     toc: true
---






# Analysing Kraken2 data with phyloseq. 





## Introduction
In this script, we will explore the imported dataset from Kraken2 using 
phyloseq.

We will use the data from 002-phyloseq-kraken2.R script (ps.q.agg
agglomerated tables at phylum, family, genus, and Species level).




## Load necessary libraries and scripts.


``` r
# install.packages(c("tidyverse","ggtext"))
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
library(ggtext)
library(vegan)
```

```
## Loading required package: permute
## Loading required package: lattice
## This is vegan 2.6-4
```



## Import the config file.


``` r
source(here::here("config/R/config.R"))# config file with global variables
source(here::here("config/R/themes.R"))# config file with themes
```

Load necessary scripts.


``` r
source(file.path(util.functions.r,"get_mean_and_sd_relab_from_phyloseq.R"))
source(file.path(util.functions.r,"create_summary_stats.R"))
source(file.path(util.functions.r,"multi_rarefy_phyloseq.R"))
source(file.path(util.functions.r,"get_average_summary_stats.R"))

source(file.path(util.functions.r,"get_unclassified_summary_stats.R"))
source(file.path(util.functions.r,"get_average_unclassified_summary_stats.R"))
source(file.path(util.functions.r,"ggplot_species.R"))
```



Import raw datasets as rds files.


``` r
ps.q <- readRDS(file = ps.q.raw.fname)
ps.q.rel <- readRDS(file = ps.q.rel.raw.fname)
custom.md <- readRDS(custom.md.path)%>%
  filter(sequencing_type =="Naked mole-rat whole metagenome sequencing")
gg.labs.name <- "Age group"
```



## Calculating summary statistics.

Load the NMR age metadata.


``` r
custom.md.ages<-readRDS(file.path(metadatadir,"custom.md.ages.rds"))%>%
  filter(sequencing_type == "Naked mole-rat whole metagenome sequencing")
pretty.level.names <-names(table(custom.md$agegroup))
names(pretty.level.names)<-names(table(custom.md$agegroup))
```

Setup sample levels for NMR for barplots.


``` r
sample.levels<-custom.md.ages%>%
  ungroup()%>%
  select(Sample,age)%>%
  arrange(age)%>%
  distinct()%>%
  mutate(NewSample=paste0(Sample," (",age,")"))%>%
  mutate(Sample=factor(Sample,levels=Sample),
         NewSample=factor(NewSample,levels=NewSample))
```

Agglomerate the `ps.q.filtered` object and convert it into a dataframe with `psmelt()` 
function. Then, save it in tab-separated format and as an rds object. 


``` r
melt_phyloseq_and_save <-function (ps, ps.index){
  ps.indexes <- c("Kingdom","Phylum","Class","Order", "Family", "Genus","Species")
  
  if(ps.index == "OTU"){
    ps.df <- ps %>%
      psmelt()
  }else{
    ps.df <- ps %>%
      tax_glom(ps.index, NArm=FALSE)%>%
      psmelt()
  }
  
  if(ps.index == "OTU"){
    ps.df <- ps.df %>%
      dplyr::select(all_of(c("OTU","Sample", "Abundance", "class", ps.indexes)))
  }else{
    selected.indexes <- ps.indexes[1:match(ps.index, ps.indexes)] # vector of taxonomic ranks to select
    ps.df <- ps.df %>%
      dplyr::select(all_of(c("Sample", "Abundance", "class", selected.indexes)))
  }
  # Save the tables in TSV format and as an RDS object
  ps.df.fname <- case_when(ps.index =="OTU" ~ ps.q.agg.asv.fname.no_ext,
                           ps.index =="Species" ~ ps.q.agg.species.fname.no_ext,
                           ps.index =="Genus" ~ ps.q.agg.genus.fname.no_ext,
                           ps.index =="Family" ~ps.q.agg.family.fname.no_ext,
                           ps.index =="Phylum" ~ ps.q.agg.phylum.fname.no_ext
  ) 
  tsv.fname <- file.path(community.composition.tables, paste0(ps.df.fname,".tsv"))
  rda.fname <- file.path(community.composition.rdafiles, paste0(ps.df.fname,".rds"))
  if(!file.exists(tsv.fname)){
    write.table(ps.df,
                file = tsv.fname,
                row.names = F,sep = "\t")
  }
  if(!file.exists(rda.fname)){
    saveRDS(ps.df,
            file = rda.fname)
  }
  return(ps.df)
}

ps.q.agg.asv<-melt_phyloseq_and_save(ps.q, ps.index = "OTU")%>%
  filter(Abundance!=0)%>%
  as_tibble()
ps.q.rel.agg.asv<- ps.q.rel %>%
  psmelt() %>%
  dplyr::select(all_of(c("OTU","Sample", "Abundance", "class", 
                         c("Kingdom","Phylum","Class","Order", "Family", "Genus","Species"))))%>%
  filter(Abundance!=0)%>%
  as_tibble()%>%
  rename("RelativeAbundance" = "Abundance")
ps.q.agg.species<-melt_phyloseq_and_save(ps.q, ps.index = "Species")%>%
  filter(Abundance!=0)%>%
  as_tibble()
ps.q.agg.genus <- melt_phyloseq_and_save(ps.q, ps.index = "Genus")%>%
  filter(Abundance!=0)%>%
  as_tibble()
ps.q.agg.family <- melt_phyloseq_and_save(ps.q, ps.index = "Family")%>%
  filter(Abundance!=0)%>%
  as_tibble()
ps.q.agg.phylum<-melt_phyloseq_and_save(ps.q, ps.index = "Phylum")%>%
  filter(Abundance!=0)%>%
  as_tibble()
```

Save filtered phyloseq objects and the contents (OTU table, taxonomy, tree).


``` r
if(!file.exists(ps.q.raw.fname)){
  saveRDS(ps.q,
          file = ps.q.raw.fname)
}

if(!file.exists(ps.q.rel.raw.fname)){
  saveRDS(ps.q.rel,
          file = ps.q.rel.raw.fname)
}

if(!file.exists(ps.q.raw.otu_table.fname)){
  saveRDS(ps.q@otu_table %>% as.matrix() %>%as.data.frame(),
          file = ps.q.raw.otu_table.fname)
}
if(!file.exists(ps.q.raw.tax_table.fname)){
  saveRDS(ps.q@tax_table%>% as.matrix() %>%as.data.frame(),
          file = ps.q.raw.tax_table.fname)
}
ps.q.raw.otu_table.fname.tsv <- 
  file.path(community.composition.tables, 
            gsub("rds", "tsv",basename(ps.q.raw.otu_table.fname))
)
if(!file.exists(ps.q.raw.otu_table.fname.tsv)){
  write.table(ps.q@otu_table%>% as.matrix() %>%as.data.frame() %>% rownames_to_column("OTU"),
              file = ps.q.raw.otu_table.fname.tsv,
              sep = "\t", 
              row.names = FALSE)
  
}
ps.q.raw.tax_table.fname.tsv <- 
  file.path(community.composition.tables, 
            gsub("rds", "tsv",basename(ps.q.raw.tax_table.fname)))
if(!file.exists(ps.q.raw.tax_table.fname.tsv)){
  write.table(ps.q@tax_table%>% as.matrix() %>%as.data.frame()%>% rownames_to_column("OTU"),
              file = ps.q.raw.tax_table.fname.tsv,
              sep = "\t", 
              row.names = FALSE)
  
}
```



### Calculating summary statistics per host.




The columns are:  
- Total reads  
- Library size (mean Abundance ± SD)  
- Number of ASV per host  
- Number of phyla per host  
- Number of families per host  
- Number of genera per host  
Here, Species are called ASV because the same function is used in the 
16S rRNA gene sequencing analysis.  


``` r
summary.stats <- create_summary_stats(ps.q.agg.asv)
knitr::kable(as.data.frame(summary.stats),format = "simple")
```



class    TotalSamplesPerHost   TotalReadsPerHost   MeanLibrarySize   SDLibrarySize   ASVPerHost   PhylaPerHost   FamiliesPerHost   GeneraPerHost
------  --------------------  ------------------  ----------------  --------------  -----------  -------------  ----------------  --------------
NMR                       11           430634351          39148577         3492314          619             34               185             343

``` r
summary.stats.fname<-file.path(community.composition.tables,
                               "summary-table.tsv")
if(!file.exists(summary.stats.fname)){
  write.table(summary.stats,
              file=summary.stats.fname,
              row.names = F,sep = "\t")
}
```



### Create a summary table for each sample.


``` r
n.taxa.per_sample<-ps.q.agg.asv%>%
  group_by(Sample)%>%
  summarise(n_species=n_distinct(Species),
            n_genera=n_distinct(Genus),
            n_families=n_distinct(Family),
            n_phyla=n_distinct(Phylum))
knitr::kable(n.taxa.per_sample,format = "simple")
```



Sample    n_species   n_genera   n_families   n_phyla
-------  ----------  ---------  -----------  --------
2D10            287        161           83        16
2D14            321        183          102        23
G14             308        177           91        19
G18             346        186          102        25
H15             365        203          114        22
H21             319        181          102        23
H3              329        192          104        22
H4              317        174          101        21
O15             356        188          103        23
Y51b            322        195          102        21
Y66b            351        207          108        21

``` r
n.taxa.per_sample.fname<-file.path(community.composition.tables,
                               "kraken2-n_taxa-per-sample-table.tsv")
if(!file.exists(n.taxa.per_sample.fname)){
  write.table(n.taxa.per_sample,
              file=n.taxa.per_sample.fname,
              row.names = F,sep = "\t")
}
```



## Calculate average, min, max, and SD relative abundances from phyloseq object.
For phyla, families, and genera, no grouping. For species 
(here, we call them OTU), create a dataframe by age group and total (no grouping).


``` r
mean_sd_relab.nmr_no_groups.phyla <- get_mean_and_sd_relab_from_phyloseq(ps.q.rel, "class", "Phylum")
mean_sd_relab.nmr_no_groups.families <- get_mean_and_sd_relab_from_phyloseq(ps.q.rel, "class", "Family")
mean_sd_relab.nmr_no_groups.genera <- get_mean_and_sd_relab_from_phyloseq(ps.q.rel, "class", "Genus")
```

By age group.


``` r
mean_sd_relab.nmr_ages.asv <- 
  get_mean_and_sd_relab_from_phyloseq(ps.q.rel, "agegroup", "OTU")
```

In total (no grouping).


``` r
mean_sd_relab.nmr_no_groups.asv <- 
  get_mean_and_sd_relab_from_phyloseq(ps.q.rel, "class", "OTU")

knitr::kable(head(mean_sd_relab.nmr_no_groups.phyla),format = "simple")
```



class   Kingdom     Phylum            MeanRelativeAbundance   sdRelativeAbundance   minRelativeAbundance   maxRelativeAbundance   n_samples
------  ----------  ---------------  ----------------------  --------------------  ---------------------  ---------------------  ----------
NMR     Bacteria    Bacillota                     40.682066            12.2238057             20.8003012              60.132657          11
NMR     Eukaryota   Streptophyta                  23.363548             3.2931779             18.3233954              29.118435          11
NMR     Bacteria    Bacteroidota                  22.856495             9.8706711              7.5881662              41.637210          11
NMR     Bacteria    Actinomycetota                 4.870373             1.2582191              3.1767290               8.269697          11
NMR     Bacteria    Pseudomonadota                 3.670396             0.3030072              3.1864896               4.249930          11
NMR     Archaea     Euryarchaeota                  2.490599             3.3789754              0.0265766              12.065008          11

``` r
knitr::kable(head(mean_sd_relab.nmr_no_groups.families),format = "simple")
```



class   Kingdom     Phylum         Family                MeanRelativeAbundance   sdRelativeAbundance   minRelativeAbundance   maxRelativeAbundance   n_samples
------  ----------  -------------  -------------------  ----------------------  --------------------  ---------------------  ---------------------  ----------
NMR     Bacteria    Bacteroidota   Bacteroidaceae                    16.011616              8.441047               5.060626               30.10842          11
NMR     Bacteria    Bacillota      Lachnospiraceae                   14.376708              5.547365               4.885386               26.02276          11
NMR     Bacteria    Bacillota      Enterococcaceae                    7.427814              2.656070               3.408864               11.09731          11
NMR     Bacteria    Bacteroidota   Tannerellaceae                     5.368870              3.116431               0.905413               11.02561          11
NMR     Bacteria    Bacillota      Acidaminococcaceae                 5.159374             11.496351               0.000000               34.31199           8
NMR     Eukaryota   Streptophyta   Fabaceae                           4.998068              2.979583               1.860610               10.74997          11

``` r
knitr::kable(head(mean_sd_relab.nmr_no_groups.genera),format = "simple")
```



class   Kingdom     Phylum         Family               Genus                    MeanRelativeAbundance   sdRelativeAbundance   minRelativeAbundance   maxRelativeAbundance   n_samples
------  ----------  -------------  -------------------  ----------------------  ----------------------  --------------------  ---------------------  ---------------------  ----------
NMR     Bacteria    Bacteroidota   Bacteroidaceae       Bacteroides                          15.975296              8.432028              5.0466012               30.04948          11
NMR     Bacteria    Bacillota      Enterococcaceae      Enterococcus                          7.425959              2.655508              3.4088639               11.09731          11
NMR     Bacteria    Bacillota      Lachnospiraceae      Anaerostipes                          6.250690              3.879774              0.1485352               14.12265          11
NMR     Bacteria    Bacteroidota   Tannerellaceae       Parabacteroides                       5.367006              3.116180              0.9049122               11.02302          11
NMR     Bacteria    Bacillota      Acidaminococcaceae   Phascolarctobacterium                 5.159374             11.496351              0.0000000               34.31199           8
NMR     Eukaryota   Streptophyta   Convolvulaceae       Ipomoea                               4.394285              2.679296              1.7091821               10.52671          11

``` r
knitr::kable(head(mean_sd_relab.nmr_ages.asv),format = "simple")
```



agegroup       Kingdom    Phylum         Family                  Genus            Species                            OTU                                                                                                                                      MeanRelativeAbundance   sdRelativeAbundance   minRelativeAbundance   maxRelativeAbundance   n_samples
-------------  ---------  -------------  ----------------------  ---------------  ---------------------------------  --------------------------------------------------------------------------------------------------------------------------------------  ----------------------  --------------------  ---------------------  ---------------------  ----------
agegroup0_10   Bacteria   Bacillota      Enterococcaceae         Enterococcus     Enterococcus_faecium               k__Bacteria|p__Bacillota|c__Bacilli|o__Lactobacillales|f__Enterococcaceae|g__Enterococcus|s__Enterococcus_faecium                                     8.916690              1.481384              6.9493215              10.539180           4
agegroup0_10   Bacteria   Bacillota      Lachnospiraceae         Anaerostipes     Anaerostipes_caccae                k__Bacteria|p__Bacillota|c__Clostridia|o__Lachnospirales|f__Lachnospiraceae|g__Anaerostipes|s__Anaerostipes_caccae                                    8.636014              3.925037              4.9429860              14.082276           4
agegroup0_10   Bacteria   Bacillota      Peptostreptococcaceae   Clostridioides   Clostridioides_difficile           k__Bacteria|p__Bacillota|c__Clostridia|o__Peptostreptococcales|f__Peptostreptococcaceae|g__Clostridioides|s__Clostridioides_difficile                 6.289186              2.746359              2.8370679               9.550888           4
agegroup0_10   Bacteria   Bacteroidota   Bacteroidaceae          Bacteroides      Bacteroides_thetaiotaomicron       k__Bacteria|p__Bacteroidota|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides_thetaiotaomicron                          4.954287              3.998748              0.3264551              10.061379           4
agegroup0_10   Bacteria   Bacteroidota   Bacteroidaceae          Bacteroides      Bacteroides_xylanisolvens          k__Bacteria|p__Bacteroidota|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides_xylanisolvens                             4.390999              3.408020              0.0288592               7.456517           4
agegroup0_10   Bacteria   Bacillota      Eubacteriales Order     Intestinimonas   Intestinimonas_butyriciproducens   k__Bacteria|p__Bacillota|c__Clostridia|o__Eubacteriales|NA|g__Intestinimonas|s__Intestinimonas_butyriciproducens                                      3.554697              1.880238              1.4144947               5.527708           4

``` r
knitr::kable(head(mean_sd_relab.nmr_no_groups.asv),format = "simple")
```



class   Kingdom     Phylum         Family               Genus                   Species                         OTU                                                                                                                                               MeanRelativeAbundance   sdRelativeAbundance   minRelativeAbundance   maxRelativeAbundance   n_samples
------  ----------  -------------  -------------------  ----------------------  ------------------------------  -----------------------------------------------------------------------------------------------------------------------------------------------  ----------------------  --------------------  ---------------------  ---------------------  ----------
NMR     Bacteria    Bacteroidota   Bacteroidaceae       Bacteroides             Bacteroides_xylanisolvens       k__Bacteria|p__Bacteroidota|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides_xylanisolvens                                      8.647850              8.297227              0.0288592               24.80239          11
NMR     Bacteria    Bacillota      Enterococcaceae      Enterococcus            Enterococcus_faecium            k__Bacteria|p__Bacillota|c__Bacilli|o__Lactobacillales|f__Enterococcaceae|g__Enterococcus|s__Enterococcus_faecium                                              7.080774              2.698159              3.1880824               10.72238          11
NMR     Bacteria    Bacillota      Lachnospiraceae      Anaerostipes            Anaerostipes_caccae             k__Bacteria|p__Bacillota|c__Clostridia|o__Lachnospirales|f__Lachnospiraceae|g__Anaerostipes|s__Anaerostipes_caccae                                             5.654391              3.674711              0.1079273               14.08228          11
NMR     Bacteria    Bacillota      Acidaminococcaceae   Phascolarctobacterium   Phascolarctobacterium_faecium   k__Bacteria|p__Bacillota|c__Negativicutes|o__Acidaminococcales|f__Acidaminococcaceae|g__Phascolarctobacterium|s__Phascolarctobacterium_faecium                 5.155998             11.497636              0.0000000               34.31199           3
NMR     Bacteria    Bacteroidota   Bacteroidaceae       Bacteroides             Bacteroides_thetaiotaomicron    k__Bacteria|p__Bacteroidota|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides_thetaiotaomicron                                   4.563492              4.594072              0.3264551               15.08105          11
NMR     Eukaryota   Streptophyta   Convolvulaceae       Ipomoea                 Ipomoea_triloba                 k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Solanales|f__Convolvulaceae|g__Ipomoea|s__Ipomoea_triloba                                                     4.394285              2.679296              1.7091821               10.52671          11

Save tables.


``` r
mean_sd_relab.nmr_no_groups.phyla.fname <- 
  file.path(community.composition.tables,
            "mean-max-min-sd-relab-nmr_no_groups-Phylum.tsv")
if(!file.exists(mean_sd_relab.nmr_no_groups.phyla.fname)){
  write.table(mean_sd_relab.nmr_no_groups.phyla,
              file = mean_sd_relab.nmr_no_groups.phyla.fname,
              row.names = F,sep = "\t")
}

mean_sd_relab.nmr_no_groups.families.fname <-
  file.path(community.composition.tables,
            "mean-max-min-sd-relab-nmr_no_groups-Family.tsv")
if(!file.exists (mean_sd_relab.nmr_no_groups.families.fname)){
  write.table(mean_sd_relab.nmr_no_groups.families,
              file = mean_sd_relab.nmr_no_groups.families.fname,
              row.names = F,sep = "\t")
}

mean_sd_relab.nmr_no_groups.genera.fname <-
  file.path(community.composition.tables,
            "mean-max-min-sd-relab-nmr_no_groups-Genus.tsv")
if(! file.exists (mean_sd_relab.nmr_no_groups.genera.fname)){
  write.table(mean_sd_relab.nmr_no_groups.genera,
              file = mean_sd_relab.nmr_no_groups.genera.fname,
              row.names = F,sep = "\t")
}

mean_sd_relab.nmr_ages.asv.fname <-
  file.path(community.composition.tables,
            "mean-max-min-sd-relab-nmr_ages-Species.tsv")
if(! file.exists (mean_sd_relab.nmr_ages.asv.fname)){
  write.table(mean_sd_relab.nmr_ages.asv,
              file = mean_sd_relab.nmr_ages.asv.fname,
              row.names = F,sep = "\t")
}
mean_sd_relab.nmr_no_groups.asv.fname <-
  file.path(community.composition.tables,
            "mean-max-min-sd-relab-nmr_no_groups-Species.tsv")
if(! file.exists (mean_sd_relab.nmr_no_groups.asv.fname)){
  write.table(mean_sd_relab.nmr_no_groups.asv,
              file = mean_sd_relab.nmr_no_groups.asv.fname,
              row.names = F,sep = "\t")
}
```



## Rarefy the table and check the percentage of unclassified taxa. ####
First, perform multiple rarefaction:


``` r
# print(paste("Rarefying data with", rare.num_samples, "subsamples"))
# rarefied.fname.all <- file.path(community.composition.rdafiles,paste0(
#   paste("ps.q.rare","OTU",paste0(rare.num_samples,"_iter"),
#         paste("NMR",collapse = '-'),sep = "-"),".rds"))
# if(!file.exists(rarefied.fname.all)){
#   ps.q.agg.asv.nmr.rare<-multi_rarefy_phyloseq(ps.q,
#                                                niter = rare.num_samples,
#                                                output.filename.rds =  rarefied.fname.all
#                                                )
# }else{
#   ps.q.agg.asv.nmr.rare<- readRDS(rarefied.fname.all)
# }
```

#+ echo=FALSE
### 6.3 Create a summary stats table for the rarefied dataframe. ####
#' 
#' ### Create a summary stats table for the rarefied dataframe.
summary.stats.rare<-get_average_summary_stats(rare.list.asv = ps.q.agg.asv.nmr.rare)
knitr::kable(summary.stats.rare, format ="simple")
summary.stats.rare.fname<-file.path(community.composition.tables,
                                    "summary-table-rarefied.tsv")
if(!file.exists(summary.stats.rare.fname)){
  write.table(summary.stats.rare,
              file=summary.stats.rare.fname,
              row.names = F,sep = "\t")
}






``` r
## Check the most abundant bacterial phyla, families, genera in NMR.
```

Phyla:


``` r
mean_sd_relab.nmr_no_groups.phyla %>%
  filter(Kingdom == "Bacteria")%>%
  head()%>%
  knitr::kable(format = "simple")
```



class   Kingdom    Phylum             MeanRelativeAbundance   sdRelativeAbundance   minRelativeAbundance   maxRelativeAbundance   n_samples
------  ---------  ----------------  ----------------------  --------------------  ---------------------  ---------------------  ----------
NMR     Bacteria   Bacillota                     40.6820655            12.2238057             20.8003012             60.1326565          11
NMR     Bacteria   Bacteroidota                  22.8564945             9.8706711              7.5881662             41.6372098          11
NMR     Bacteria   Actinomycetota                 4.8703733             1.2582191              3.1767290              8.2696969          11
NMR     Bacteria   Pseudomonadota                 3.6703964             0.3030072              3.1864896              4.2499305          11
NMR     Bacteria   Fibrobacterota                 0.1651987             0.1604877              0.0017196              0.4889666          11
NMR     Bacteria   Elusimicrobiota                0.0712514             0.1641982              0.0000000              0.5600452           8

Families:


``` r
mean_sd_relab.nmr_no_groups.families %>%
  filter(Kingdom == "Bacteria")%>%
  head()%>%
  knitr::kable(format = "simple")
```



class   Kingdom    Phylum         Family                   MeanRelativeAbundance   sdRelativeAbundance   minRelativeAbundance   maxRelativeAbundance   n_samples
------  ---------  -------------  ----------------------  ----------------------  --------------------  ---------------------  ---------------------  ----------
NMR     Bacteria   Bacteroidota   Bacteroidaceae                       16.011616              8.441047              5.0606261              30.108419          11
NMR     Bacteria   Bacillota      Lachnospiraceae                      14.376708              5.547365              4.8853859              26.022758          11
NMR     Bacteria   Bacillota      Enterococcaceae                       7.427814              2.656070              3.4088639              11.097310          11
NMR     Bacteria   Bacteroidota   Tannerellaceae                        5.368870              3.116431              0.9054130              11.025607          11
NMR     Bacteria   Bacillota      Acidaminococcaceae                    5.159374             11.496351              0.0000000              34.311994           8
NMR     Bacteria   Bacillota      Peptostreptococcaceae                 3.560652              2.994174              0.3648371               9.550888          11

Genera:


``` r
mean_sd_relab.nmr_no_groups.genera %>%
  filter(Kingdom == "Bacteria")%>%
  head()%>%
  knitr::kable(format = "simple")
```



class   Kingdom    Phylum         Family                  Genus                    MeanRelativeAbundance   sdRelativeAbundance   minRelativeAbundance   maxRelativeAbundance   n_samples
------  ---------  -------------  ----------------------  ----------------------  ----------------------  --------------------  ---------------------  ---------------------  ----------
NMR     Bacteria   Bacteroidota   Bacteroidaceae          Bacteroides                          15.975296              8.432028              5.0466012              30.049475          11
NMR     Bacteria   Bacillota      Enterococcaceae         Enterococcus                          7.425959              2.655508              3.4088639              11.097310          11
NMR     Bacteria   Bacillota      Lachnospiraceae         Anaerostipes                          6.250690              3.879774              0.1485352              14.122648          11
NMR     Bacteria   Bacteroidota   Tannerellaceae          Parabacteroides                       5.367006              3.116180              0.9049122              11.023018          11
NMR     Bacteria   Bacillota      Acidaminococcaceae      Phascolarctobacterium                 5.159374             11.496351              0.0000000              34.311994           8
NMR     Bacteria   Bacillota      Peptostreptococcaceae   Clostridioides                        3.560652              2.994174              0.3648371               9.550888          11

Species:


``` r
mean_sd_relab.nmr_no_groups.asv %>%
  filter(Kingdom == "Bacteria")%>%
  head()%>%
  knitr::kable(format = "simple")
```



class   Kingdom    Phylum         Family               Genus                   Species                         OTU                                                                                                                                               MeanRelativeAbundance   sdRelativeAbundance   minRelativeAbundance   maxRelativeAbundance   n_samples
------  ---------  -------------  -------------------  ----------------------  ------------------------------  -----------------------------------------------------------------------------------------------------------------------------------------------  ----------------------  --------------------  ---------------------  ---------------------  ----------
NMR     Bacteria   Bacteroidota   Bacteroidaceae       Bacteroides             Bacteroides_xylanisolvens       k__Bacteria|p__Bacteroidota|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides_xylanisolvens                                      8.647850              8.297227              0.0288592               24.80239          11
NMR     Bacteria   Bacillota      Enterococcaceae      Enterococcus            Enterococcus_faecium            k__Bacteria|p__Bacillota|c__Bacilli|o__Lactobacillales|f__Enterococcaceae|g__Enterococcus|s__Enterococcus_faecium                                              7.080774              2.698159              3.1880824               10.72238          11
NMR     Bacteria   Bacillota      Lachnospiraceae      Anaerostipes            Anaerostipes_caccae             k__Bacteria|p__Bacillota|c__Clostridia|o__Lachnospirales|f__Lachnospiraceae|g__Anaerostipes|s__Anaerostipes_caccae                                             5.654391              3.674711              0.1079273               14.08228          11
NMR     Bacteria   Bacillota      Acidaminococcaceae   Phascolarctobacterium   Phascolarctobacterium_faecium   k__Bacteria|p__Bacillota|c__Negativicutes|o__Acidaminococcales|f__Acidaminococcaceae|g__Phascolarctobacterium|s__Phascolarctobacterium_faecium                 5.155998             11.497636              0.0000000               34.31199           3
NMR     Bacteria   Bacteroidota   Bacteroidaceae       Bacteroides             Bacteroides_thetaiotaomicron    k__Bacteria|p__Bacteroidota|c__Bacteroidia|o__Bacteroidales|f__Bacteroidaceae|g__Bacteroides|s__Bacteroides_thetaiotaomicron                                   4.563492              4.594072              0.3264551               15.08105          11
NMR     Bacteria   Bacteroidota   Tannerellaceae       Parabacteroides         Parabacteroides_distasonis      k__Bacteria|p__Bacteroidota|c__Bacteroidia|o__Bacteroidales|f__Tannerellaceae|g__Parabacteroides|s__Parabacteroides_distasonis                                 3.783470              2.833600              0.5081927               10.24094          11



## Check the most abundant non-bacterial phyla, families, genera in NMR.
Non-bacterial phyla:


``` r
mean_sd_relab.nmr_no_groups.phyla.nonbact<-mean_sd_relab.nmr_no_groups.phyla%>%
  filter(Kingdom != "Bacteria")
head(mean_sd_relab.nmr_no_groups.phyla.nonbact)%>%
  knitr::kable(format = "simple")
```



class   Kingdom     Phylum           MeanRelativeAbundance   sdRelativeAbundance   minRelativeAbundance   maxRelativeAbundance   n_samples
------  ----------  --------------  ----------------------  --------------------  ---------------------  ---------------------  ----------
NMR     Eukaryota   Streptophyta                23.3635478             3.2931779             18.3233954             29.1184346          11
NMR     Archaea     Euryarchaeota                2.4905986             3.3789754              0.0265766             12.0650075          11
NMR     Eukaryota   Ascomycota                   0.9540174             0.8389948              0.1280452              2.4811722          11
NMR     Eukaryota   Parabasalia                  0.6564723             0.5321252              0.1328118              1.7745074          11
NMR     Eukaryota   Apicomplexa                  0.0254256             0.0599956              0.0000000              0.1924194           3
NMR     Eukaryota   Basidiomycota                0.0149041             0.0122868              0.0000000              0.0285213           8

Non-bacterial families:


``` r
mean_sd_relab.nmr_no_groups.families.nonbact<-mean_sd_relab.nmr_no_groups.families%>%
  filter(Kingdom != "Bacteria")
head(mean_sd_relab.nmr_no_groups.families.nonbact)%>%
  knitr::kable(format = "simple")
```



class   Kingdom     Phylum          Family                 MeanRelativeAbundance   sdRelativeAbundance   minRelativeAbundance   maxRelativeAbundance   n_samples
------  ----------  --------------  --------------------  ----------------------  --------------------  ---------------------  ---------------------  ----------
NMR     Eukaryota   Streptophyta    Fabaceae                            4.998068             2.9795829              1.8606097              10.749966          11
NMR     Eukaryota   Streptophyta    Convolvulaceae                      4.394285             2.6792958              1.7091821              10.526715          11
NMR     Eukaryota   Streptophyta    Cupressaceae                        2.609274             0.3424485              2.3218895               3.390330          11
NMR     Archaea     Euryarchaeota   Methanobacteriaceae                 2.489866             3.3790568              0.0261897              12.065008          11
NMR     Eukaryota   Streptophyta    Proteaceae                          2.112845             0.9487269              0.9430488               3.620957          11
NMR     Eukaryota   Streptophyta    Salicaceae                          1.643107             0.8245515              0.3688666               2.984799          11

Non-bacterial genera:


``` r
mean_sd_relab.nmr_no_groups.genera.nonbact<-mean_sd_relab.nmr_no_groups.genera%>%
  filter(Kingdom != "Bacteria")
head(mean_sd_relab.nmr_no_groups.genera.nonbact)%>%
  knitr::kable(format = "simple")
```



class   Kingdom     Phylum          Family                Genus                 MeanRelativeAbundance   sdRelativeAbundance   minRelativeAbundance   maxRelativeAbundance   n_samples
------  ----------  --------------  --------------------  -------------------  ----------------------  --------------------  ---------------------  ---------------------  ----------
NMR     Eukaryota   Streptophyta    Convolvulaceae        Ipomoea                            4.394285             2.6792958              1.7091821              10.526715          11
NMR     Eukaryota   Streptophyta    Cupressaceae          Cryptomeria                        2.609274             0.3424485              2.3218895               3.390330          11
NMR     Archaea     Euryarchaeota   Methanobacteriaceae   Methanobrevibacter                 2.482948             3.3575025              0.0261897              11.988910          11
NMR     Eukaryota   Streptophyta    Proteaceae            Macadamia                          2.067874             0.9723804              0.8940231               3.620957          11
NMR     Eukaryota   Streptophyta    Salicaceae            Populus                            1.643107             0.8245515              0.3688666               2.984799          11
NMR     Eukaryota   Streptophyta    Fabaceae              Glycine                            1.210820             2.1318113              0.0000000               5.712605           5

Non-bacterial species:


``` r
mean_sd_relab.nmr_no_groups.asv.nonbact<-mean_sd_relab.nmr_no_groups.asv%>%
  filter(Kingdom != "Bacteria")
head(mean_sd_relab.nmr_no_groups.asv.nonbact)%>%
  knitr::kable(format = "simple")
```



class   Kingdom     Phylum          Family                Genus                Species                       OTU                                                                                                                                                 MeanRelativeAbundance   sdRelativeAbundance   minRelativeAbundance   maxRelativeAbundance   n_samples
------  ----------  --------------  --------------------  -------------------  ----------------------------  -------------------------------------------------------------------------------------------------------------------------------------------------  ----------------------  --------------------  ---------------------  ---------------------  ----------
NMR     Eukaryota   Streptophyta    Convolvulaceae        Ipomoea              Ipomoea_triloba               k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Solanales|f__Convolvulaceae|g__Ipomoea|s__Ipomoea_triloba                                                      4.3942845             2.6792958              1.7091821              10.526715          11
NMR     Eukaryota   Streptophyta    Cupressaceae          Cryptomeria          Cryptomeria_japonica          k__Eukaryota|p__Streptophyta|c__Pinopsida|o__Cupressales|f__Cupressaceae|g__Cryptomeria|s__Cryptomeria_japonica                                                 2.6092737             0.3424485              2.3218895               3.390330          11
NMR     Eukaryota   Streptophyta    Proteaceae            Macadamia            Macadamia_integrifolia        k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Proteales|f__Proteaceae|g__Macadamia|s__Macadamia_integrifolia                                                 2.0678735             0.9723804              0.8940231               3.620957          11
NMR     Archaea     Euryarchaeota   Methanobacteriaceae   Methanobrevibacter   Methanobrevibacter_sp._AbM4   k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanobrevibacter|s__Methanobrevibacter_sp._AbM4                1.2808634             1.6381718              0.0261897               5.897323          11
NMR     Eukaryota   Streptophyta    Fabaceae              Glycine              Glycine_soja                  k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Fabales|f__Fabaceae|g__Glycine|s__Glycine_soja                                                                 1.2108200             2.1318113              0.0000000               5.712605           5
NMR     Eukaryota   Streptophyta    Salicaceae            Populus              Populus_trichocarpa           k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Malpighiales|f__Salicaceae|g__Populus|s__Populus_trichocarpa                                                   0.9575214             0.7136691              0.3688666               2.478445          11

``` r
mean_sd_relab.nmr_no_groups.phyla.nonbact.fname <-file.path(community.composition.tables,
                                                            "mean-max-min-sd-relab-nmr_no_groups-Phylum-nonbact.tsv")

mean_sd_relab.nmr_no_groups.families.nonbact.fname <-file.path(community.composition.tables,
                                                               "mean-max-min-sd-relab-nmr_no_groups-Family-nonbact.tsv")

mean_sd_relab.nmr_no_groups.genera.nonbact.fname <- file.path(community.composition.tables,
                                                              "mean-max-min-sd-relab-nmr_no_groups-Genus-nonbact.tsv")

mean_sd_relab.nmr_no_groups.asv.nonbact.fname <- file.path(community.composition.tables,
                                                           "mean-max-min-sd-relab-nmr_no_groups-OTU-nonbact.tsv")

if(!file.exists(mean_sd_relab.nmr_no_groups.phyla.nonbact.fname)){
  write.table(mean_sd_relab.nmr_no_groups.phyla.nonbact,
              file = mean_sd_relab.nmr_no_groups.phyla.nonbact.fname,
              row.names = F,sep = "\t")
}
if(!file.exists(mean_sd_relab.nmr_no_groups.families.nonbact.fname)){
  write.table(mean_sd_relab.nmr_no_groups.families.nonbact,
              file = mean_sd_relab.nmr_no_groups.families.nonbact.fname,
              row.names = F,sep = "\t")
}
if(!file.exists(mean_sd_relab.nmr_no_groups.genera.nonbact.fname)){
  write.table(mean_sd_relab.nmr_no_groups.genera.nonbact,
              file = mean_sd_relab.nmr_no_groups.genera.nonbact.fname,
              row.names = F,sep = "\t")
}
if(!file.exists(mean_sd_relab.nmr_no_groups.asv.nonbact.fname)){
  write.table(mean_sd_relab.nmr_no_groups.asv.nonbact,
              file = mean_sd_relab.nmr_no_groups.asv.nonbact.fname,
              row.names = F,sep = "\t")
}
```



## Check how much Bacteroidaceae are in NMR.


``` r
bacteroidaceae.nmr<-mean_sd_relab.nmr_no_groups.families%>%
  filter(Family=="Bacteroidaceae")
knitr::kable(bacteroidaceae.nmr, format ="simple")
```



class   Kingdom    Phylum         Family            MeanRelativeAbundance   sdRelativeAbundance   minRelativeAbundance   maxRelativeAbundance   n_samples
------  ---------  -------------  ---------------  ----------------------  --------------------  ---------------------  ---------------------  ----------
NMR     Bacteria   Bacteroidota   Bacteroidaceae                 16.01162              8.441047               5.060626               30.10842          11

``` r
bacteroidaceae.nmr.fname <-file.path(community.composition.tables,
                                  paste("kraken2", 
                                        "bacteroidaceae-nmr-table.tsv",sep="-"))
if(! file.exists (bacteroidaceae.nmr.fname)){
  write.table(bacteroidaceae.nmr,
              file = bacteroidaceae.nmr.fname,
              row.names = F,sep = "\t")
}
```



### Check the most dominant Bacteroidota families in NMR.


``` r
bacteroidota.nmr<-mean_sd_relab.nmr_no_groups.families%>%
  filter(Phylum=="Bacteroidota")
knitr::kable(bacteroidota.nmr, format ="simple")
```



class   Kingdom    Phylum         Family                  MeanRelativeAbundance   sdRelativeAbundance   minRelativeAbundance   maxRelativeAbundance   n_samples
------  ---------  -------------  ---------------------  ----------------------  --------------------  ---------------------  ---------------------  ----------
NMR     Bacteria   Bacteroidota   Bacteroidaceae                     16.0116159             8.4410466              5.0606261             30.1084186          11
NMR     Bacteria   Bacteroidota   Tannerellaceae                      5.3688697             3.1164309              0.9054130             11.0256075          11
NMR     Bacteria   Bacteroidota   Rikenellaceae                       0.7586335             0.4449989              0.2592412              1.4949202          11
NMR     Bacteria   Bacteroidota   Prevotellaceae                      0.5325308             0.2244371              0.1940218              0.8997264          11
NMR     Bacteria   Bacteroidota   Muribaculaceae                      0.1472213             0.0974074              0.0551689              0.3881501          11
NMR     Bacteria   Bacteroidota   Odoribacteraceae                    0.0190002             0.0420774              0.0000000              0.1446391          10
NMR     Bacteria   Bacteroidota   Barnesiellaceae                     0.0062140             0.0026035              0.0020956              0.0100583          11
NMR     Bacteria   Bacteroidota   Weeksellaceae                       0.0048250             0.0068263              0.0000000              0.0203422           5
NMR     Bacteria   Bacteroidota   Flavobacteriaceae                   0.0046063             0.0076553              0.0000000              0.0215975           4
NMR     Bacteria   Bacteroidota   Reichenbachiellaceae                0.0008696             0.0003844              0.0000000              0.0015074          10
NMR     Bacteria   Bacteroidota   Porphyromonadaceae                  0.0007361             0.0002619              0.0003143              0.0010198          11
NMR     Bacteria   Bacteroidota   Paludibacteraceae                   0.0004715             0.0006139              0.0000000              0.0016192           6
NMR     Bacteria   Bacteroidota   Hymenobacteraceae                   0.0003412             0.0006157              0.0000000              0.0015611           3
NMR     Bacteria   Bacteroidota   Bacteroidales Order                 0.0002826             0.0006671              0.0000000              0.0022105           3
NMR     Bacteria   Bacteroidota   Chitinophagaceae                    0.0001516             0.0005029              0.0000000              0.0016678           1
NMR     Bacteria   Bacteroidota   Marinilabiliaceae                   0.0000574             0.0001312              0.0000000              0.0003838           2
NMR     Bacteria   Bacteroidota   Sphingobacteriaceae                 0.0000391             0.0001087              0.0000000              0.0003609           2
NMR     Bacteria   Bacteroidota   Prolixibacteraceae                  0.0000287             0.0000952              0.0000000              0.0003158           1

``` r
bacteroidota.nmr.fname <-file.path(community.composition.tables,
                                     paste("kraken2", 
                                           "bacteroidota-nmr-table.tsv",sep="-"))
if(! file.exists (bacteroidota.nmr.fname)){
  write.table(bacteroidota.nmr,
              file = bacteroidota.nmr.fname,
              row.names = F,sep = "\t")
}
```



## Check Spirochaetaceae, Spirochaetota, and Treponema in NMR. 


``` r
spirochaetaceae.nmr<-mean_sd_relab.nmr_no_groups.families%>%
  filter(Family=="Spirochaetaceae")
knitr::kable(spirochaetaceae.nmr, format ="simple")
```



class   Kingdom    Phylum          Family             MeanRelativeAbundance   sdRelativeAbundance   minRelativeAbundance   maxRelativeAbundance   n_samples
------  ---------  --------------  ----------------  ----------------------  --------------------  ---------------------  ---------------------  ----------
NMR     Bacteria   Spirochaetota   Spirochaetaceae                 2.92e-05              6.71e-05                      0              0.0001978           2

``` r
spirochaetota.nmr<-mean_sd_relab.nmr_no_groups.phyla%>%
  filter(Phylum=="Spirochaetota")
knitr::kable(spirochaetota.nmr, format ="simple")
```



class   Kingdom    Phylum           MeanRelativeAbundance   sdRelativeAbundance   minRelativeAbundance   maxRelativeAbundance   n_samples
------  ---------  --------------  ----------------------  --------------------  ---------------------  ---------------------  ----------
NMR     Bacteria   Spirochaetota                0.0224319             0.0269509                      0              0.0841713          10

``` r
spirochaetaceae.nmr.fname <- file.path(community.composition.tables,
                                       paste("kraken2", 
                                             "spirochaetaceae-nmr-table.tsv",sep="-") )
spirochaetota.nmr.fname <- file.path(community.composition.tables,
                                       paste("kraken2", 
                                             "spirochaetota-nmr-table.tsv",sep="-") )

if(!file.exists(spirochaetaceae.nmr.fname)){
  write.table(spirochaetaceae.nmr,
              file = spirochaetaceae.nmr.fname,
              row.names = F,sep = "\t")
}
if(!file.exists(spirochaetota.nmr.fname)){
  write.table(spirochaetota.nmr,
              file = spirochaetota.nmr.fname,
              row.names = F,sep = "\t")
}
```

Check the number of Species in Treponema from NMR.


``` r
mean_sd_relab.nmr_no_groups.asv%>%
  ungroup()%>%
  filter(Genus=="Treponema",class=="NMR")%>%
  distinct(OTU)%>%
  tally%>%
  knitr::kable(format = "simple")
```



   n
 ---
  10

``` r
treponema.nmr<-mean_sd_relab.nmr_no_groups.genera%>%
  filter(Genus=="Treponema")
knitr::kable(treponema.nmr, format ="simple")
```



class   Kingdom    Phylum          Family            Genus        MeanRelativeAbundance   sdRelativeAbundance   minRelativeAbundance   maxRelativeAbundance   n_samples
------  ---------  --------------  ----------------  ----------  ----------------------  --------------------  ---------------------  ---------------------  ----------
NMR     Bacteria   Spirochaetota   Treponemataceae   Treponema                 0.021809             0.0268805                      0              0.0841713           9

``` r
treponema.nmr.names<-mean_sd_relab.nmr_no_groups.asv%>%
  filter(Genus=="Treponema")
knitr::kable(treponema.nmr.names, format ="simple")
```



class   Kingdom    Phylum          Family            Genus       Species                    OTU                                                                                                                           MeanRelativeAbundance   sdRelativeAbundance   minRelativeAbundance   maxRelativeAbundance   n_samples
------  ---------  --------------  ----------------  ----------  -------------------------  ---------------------------------------------------------------------------------------------------------------------------  ----------------------  --------------------  ---------------------  ---------------------  ----------
NMR     Bacteria   Spirochaetota   Treponemataceae   Treponema   Treponema_denticola        k__Bacteria|p__Spirochaetota|c__Spirochaetia|o__Spirochaetales|f__Treponemataceae|g__Treponema|s__Treponema_denticola                     0.0122255             0.0152929                      0              0.0383772           6
NMR     Bacteria   Spirochaetota   Treponemataceae   Treponema   Treponema_succinifaciens   k__Bacteria|p__Spirochaetota|c__Spirochaetia|o__Spirochaetales|f__Treponemataceae|g__Treponema|s__Treponema_succinifaciens                0.0035980             0.0119331                      0              0.0395775           1
NMR     Bacteria   Spirochaetota   Treponemataceae   Treponema   Treponema_socranskii       k__Bacteria|p__Spirochaetota|c__Spirochaetia|o__Spirochaetales|f__Treponemataceae|g__Treponema|s__Treponema_socranskii                    0.0011864             0.0022392                      0              0.0059638           3
NMR     Bacteria   Spirochaetota   Treponemataceae   Treponema   Treponema_brennaborense    k__Bacteria|p__Spirochaetota|c__Spirochaetia|o__Spirochaetales|f__Treponemataceae|g__Treponema|s__Treponema_brennaborense                 0.0010383             0.0027974                      0              0.0092385           2
NMR     Bacteria   Spirochaetota   Treponemataceae   Treponema   Treponema_pedis            k__Bacteria|p__Spirochaetota|c__Spirochaetia|o__Spirochaetales|f__Treponemataceae|g__Treponema|s__Treponema_pedis                         0.0008930             0.0020462                      0              0.0060069           2
NMR     Bacteria   Spirochaetota   Treponemataceae   Treponema   Treponema_ruminis          k__Bacteria|p__Spirochaetota|c__Spirochaetia|o__Spirochaetales|f__Treponemataceae|g__Treponema|s__Treponema_ruminis                       0.0008220             0.0014180                      0              0.0049113           7
NMR     Bacteria   Spirochaetota   Treponemataceae   Treponema   Treponema_peruense         k__Bacteria|p__Spirochaetota|c__Spirochaetia|o__Spirochaetales|f__Treponemataceae|g__Treponema|s__Treponema_peruense                      0.0007219             0.0011220                      0              0.0031922           4
NMR     Bacteria   Spirochaetota   Treponemataceae   Treponema   Treponema_putidum          k__Bacteria|p__Spirochaetota|c__Spirochaetia|o__Spirochaetales|f__Treponemataceae|g__Treponema|s__Treponema_putidum                       0.0006442             0.0021367                      0              0.0070867           1
NMR     Bacteria   Spirochaetota   Treponemataceae   Treponema   Treponema_parvum           k__Bacteria|p__Spirochaetota|c__Spirochaetia|o__Spirochaetales|f__Treponemataceae|g__Treponema|s__Treponema_parvum                        0.0005244             0.0010686                      0              0.0033756           3
NMR     Bacteria   Spirochaetota   Treponemataceae   Treponema   Treponema_vincentii        k__Bacteria|p__Spirochaetota|c__Spirochaetia|o__Spirochaetales|f__Treponemataceae|g__Treponema|s__Treponema_vincentii                     0.0001554             0.0005154                      0              0.0017094           1

``` r
treponema.nmr.fname <- file.path(community.composition.tables,
                                 paste("kraken2", 
                                       "treponema-nmr-table.tsv",sep="-") )

treponema.nmr.names.fname <- file.path(community.composition.tables,
                                       paste("kraken2", 
                                             "treponema-nmr-names.tsv",sep="-") )
if(!file.exists(treponema.nmr.fname)){
  write.table(treponema.nmr,
              file = treponema.nmr.fname,
              row.names = F,sep = "\t")
}
if(!file.exists(treponema.nmr.names.fname)){
  write.table(treponema.nmr.names,
              file = treponema.nmr.names.fname,
              row.names = F,sep = "\t")
}
```



## Analyse sulfur-metabolising bacteria in NMR. 
We don't know where the "sulf" pattern is (Phylum, Order, Genus, or anywhere 
else), so we use sapply with grepl.
https://stackoverflow.com/questions/47941680/grepl-across-multiple-specified-columns


``` r
all.ranks<-c("Kingdom", "Phylum", "Class", "Order", "Family","Genus","Species")
sulfur.bact.nmr<-
  ps.q.agg.species[!!rowSums(sapply(ps.q.agg.species[,all.ranks], grepl, 
                            pattern = "sulf|thio") ),]%>%
  distinct(Species,.keep_all = T)%>%
  select(all_of(all.ranks))%>%
  inner_join(mean_sd_relab.nmr_no_groups.asv)
```

```
## Joining with `by = join_by(Kingdom, Phylum, Family, Genus, Species)`
```

``` r
head(sulfur.bact.nmr)%>%
  knitr::kable(format = "simple")
```



Kingdom    Phylum                    Class                 Order                Family                    Genus              Species                        class   OTU                                                                                                                                                MeanRelativeAbundance   sdRelativeAbundance   minRelativeAbundance   maxRelativeAbundance   n_samples
---------  ------------------------  --------------------  -------------------  ------------------------  -----------------  -----------------------------  ------  ------------------------------------------------------------------------------------------------------------------------------------------------  ----------------------  --------------------  ---------------------  ---------------------  ----------
Bacteria   Pseudomonadota            Gammaproteobacteria   Chromatiales         Ectothiorhodospiraceae    Thioalkalivibrio   Thioalkalivibrio_versutus      NMR     k__Bacteria|p__Pseudomonadota|c__Gammaproteobacteria|o__Chromatiales|f__Ectothiorhodospiraceae|g__Thioalkalivibrio|s__Thioalkalivibrio_versutus                0.0086914             0.0152498                      0              0.0374657           3
Bacteria   Thermodesulfobacteriota   Desulfovibrionia      Desulfovibrionales   Desulfovibrionaceae       Desulfovibrio      Desulfovibrio_piger            NMR     k__Bacteria|p__Thermodesulfobacteriota|c__Desulfovibrionia|o__Desulfovibrionales|f__Desulfovibrionaceae|g__Desulfovibrio|s__Desulfovibrio_piger                0.0030931             0.0032614                      0              0.0097499           7
Bacteria   Thermodesulfobacteriota   Desulfovibrionia      Desulfovibrionales   Desulfovibrionaceae       Taurinivorans      Taurinivorans_muris            NMR     k__Bacteria|p__Thermodesulfobacteriota|c__Desulfovibrionia|o__Desulfovibrionales|f__Desulfovibrionaceae|g__Taurinivorans|s__Taurinivorans_muris                0.0008042             0.0022070                      0              0.0073123           2
Bacteria   Synergistota              Synergistia           Synergistales        Dethiosulfovibrionaceae   Jonquetella        Jonquetella_anthropi           NMR     k__Bacteria|p__Synergistota|c__Synergistia|o__Synergistales|f__Dethiosulfovibrionaceae|g__Jonquetella|s__Jonquetella_anthropi                                  0.0006858             0.0008508                      0              0.0024127           5
Bacteria   Thermodesulfobacteriota   Desulfobulbia         Desulfobulbales      Desulfobulbaceae          Desulfobulbus      Desulfobulbus_oligotrophicus   NMR     k__Bacteria|p__Thermodesulfobacteriota|c__Desulfobulbia|o__Desulfobulbales|f__Desulfobulbaceae|g__Desulfobulbus|s__Desulfobulbus_oligotrophicus                0.0001801             0.0005973                      0              0.0019809           1
Archaea    Thermoproteota            Thermoprotei          Desulfurococcales    Desulfurococcaceae        Aeropyrum          Aeropyrum_camini               NMR     k__Archaea|p__Thermoproteota|c__Thermoprotei|o__Desulfurococcales|f__Desulfurococcaceae|g__Aeropyrum|s__Aeropyrum_camini                                       0.0001182             0.0003693                      0              0.0012309           3

``` r
sulfur.bact.nmr.fname <- file.path(community.composition.tables,
                                   paste("kraken2", 
                                         "sulfur.bact.nmr-nmr-table.tsv",sep="-") )
if(!file.exists(sulfur.bact.nmr.fname)){
  write.table(sulfur.bact.nmr,
              file = sulfur.bact.nmr.fname,
              row.names = F,sep = "\t")
}
```



### Plot sulfur metabolising bacteria.


``` r
sulfur.bact.plot <- ps.q.rel.agg.asv%>%
  filter(Species %in% sulfur.bact.nmr$Species)%>%
  mutate(Species = factor(Species, levels = sulfur.bact.nmr$Species),
         Species = gsub("_"," ",Species),
         Species = paste0("<i>",Species,"</i>")
  )%>%
  left_join(custom.md[c("Sample","agegroup")])%>%
  mutate(Sample=factor(Sample,levels=sample.levels$Sample))%>%
  group_by_at(c("class","Species"))%>%
  ggplot(aes(x=Sample,
             y=RelativeAbundance,
             fill=factor(agegroup)))+
  geom_bar(stat="identity")+
  facet_wrap(~Species,
             scales = "free",
             ncol = 3)+
  theme_bw()+
  labs(x="",
       y="Relative abundance (%)",
       fill=gg.labs.name)+
  coord_cartesian(expand = c("bottom" = F))+
  scale_color_manual(breaks = unname(pretty.level.names),
                     labels=unname(pretty.level.names))+
  scale_x_discrete(labels=sample.levels$Sample,
                   limits=sample.levels$Sample)+
  # scale_fill_manual(labels=pretty.level.names)+
  scale_fill_viridis_d(option = "C")+
  theme(axis.title.y = element_text(size = 15),
        axis.title = element_text(size = 10),
        axis.text.y = ggtext::element_markdown(size=10),
        axis.text.x = element_text(size=6),
        strip.text.x = ggtext::element_markdown(size=10),
        plot.title = element_text(size = 17),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 15),
        legend.position = "right",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
```

```
## Joining with `by = join_by(Sample)`
```

``` r
print(sulfur.bact.plot +
        theme(axis.title.y = element_text(size = 5),
              axis.title = element_text(size = 5),
              axis.text.y = ggtext::element_markdown(size=5),
              axis.text.x = element_text(size=6, angle=45, hjust=1),
              strip.text.x = ggtext::element_markdown(size=6),
              plot.title = element_text(size = 10),
              legend.text = element_text(size = 8),
              legend.title = element_text(size = 10),
              legend.position = "right",
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank())+
        ggtitle(paste0("Relative abundance of sulfur-utilizing bacteria in different naked mole-rat age groups")))
```

<img src="003-summary-stats-kraken2_files/figure-html/unnamed-chunk-49-1.png" alt="" width="672" />

``` r
for (image.format in c("png","tiff")){
  ggsave(filename = paste0("sulfur-bacteria-nmr.", image.format),
         path = community.composition.figures,
         plot=sulfur.bact.plot,
         width=11, height=8,units="in",
         dpi=300,device = image.format)
}

treponema.plot<-ps.q.rel.agg.asv%>%
  filter(Genus=="Treponema")%>%
  arrange(Species)%>%
  mutate(Species = factor(Species, levels = unique(Species)),
         Species = gsub("_"," ",Species),
         Species = paste0("<i>",Species,"</i>")
  )%>%
  left_join(custom.md[c("Sample","agegroup")])%>%
  mutate(Sample=factor(Sample,levels=sample.levels$Sample))%>%
  group_by_at(c("class","Species"))%>%
  ggplot(aes(x=Sample,
             y=RelativeAbundance,
             fill=factor(agegroup)))+
  geom_bar(stat="identity")+
  facet_wrap(~Species,
             scales = "free",
             ncol = 3)+
  theme_bw()+
  labs(x="",
       y="Relative abundance (%)",
       fill=gg.labs.name)+
  coord_cartesian(expand = c("bottom" = F))+
  scale_color_manual(breaks = unname(pretty.level.names),
                     labels=unname(pretty.level.names))+
  scale_x_discrete(labels=sample.levels$Sample,
                   limits=sample.levels$Sample)+
  scale_fill_viridis_d(option = "C")+
  theme(axis.title.y = element_text(size = 15),
        axis.title = element_text(size = 10),
        axis.text.y = ggtext::element_markdown(size=10),
        axis.text.x = element_text(size=6),
        strip.text.x = ggtext::element_markdown(size=10),
        plot.title = element_text(size = 17),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 15),
        legend.position = "right",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())
```

```
## Joining with `by = join_by(Sample)`
```

``` r
print(treponema.plot+
        theme(axis.title.y = element_text(size = 5),
              axis.title = element_text(size = 5),
              axis.text.y = ggtext::element_markdown(size=5),
              axis.text.x = element_text(size=6, angle=45, hjust=1),
              strip.text.x = ggtext::element_markdown(size=8),
              plot.title = element_text(size = 10),
              legend.text = element_text(size = 8),
              legend.title = element_text(size = 10),
              legend.position = "right",
              panel.grid.minor = element_blank(),
              panel.grid.major = element_blank())+
        ggtitle(paste0("Relative abundance of Treponemas in naked mole-rats")))
```

<img src="003-summary-stats-kraken2_files/figure-html/unnamed-chunk-50-1.png" alt="" width="672" />

``` r
for (image.format in c("png","tiff")){
  ggsave(filename = paste0("treponema.sp-nmr.", image.format),
         path = community.composition.figures,
         plot=treponema.plot,
         width=11, height=8,units="in",
         dpi=300,device = image.format)
}
```



## Groups 1, 2, and 3 (**Supplementary figure**).


``` r
group1.genera<-c("Faecalibacterium",
                 "Roseburia",
                 "Coprococcus",
                 "Prevotella")
group1.species<-c("Eubacterium_rectale") # Not in our data
group2.increased.genera<-c("Eggerthella",
                           "Bilophila",
                           "Desulfovibrio",
                           "Fusobacterium",
                           "Anaerotruncus",
                           "Streptococcus",
                           "Escherichia")
group2.unhealthy.genera<-c("Eggerthella",
                           "Coprobacillus",
                           "Streptococcus",
                           "Bilophila",
                           "Actinomyces",
                           "Desulfovibrio",
                           "Campylobacter",
                           "Veillonella",
                           "Enterococcus")
group2.unhealthy.families<-c("Atopobiaceae",
                             "Enterobacteriaceae")
group2.unhealthy.species<-c("Bacteroides_fragilis",
                            "Clostridium_hathewayi",
                            "Enterocloster_clostridioformis", #Clostridium_clostridioforme
                            "Enterocloster_bolteae", # Clostridium_bolteae
                            "Clostridium_cindens",
                            "[Ruminococcus]_torques", #Ruminococcus_torques
                            "Mediterraneibacter_gnavus", # Ruminococcus_gnavus
                            "Clostridioides_difficile")

group3.genera<-c("Akkermansia",
                 "Odoribacter",
                 "Butyricimonas",
                 "Butyrivibrio",
                 "Barnesiella",
                 "Oscillospira")
group3.families<-c("Christensenellaceae")
```

Abundances of group 1 genera (decreased with age and associated with 
healthy aging in humans):


``` r
group1.genera
```

```
## [1] "Faecalibacterium" "Roseburia"        "Coprococcus"      "Prevotella"
```

``` r
group1.genera.plot<-ggplot_species(group1.genera,
                                   ps.q.rel.agg.asv,
                                   tax.rank="Genus",
                                   sample.order=sample.levels,
                                   group.names = pretty.level.names,
                                   grouping.variable = "agegroup",
                                   metadata.df = custom.md,
                                   ggplot.fill.name = gg.labs.name)
```

```
## Joining with `by = join_by(Sample)`
```

``` r
print(group1.genera.plot+
        theme(axis.text.x = element_text(size=6, angle=45, hjust=1))+
        ggtitle(paste0("Relative abundance of Group 1 members"))) # 4 plots
```

<img src="003-summary-stats-kraken2_files/figure-html/unnamed-chunk-54-1.png" alt="" width="672" />

``` r
group1.species%in%ps.q.rel.agg.asv$Species # no need to plot
```

```
## [1] FALSE
```

Abundances of group 2 genera (increased with age):


``` r
group2.increased.genera
```

```
## [1] "Eggerthella"   "Bilophila"     "Desulfovibrio" "Fusobacterium"
## [5] "Anaerotruncus" "Streptococcus" "Escherichia"
```

``` r
group2.increased.genera.plot<-ggplot_species(taxa.to.plot =group2.increased.genera,
                                             tax.df =ps.q.rel.agg.asv,
                                             tax.rank="Genus",
                                             sample.order=sample.levels,
                                             group.names = pretty.level.names,
                                             grouping.variable = "agegroup",
                                             metadata.df = custom.md,
                                             ggplot.fill.name = gg.labs.name)
```

```
## Joining with `by = join_by(Sample)`
```

``` r
print(group2.increased.genera.plot+
        theme(plot.title = element_text(size = 13),
              axis.text.x = element_text(size=6, angle=45, hjust=1))+
        ggtitle(paste0("Relative abundance of Group 2 members increased with age")))
```

<img src="003-summary-stats-kraken2_files/figure-html/unnamed-chunk-56-1.png" alt="" width="672" />

Abundances of group 2 genera (associated with unhealthy aging):


``` r
group2.unhealthy.genera
```

```
## [1] "Eggerthella"   "Coprobacillus" "Streptococcus" "Bilophila"    
## [5] "Actinomyces"   "Desulfovibrio" "Campylobacter" "Veillonella"  
## [9] "Enterococcus"
```

``` r
group2.unhealthy.genera.plot<-ggplot_species(taxa.to.plot=group2.unhealthy.genera,
                                              tax.df=ps.q.rel.agg.asv,
                                              tax.rank="Genus",
                                              sample.order=sample.levels,
                                              group.names = pretty.level.names,
                                              grouping.variable = "agegroup",
                                              metadata.df = custom.md,
                                              ggplot.fill.name = gg.labs.name)
```

```
## Joining with `by = join_by(Sample)`
```

``` r
print(group2.unhealthy.genera.plot+
        theme(plot.title = element_text(size = 13),
              axis.text.x = element_text(size=6, angle=45, hjust=1))+
  ggtitle(paste0("Relative abundance of Group 2 members associated with unhealthy aging")))
```

<img src="003-summary-stats-kraken2_files/figure-html/unnamed-chunk-58-1.png" alt="" width="672" />

Abundances of group 2 families (associated with unhealthy aging):


``` r
group2.unhealthy.families
```

```
## [1] "Atopobiaceae"       "Enterobacteriaceae"
```

``` r
group2.unhealthy.families.plot<-ggplot_species(taxa.to.plot=group2.unhealthy.families,
                                              tax.df=ps.q.rel.agg.asv,
                                              tax.rank="Family",
                                              sample.order=sample.levels,
                                              group.names = pretty.level.names,
                                              grouping.variable = "agegroup",
                                              metadata.df = custom.md,
                                              ggplot.fill.name = gg.labs.name)
```

```
## Joining with `by = join_by(Sample)`
```

``` r
print(group2.unhealthy.families.plot+
        theme(plot.title = element_text(size = 13),
              axis.text.x = element_text(size=6, angle=45, hjust=1))+
        ggtitle(paste0("Relative abundance of Group 2 members associated with unhealthy aging")))
```

<img src="003-summary-stats-kraken2_files/figure-html/unnamed-chunk-60-1.png" alt="" width="672" />

Abundances of group 2 species (associated with unhealthy aging):


``` r
group2.unhealthy.species
```

```
## [1] "Bacteroides_fragilis"           "Clostridium_hathewayi"         
## [3] "Enterocloster_clostridioformis" "Enterocloster_bolteae"         
## [5] "Clostridium_cindens"            "[Ruminococcus]_torques"        
## [7] "Mediterraneibacter_gnavus"      "Clostridioides_difficile"
```

``` r
group2.unhealthy.species.plot<-ggplot_species(taxa.to.plot =group2.unhealthy.species,
                                             tax.df = ps.q.rel.agg.asv,
                                              tax.rank="Species",
                                             sample.order=sample.levels,
                                             group.names = pretty.level.names,
                                             grouping.variable = "agegroup",
                                             metadata.df = custom.md,
                                             ggplot.fill.name = gg.labs.name)
```

```
## Joining with `by = join_by(Sample)`
```

``` r
print(group2.unhealthy.species.plot+
        theme(plot.title = element_text(size = 13),
              axis.text.x = element_text(size=6, angle=45, hjust=1))+
        ggtitle(paste0("Relative abundance of Group 2 members associated with unhealthy aging")))
```

<img src="003-summary-stats-kraken2_files/figure-html/unnamed-chunk-62-1.png" alt="" width="672" />

Abundances of group 3 genera (increased with age and associated with 
healthy aging):


``` r
group3.genera
```

```
## [1] "Akkermansia"   "Odoribacter"   "Butyricimonas" "Butyrivibrio" 
## [5] "Barnesiella"   "Oscillospira"
```

``` r
group3.genera.plot<-ggplot_species(taxa.to.plot =group3.genera,
                                    tax.df = ps.q.rel.agg.asv,
                                    tax.rank = "Genus",
                                    sample.order=sample.levels,
                                    group.names = pretty.level.names,
                                    grouping.variable = "agegroup",
                                    metadata.df = custom.md,
                                    ggplot.fill.name = gg.labs.name)
```

```
## Joining with `by = join_by(Sample)`
```

``` r
print(group3.genera.plot+
        theme(axis.text.x = element_text(size=6, angle=45, hjust=1))+
        ggtitle(paste0("Relative abundance of Group 3 members")))
```

<img src="003-summary-stats-kraken2_files/figure-html/unnamed-chunk-64-1.png" alt="" width="672" />

Abundances of group 3 families (increased with age and associated 
with healthy aging):


``` r
group3.families.plot<-ggplot_species(taxa.to.plot = group3.families,
                                      tax.df = ps.q.rel.agg.asv,
                                    tax.rank = "Family",
                                    sample.order=sample.levels,
                                    group.names = pretty.level.names,
                                    grouping.variable = "agegroup",
                                    metadata.df = custom.md,
                                    ggplot.fill.name = gg.labs.name)
```

```
## Joining with `by = join_by(Sample)`
```

``` r
print(group3.families.plot+
        theme(axis.text.x = element_text(size=6, angle=45, hjust=1))+
        ggtitle(paste0("Relative abundance of Group 3 members")))
```

<img src="003-summary-stats-kraken2_files/figure-html/unnamed-chunk-66-1.png" alt="" width="672" />

``` r
for (image.format in c("png","tiff")){
  ggsave(filename = paste0("group1.genera-nmr.", image.format),
         path = community.composition.figures,
         plot=group1.genera.plot,
         width=10, height=7,units="in",
         dpi=300,device = image.format) # 4 plots
  ggsave(filename = paste0("group2.increased.genera-nmr.",image.format),
         path = community.composition.figures,
         plot=group2.increased.genera.plot,
         width=10, height=7,units="in",
         dpi=300,device = image.format) # 4 plots
  ggsave(filename = paste0("group2.unhealthy.genera-nmr.",image.format),
         path = community.composition.figures,
         plot=group2.unhealthy.genera.plot,
         width=10, height=7,units="in",
         dpi=300,device = image.format) # 7 plots
  ggsave(filename = paste0("group2.unhealthy.families-nmr.",image.format),
         path = community.composition.figures,
         plot = group2.unhealthy.families.plot,
         width=10, height=5,units="in",
         dpi=300,device = image.format)    # 2 plots
  ggsave(filename = paste0("group2.unhealthy.species-nmr.",image.format),
         path = community.composition.figures,
         plot = group2.unhealthy.species.plot,
         width=10, height=8,units="in",
         dpi=300,device = image.format)    # 6 plots
  ggsave(paste0("group3.genera-nmr.", image.format), 
         path = community.composition.figures,
         plot= group3.genera.plot ,
         width=10, height=8,units="in",
         dpi=300,device = image.format)    # 5 plots
  ggsave(paste0("group3.families-nmr.", image.format),
         path = community.composition.figures,
         plot = group3.families.plot ,
         width=7, height=5,units="in",
         dpi=300,device = image.format)    # 1 plot

}
```



## Analyse each kingdom in detail. 

Number of bacterial species


``` r
mean_sd_relab.nmr_no_groups.asv%>%
  filter(Kingdom=="Bacteria")%>%
  nrow
```

```
## [1] 492
```



### Eukaryotes
How many eukaryote species do we have


``` r
mean_sd_relab.nmr_no_groups.asv.nonbact%>%
  filter(Kingdom=="Eukaryota")%>%
  distinct(Species)%>%
  tally()
```

```
## # A tibble: 1 × 1
##       n
##   <int>
## 1   100
```

Total abundance of eukaryotes in classified data


``` r
mean_sd_relab.nmr_no_groups.asv.nonbact%>%
  ungroup%>%
  filter(Kingdom=="Eukaryota")%>%
  distinct() %>%
  knitr::kable(format = "simple")
```



class   Kingdom     Phylum          Family                     Genus                 Species                              OTU                                                                                                                                                               MeanRelativeAbundance   sdRelativeAbundance   minRelativeAbundance   maxRelativeAbundance   n_samples
------  ----------  --------------  -------------------------  --------------------  -----------------------------------  ---------------------------------------------------------------------------------------------------------------------------------------------------------------  ----------------------  --------------------  ---------------------  ---------------------  ----------
NMR     Eukaryota   Streptophyta    Convolvulaceae             Ipomoea               Ipomoea_triloba                      k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Solanales|f__Convolvulaceae|g__Ipomoea|s__Ipomoea_triloba                                                                    4.3942845             2.6792958              1.7091821             10.5267146          11
NMR     Eukaryota   Streptophyta    Cupressaceae               Cryptomeria           Cryptomeria_japonica                 k__Eukaryota|p__Streptophyta|c__Pinopsida|o__Cupressales|f__Cupressaceae|g__Cryptomeria|s__Cryptomeria_japonica                                                               2.6092737             0.3424485              2.3218895              3.3903297          11
NMR     Eukaryota   Streptophyta    Proteaceae                 Macadamia             Macadamia_integrifolia               k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Proteales|f__Proteaceae|g__Macadamia|s__Macadamia_integrifolia                                                               2.0678735             0.9723804              0.8940231              3.6209573          11
NMR     Eukaryota   Streptophyta    Fabaceae                   Glycine               Glycine_soja                         k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Fabales|f__Fabaceae|g__Glycine|s__Glycine_soja                                                                               1.2108200             2.1318113              0.0000000              5.7126046           5
NMR     Eukaryota   Streptophyta    Salicaceae                 Populus               Populus_trichocarpa                  k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Malpighiales|f__Salicaceae|g__Populus|s__Populus_trichocarpa                                                                 0.9575214             0.7136691              0.3688666              2.4784448          11
NMR     Eukaryota   Streptophyta    Fabaceae                   Arachis               Arachis_hypogaea                     k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Fabales|f__Fabaceae|g__Arachis|s__Arachis_hypogaea                                                                           0.7151611             1.0374879              0.0000000              3.3176101           5
NMR     Eukaryota   Streptophyta    Salicaceae                 Populus               Populus_nigra                        k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Malpighiales|f__Salicaceae|g__Populus|s__Populus_nigra                                                                       0.6855855             0.7582252              0.0000000              1.9464763           6
NMR     Eukaryota   Parabasalia     Trichomonadidae            Trichomonas           Trichomonas_vaginalis                k__Eukaryota|p__Parabasalia|NA|o__Trichomonadida|f__Trichomonadidae|g__Trichomonas|s__Trichomonas_vaginalis                                                                   0.6564723             0.5321252              0.1328118              1.7745074          11
NMR     Eukaryota   Streptophyta    Cucurbitaceae              Cucurbita             Cucurbita_pepo                       k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Cucurbitales|f__Cucurbitaceae|g__Cucurbita|s__Cucurbita_pepo                                                                 0.6400268             0.6072243              0.0000000              1.9084192           9
NMR     Eukaryota   Streptophyta    Solanaceae                 Capsicum              Capsicum_annuum                      k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Solanales|f__Solanaceae|g__Capsicum|s__Capsicum_annuum                                                                       0.6317124             0.3555059              0.2183776              1.2555887          11
NMR     Eukaryota   Streptophyta    Fabaceae                   Medicago              Medicago_truncatula                  k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Fabales|f__Fabaceae|g__Medicago|s__Medicago_truncatula                                                                       0.5781865             0.3389361              0.1980142              1.2346710          11
NMR     Eukaryota   Streptophyta    Poaceae                    Triticum              Triticum_dicoccoides                 k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Poales|f__Poaceae|g__Triticum|s__Triticum_dicoccoides                                                                        0.5543725             0.4327733              0.0000000              1.1225333           8
NMR     Eukaryota   Ascomycota      Pichiaceae                 Pichia                Pichia_kudriavzevii                  k__Eukaryota|p__Ascomycota|c__Saccharomycetes|o__Saccharomycetales|f__Pichiaceae|g__Pichia|s__Pichia_kudriavzevii                                                             0.5206977             0.9007816              0.0000000              2.3987543           8
NMR     Eukaryota   Streptophyta    Fabaceae                   Vicia                 Vicia_villosa                        k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Fabales|f__Fabaceae|g__Vicia|s__Vicia_villosa                                                                                0.4741428             0.3028619              0.0000000              0.9873599          10
NMR     Eukaryota   Streptophyta    Ximeniaceae                Malania               Malania_oleifera                     k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Santalales|f__Ximeniaceae|g__Malania|s__Malania_oleifera                                                                     0.4353363             0.1918615              0.2008671              0.7919929          11
NMR     Eukaryota   Streptophyta    Fabaceae                   Cicer                 Cicer_arietinum                      k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Fabales|f__Fabaceae|g__Cicer|s__Cicer_arietinum                                                                              0.3916369             0.4878842              0.0000000              1.5989451           7
NMR     Eukaryota   Streptophyta    Fabaceae                   Pisum                 Pisum_sativum                        k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Fabales|f__Fabaceae|g__Pisum|s__Pisum_sativum                                                                                0.3628890             0.2299912              0.1329909              0.8427452          11
NMR     Eukaryota   Ascomycota      Aspergillaceae             Aspergillus           Aspergillus_puulaauensis             k__Eukaryota|p__Ascomycota|c__Eurotiomycetes|o__Eurotiales|f__Aspergillaceae|g__Aspergillus|s__Aspergillus_puulaauensis                                                       0.3496001             0.5269957              0.0000000              1.7242715           8
NMR     Eukaryota   Streptophyta    Asteraceae                 Lactuca               Lactuca_sativa                       k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Asterales|f__Asteraceae|g__Lactuca|s__Lactuca_sativa                                                                         0.3475837             0.2163442              0.0768528              0.7642434          11
NMR     Eukaryota   Streptophyta    Fabaceae                   Cajanus               Cajanus_cajan                        k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Fabales|f__Fabaceae|g__Cajanus|s__Cajanus_cajan                                                                              0.3137126             0.4832510              0.0000000              1.5356188           5
NMR     Eukaryota   Streptophyta    Dioscoreaceae              Dioscorea             Dioscorea_cayenensis                 k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Dioscoreales|f__Dioscoreaceae|g__Dioscorea|s__Dioscorea_cayenensis                                                           0.3083730             0.2427848              0.0000000              0.7280405           9
NMR     Eukaryota   Streptophyta    Oleaceae                   Olea                  Olea_europaea                        k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Lamiales|f__Oleaceae|g__Olea|s__Olea_europaea                                                                                0.2996996             0.2736460              0.1179943              1.0917522          11
NMR     Eukaryota   Streptophyta    Fabaceae                   Trifolium             Trifolium_pratense                   k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Fabales|f__Fabaceae|g__Trifolium|s__Trifolium_pratense                                                                       0.2853775             0.3922929              0.0000000              1.0114551           5
NMR     Eukaryota   Streptophyta    Fabaceae                   Lotus                 Lotus_japonicus                      k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Fabales|f__Fabaceae|g__Lotus|s__Lotus_japonicus                                                                              0.2657827             0.2137448              0.0000000              0.6026839           8
NMR     Eukaryota   Streptophyta    Cannabaceae                Humulus               Humulus_lupulus                      k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Rosales|f__Cannabaceae|g__Humulus|s__Humulus_lupulus                                                                         0.2544362             0.2327878              0.0000000              0.7988457           8
NMR     Eukaryota   Streptophyta    Solanaceae                 Solanum               Solanum_stenotomum                   k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Solanales|f__Solanaceae|g__Solanum|s__Solanum_stenotomum                                                                     0.2361178             0.2708437              0.0000000              0.8633769           7
NMR     Eukaryota   Streptophyta    Fagaceae                   Quercus               Quercus_robur                        k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Fagales|f__Fagaceae|g__Quercus|s__Quercus_robur                                                                              0.2341839             0.2978310              0.0000000              0.8649622           6
NMR     Eukaryota   Streptophyta    Chenopodiaceae             Spinacia              Spinacia_oleracea                    k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Caryophyllales|f__Chenopodiaceae|g__Spinacia|s__Spinacia_oleracea                                                            0.2134506             0.3845052              0.0000000              1.0842872           3
NMR     Eukaryota   Streptophyta    Ebenaceae                  Diospyros             Diospyros_lotus                      k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Ericales|f__Ebenaceae|g__Diospyros|s__Diospyros_lotus                                                                        0.1954061             0.2473246              0.0000000              0.8612026           8
NMR     Eukaryota   Streptophyta    Fabaceae                   Vigna                 Vigna_radiata                        k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Fabales|f__Fabaceae|g__Vigna|s__Vigna_radiata                                                                                0.1779137             0.2564518              0.0000000              0.8706590           6
NMR     Eukaryota   Streptophyta    Nymphaeaceae               Nymphaea              Nymphaea_colorata                    k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Nymphaeales|f__Nymphaeaceae|g__Nymphaea|s__Nymphaea_colorata                                                                 0.1656028             0.1169881              0.0000000              0.3798656           9
NMR     Eukaryota   Streptophyta    Papaveraceae               Papaver               Papaver_somniferum                   k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Ranunculales|f__Papaveraceae|g__Papaver|s__Papaver_somniferum                                                                0.1434742             0.0811289              0.0425880              0.3123466          11
NMR     Eukaryota   Streptophyta    Asteraceae                 Cynara                Cynara_cardunculus                   k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Asterales|f__Asteraceae|g__Cynara|s__Cynara_cardunculus                                                                      0.1426015             0.2742980              0.0000000              0.9419459           6
NMR     Eukaryota   Streptophyta    Zingiberaceae              Zingiber              Zingiber_officinale                  k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Zingiberales|f__Zingiberaceae|g__Zingiber|s__Zingiber_officinale                                                             0.1361383             0.0806995              0.0000000              0.2242024           9
NMR     Eukaryota   Streptophyta    Solanaceae                 Solanum               Solanum_dulcamara                    k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Solanales|f__Solanaceae|g__Solanum|s__Solanum_dulcamara                                                                      0.1239114             0.1051990              0.0000000              0.2714040           7
NMR     Eukaryota   Streptophyta    Brassicaceae               Raphanus              Raphanus_sativus                     k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Brassicales|f__Brassicaceae|g__Raphanus|s__Raphanus_sativus                                                                  0.1205142             0.2797971              0.0000000              0.8416333           2
NMR     Eukaryota   Streptophyta    Poaceae                    Hordeum               Hordeum_vulgare                      k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Poales|f__Poaceae|g__Hordeum|s__Hordeum_vulgare                                                                              0.1177818             0.0883130              0.0000000              0.2764294           9
NMR     Eukaryota   Streptophyta    Anacardiaceae              Mangifera             Mangifera_indica                     k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Sapindales|f__Anacardiaceae|g__Mangifera|s__Mangifera_indica                                                                 0.1169308             0.1426927              0.0000000              0.3665450           5
NMR     Eukaryota   Streptophyta    Poaceae                    Sorghum               Sorghum_bicolor                      k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Poales|f__Poaceae|g__Sorghum|s__Sorghum_bicolor                                                                              0.1157767             0.0866204              0.0000000              0.2567593           9
NMR     Eukaryota   Streptophyta    Euphorbiaceae              Manihot               Manihot_esculenta                    k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Malpighiales|f__Euphorbiaceae|g__Manihot|s__Manihot_esculenta                                                                0.1099878             0.2333774              0.0000000              0.7547693           3
NMR     Eukaryota   Streptophyta    Asparagaceae               Asparagus             Asparagus_officinalis                k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Asparagales|f__Asparagaceae|g__Asparagus|s__Asparagus_officinalis                                                            0.1059029             0.1407302              0.0000000              0.3848442           6
NMR     Eukaryota   Streptophyta    Betulaceae                 Alnus                 Alnus_glutinosa                      k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Fagales|f__Betulaceae|g__Alnus|s__Alnus_glutinosa                                                                            0.1057659             0.1634050              0.0000000              0.4821699           4
NMR     Eukaryota   Streptophyta    Euphorbiaceae              Hevea                 Hevea_brasiliensis                   k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Malpighiales|f__Euphorbiaceae|g__Hevea|s__Hevea_brasiliensis                                                                 0.1025399             0.1272880              0.0000000              0.3228165           5
NMR     Eukaryota   Streptophyta    Arecaceae                  Elaeis                Elaeis_guineensis                    k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Arecales|f__Arecaceae|g__Elaeis|s__Elaeis_guineensis                                                                         0.1001719             0.1229870              0.0000000              0.3964505           6
NMR     Eukaryota   Streptophyta    Fabaceae                   Lupinus               Lupinus_angustifolius                k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Fabales|f__Fabaceae|g__Lupinus|s__Lupinus_angustifolius                                                                      0.0962704             0.1455330              0.0000000              0.4966562           6
NMR     Eukaryota   Streptophyta    Euphorbiaceae              Mercurialis           Mercurialis_annua                    k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Malpighiales|f__Euphorbiaceae|g__Mercurialis|s__Mercurialis_annua                                                            0.0872185             0.2260812              0.0000000              0.7391022           2
NMR     Eukaryota   Streptophyta    Solanaceae                 Nicotiana             Nicotiana_attenuata                  k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Solanales|f__Solanaceae|g__Nicotiana|s__Nicotiana_attenuata                                                                  0.0827604             0.1110623              0.0000000              0.2847000           6
NMR     Eukaryota   Streptophyta    Vitaceae                   Vitis                 Vitis_riparia                        k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Vitales|f__Vitaceae|g__Vitis|s__Vitis_riparia                                                                                0.0823062             0.1832309              0.0000000              0.4669490           2
NMR     Eukaryota   Streptophyta    Poaceae                    Phragmites            Phragmites_australis                 k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Poales|f__Poaceae|g__Phragmites|s__Phragmites_australis                                                                      0.0791282             0.0607035              0.0000000              0.1646163           8
NMR     Eukaryota   Streptophyta    Poaceae                    Setaria               Setaria_viridis                      k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Poales|f__Poaceae|g__Setaria|s__Setaria_viridis                                                                              0.0772970             0.1197944              0.0000000              0.3963346           6
NMR     Eukaryota   Streptophyta    Apiaceae                   Daucus                Daucus_carota                        k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Apiales|f__Apiaceae|g__Daucus|s__Daucus_carota                                                                               0.0771202             0.1357514              0.0000000              0.3582752           3
NMR     Eukaryota   Streptophyta    Magnoliaceae               Magnolia              Magnolia_sinica                      k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Magnoliales|f__Magnoliaceae|g__Magnolia|s__Magnolia_sinica                                                                   0.0769083             0.0692871              0.0000000              0.2174395           8
NMR     Eukaryota   Streptophyta    Euphorbiaceae              Ricinus               Ricinus_communis                     k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Malpighiales|f__Euphorbiaceae|g__Ricinus|s__Ricinus_communis                                                                 0.0758454             0.1120358              0.0000000              0.2909059           4
NMR     Eukaryota   Streptophyta    Asteraceae                 Helianthus            Helianthus_annuus                    k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Asterales|f__Asteraceae|g__Helianthus|s__Helianthus_annuus                                                                   0.0728783             0.0650898              0.0000000              0.1926630           8
NMR     Eukaryota   Streptophyta    Malvaceae                  Gossypium             Gossypium_hirsutum                   k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Malvales|f__Malvaceae|g__Gossypium|s__Gossypium_hirsutum                                                                     0.0727589             0.1087027              0.0000000              0.2869660           4
NMR     Eukaryota   Streptophyta    Fabaceae                   Vigna                 Vigna_angularis                      k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Fabales|f__Fabaceae|g__Vigna|s__Vigna_angularis                                                                              0.0650134             0.2156251              0.0000000              0.7151475           1
NMR     Eukaryota   Streptophyta    Cucurbitaceae              Cucumis               Cucumis_melo                         k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Cucurbitales|f__Cucurbitaceae|g__Cucumis|s__Cucumis_melo                                                                     0.0581199             0.1927619              0.0000000              0.6393189           1
NMR     Eukaryota   Streptophyta    Cucurbitaceae              Benincasa             Benincasa_hispida                    k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Cucurbitales|f__Cucurbitaceae|g__Benincasa|s__Benincasa_hispida                                                              0.0549863             0.0658992              0.0000000              0.1664861           5
NMR     Eukaryota   Streptophyta    Fabaceae                   Vigna                 Vigna_unguiculata                    k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Fabales|f__Fabaceae|g__Vigna|s__Vigna_unguiculata                                                                            0.0518170             0.1174190              0.0000000              0.3348174           2
NMR     Eukaryota   Streptophyta    Amaranthaceae              Amaranthus            Amaranthus_tricolor                  k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Caryophyllales|f__Amaranthaceae|g__Amaranthus|s__Amaranthus_tricolor                                                         0.0495862             0.1119620              0.0000000              0.3154145           2
NMR     Eukaryota   Streptophyta    Poaceae                    Triticum              Triticum_aestivum                    k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Poales|f__Poaceae|g__Triticum|s__Triticum_aestivum                                                                           0.0458072             0.0498172              0.0000000              0.1249015           6
NMR     Eukaryota   Streptophyta    Poaceae                    Oryza                 Oryza_sativa                         k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Poales|f__Poaceae|g__Oryza|s__Oryza_sativa                                                                                   0.0456287             0.0715483              0.0000000              0.1990071           4
NMR     Eukaryota   Streptophyta    Juglandaceae               Juglans               Juglans_microcarpa_x_Juglans_regia   k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Fagales|f__Juglandaceae|g__Juglans|s__Juglans_microcarpa_x_Juglans_regia                                                     0.0452890             0.0682368              0.0000000              0.1967393           4
NMR     Eukaryota   Streptophyta    Proteaceae                 Telopea               Telopea_speciosissima                k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Proteales|f__Proteaceae|g__Telopea|s__Telopea_speciosissima                                                                  0.0449717             0.0302270              0.0000000              0.0797439           8
NMR     Eukaryota   Streptophyta    Ericaceae                  Rhododendron          Rhododendron_vialii                  k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Ericales|f__Ericaceae|g__Rhododendron|s__Rhododendron_vialii                                                                 0.0421433             0.0656541              0.0000000              0.1786530           4
NMR     Eukaryota   Streptophyta    Funariaceae                Physcomitrium         Physcomitrium_patens                 k__Eukaryota|p__Streptophyta|c__Bryopsida|o__Funariales|f__Funariaceae|g__Physcomitrium|s__Physcomitrium_patens                                                               0.0370789             0.0700601              0.0000000              0.2106432           3
NMR     Eukaryota   Streptophyta    Rubiaceae                  Coffea                Coffea_arabica                       k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Gentianales|f__Rubiaceae|g__Coffea|s__Coffea_arabica                                                                         0.0368507             0.0588294              0.0000000              0.1583226           4
NMR     Eukaryota   Streptophyta    Rosaceae                   Malus                 Malus_sylvestris                     k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Rosales|f__Rosaceae|g__Malus|s__Malus_sylvestris                                                                             0.0363499             0.0837360              0.0000000              0.2484635           2
NMR     Eukaryota   Ascomycota      Debaryomycetaceae          Candida               Candida_albicans                     k__Eukaryota|p__Ascomycota|c__Saccharomycetes|o__Saccharomycetales|f__Debaryomycetaceae|g__Candida|s__Candida_albicans                                                        0.0329433             0.1092607              0.0000000              0.3623767           1
NMR     Eukaryota   Streptophyta    Poaceae                    Poaceae Family        Zea_mays                             k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Poales|f__Poaceae|g__Zea|s__Zea_mays                                                                                         0.0316195             0.0386593              0.0000000              0.0933806           5
NMR     Eukaryota   Streptophyta    Myrtaceae                  Eucalyptus            Eucalyptus_grandis                   k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Myrtales|f__Myrtaceae|g__Eucalyptus|s__Eucalyptus_grandis                                                                    0.0315910             0.0415178              0.0000000              0.1199416           5
NMR     Eukaryota   Streptophyta    Solanaceae                 Lycium                Lycium_ferocissimum                  k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Solanales|f__Solanaceae|g__Lycium|s__Lycium_ferocissimum                                                                     0.0297591             0.0407409              0.0000000              0.1240966           5
NMR     Eukaryota   Streptophyta    Musaceae                   Musa                  Musa_acuminata                       k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Zingiberales|f__Musaceae|g__Musa|s__Musa_acuminata                                                                           0.0272328             0.0408220              0.0000000              0.1037387           4
NMR     Eukaryota   Ascomycota      Sclerotiniaceae            Botrytis              Botrytis_cinerea                     k__Eukaryota|p__Ascomycota|c__Leotiomycetes|o__Helotiales|f__Sclerotiniaceae|g__Botrytis|s__Botrytis_cinerea                                                                  0.0263571             0.0604909              0.0000000              0.1781600           2
NMR     Eukaryota   Streptophyta    Rosaceae                   Malus                 Malus_domestica                      k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Rosales|f__Rosaceae|g__Malus|s__Malus_domestica                                                                              0.0246209             0.0555188              0.0000000              0.1556262           2
NMR     Eukaryota   Apicomplexa     Plasmodiidae               Plasmodium            Plasmodium_knowlesi                  k__Eukaryota|p__Apicomplexa|c__Aconoidasida|o__Haemosporida|f__Plasmodiidae|g__Plasmodium|s__Plasmodium_knowlesi                                                              0.0245241             0.0603256              0.0000000              0.1924194           2
NMR     Eukaryota   Streptophyta    Lythraceae                 Punica                Punica_granatum                      k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Myrtales|f__Lythraceae|g__Punica|s__Punica_granatum                                                                          0.0244833             0.0544891              0.0000000              0.1377269           2
NMR     Eukaryota   Streptophyta    Lamiaceae                  Salvia                Salvia_splendens                     k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Lamiales|f__Lamiaceae|g__Salvia|s__Salvia_splendens                                                                          0.0216213             0.0375662              0.0000000              0.0955285           3
NMR     Eukaryota   Streptophyta    Brassicaceae               Camelina              Camelina_sativa                      k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Brassicales|f__Brassicaceae|g__Camelina|s__Camelina_sativa                                                                   0.0207656             0.0473141              0.0000000              0.1370294           2
NMR     Eukaryota   Streptophyta    Chenopodiaceae             Beta                  Beta_vulgaris                        k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Caryophyllales|f__Chenopodiaceae|g__Beta|s__Beta_vulgaris                                                                    0.0199308             0.0361458              0.0000000              0.1016431           3
NMR     Eukaryota   Streptophyta    Poaceae                    Lolium                Lolium_rigidum                       k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Poales|f__Poaceae|g__Lolium|s__Lolium_rigidum                                                                                0.0193413             0.0331825              0.0000000              0.1114278           5
NMR     Eukaryota   Ascomycota      Aspergillaceae             Aspergillus           Aspergillus_chevalieri               k__Eukaryota|p__Ascomycota|c__Eurotiomycetes|o__Eurotiales|f__Aspergillaceae|g__Aspergillus|s__Aspergillus_chevalieri                                                         0.0173051             0.0397844              0.0000000              0.1280452           3
NMR     Eukaryota   Streptophyta    Fagaceae                   Quercus               Quercus_lobata                       k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Fagales|f__Fagaceae|g__Quercus|s__Quercus_lobata                                                                             0.0167803             0.0556539              0.0000000              0.1845829           1
NMR     Eukaryota   Streptophyta    Solanaceae                 Lycium                Lycium_barbarum                      k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Solanales|f__Solanaceae|g__Lycium|s__Lycium_barbarum                                                                         0.0154722             0.0189632              0.0000000              0.0433267           5
NMR     Eukaryota   Basidiomycota   Pucciniaceae               Puccinia              Puccinia_striiformis                 k__Eukaryota|p__Basidiomycota|c__Pucciniomycetes|o__Pucciniales|f__Pucciniaceae|g__Puccinia|s__Puccinia_striiformis                                                           0.0106036             0.0117712              0.0000000              0.0285213           6
NMR     Eukaryota   Streptophyta    Fabaceae                   Phaseolus             Phaseolus_vulgaris                   k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Fabales|f__Fabaceae|g__Phaseolus|s__Phaseolus_vulgaris                                                                       0.0093443             0.0309914              0.0000000              0.1027868           1
NMR     Eukaryota   Streptophyta    Poaceae                    Oryza                 Oryza_glaberrima                     k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Poales|f__Poaceae|g__Oryza|s__Oryza_glaberrima                                                                               0.0084996             0.0281901              0.0000000              0.0934959           1
NMR     Eukaryota   Fornicata       Hexamitidae                Giardia               Giardia_intestinalis                 k__Eukaryota|p__Fornicata|NA|o__Diplomonadida|f__Hexamitidae|g__Giardia|s__Giardia_intestinalis                                                                               0.0084750             0.0184802              0.0000000              0.0634814           7
NMR     Eukaryota   Chlorophyta     Bathycoccaceae             Ostreococcus          Ostreococcus_sp._'lucimarinus'       k__Eukaryota|p__Chlorophyta|c__Mamiellophyceae|o__Mamiellales|f__Bathycoccaceae|g__Ostreococcus|s__Ostreococcus_sp._lucimarinus                                               0.0068757             0.0228042              0.0000000              0.0756331           1
NMR     Eukaryota   Streptophyta    Poaceae                    Panicum               Panicum_hallii                       k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Poales|f__Poaceae|g__Panicum|s__Panicum_hallii                                                                               0.0067167             0.0136447              0.0000000              0.0435503           3
NMR     Eukaryota   Streptophyta    Rosaceae                   Rosa                  Rosa_chinensis                       k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Rosales|f__Rosaceae|g__Rosa|s__Rosa_chinensis                                                                                0.0058685             0.0136958              0.0000000              0.0415238           2
NMR     Eukaryota   Basidiomycota   Pucciniaceae               Puccinia              Puccinia_triticina                   k__Eukaryota|p__Basidiomycota|c__Pucciniomycetes|o__Pucciniales|f__Pucciniaceae|g__Puccinia|s__Puccinia_triticina                                                             0.0043004             0.0066941              0.0000000              0.0180182           4
NMR     Eukaryota   Streptophyta    Arecaceae                  Phoenix               Phoenix_dactylifera                  k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Arecales|f__Arecaceae|g__Phoenix|s__Phoenix_dactylifera                                                                      0.0040302             0.0133667              0.0000000              0.0443324           1
NMR     Eukaryota   Ascomycota      Schizosaccharomycetaceae   Schizosaccharomyces   Schizosaccharomyces_pombe            k__Eukaryota|p__Ascomycota|c__Schizosaccharomycetes|o__Schizosaccharomycetales|f__Schizosaccharomycetaceae|g__Schizosaccharomyces|s__Schizosaccharomyces_pombe                0.0034099             0.0060928              0.0000000              0.0169310           3
NMR     Eukaryota   Streptophyta    Poaceae                    Lolium                Lolium_perenne                       k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Poales|f__Poaceae|g__Lolium|s__Lolium_perenne                                                                                0.0027017             0.0089605              0.0000000              0.0297187           1
NMR     Eukaryota   Ascomycota      Mycosphaerellaceae         Zymoseptoria          Zymoseptoria_tritici                 k__Eukaryota|p__Ascomycota|c__Dothideomycetes|o__Mycosphaerellales|f__Mycosphaerellaceae|g__Zymoseptoria|s__Zymoseptoria_tritici                                              0.0025170             0.0056850              0.0000000              0.0160337           2
NMR     Eukaryota   Streptophyta    Poaceae                    Brachypodium          Brachypodium_distachyon              k__Eukaryota|p__Streptophyta|c__Magnoliopsida|o__Poales|f__Poaceae|g__Brachypodium|s__Brachypodium_distachyon                                                                 0.0010741             0.0035624              0.0000000              0.0118152           1
NMR     Eukaryota   Apicomplexa     Sarcocystidae              Neospora              Neospora_caninum                     k__Eukaryota|p__Apicomplexa|c__Conoidasida|o__Eucoccidiorida|f__Sarcocystidae|g__Neospora|s__Neospora_caninum                                                                 0.0009016             0.0029902              0.0000000              0.0099174           1
NMR     Eukaryota   Ascomycota      Debaryomycetaceae          Debaryomyces          Debaryomyces_hansenii                k__Eukaryota|p__Ascomycota|c__Saccharomycetes|o__Saccharomycetales|f__Debaryomycetaceae|g__Debaryomyces|s__Debaryomyces_hansenii                                              0.0008571             0.0028426              0.0000000              0.0094279           1
NMR     Eukaryota   Ascomycota      Dipodascaceae              Yarrowia              Yarrowia_lipolytica                  k__Eukaryota|p__Ascomycota|c__Saccharomycetes|o__Saccharomycetales|f__Dipodascaceae|g__Yarrowia|s__Yarrowia_lipolytica                                                        0.0003299             0.0010942              0.0000000              0.0036289           1

Eukaryotic phyla


``` r
mean_sd_relab.nmr_no_groups.asv.nonbact%>%
  filter(Kingdom=="Eukaryota")%>%
  group_by(Phylum)%>%
  distinct(Species,.keep_all = T)%>%
  tally()%>%
  knitr::kable(format = "simple")
```



Phylum            n
--------------  ---
Apicomplexa       2
Ascomycota        9
Basidiomycota     2
Chlorophyta       1
Fornicata         1
Parabasalia       1
Streptophyta     84

Eukaryotic species that aren't Streptophyta


``` r
mean_sd_relab.nmr_no_groups.asv.nonbact%>%
  filter(Kingdom=="Eukaryota",Phylum!="Streptophyta")%>%
  knitr::kable(format = "simple")
```



class   Kingdom     Phylum          Family                     Genus                 Species                          OTU                                                                                                                                                               MeanRelativeAbundance   sdRelativeAbundance   minRelativeAbundance   maxRelativeAbundance   n_samples
------  ----------  --------------  -------------------------  --------------------  -------------------------------  ---------------------------------------------------------------------------------------------------------------------------------------------------------------  ----------------------  --------------------  ---------------------  ---------------------  ----------
NMR     Eukaryota   Parabasalia     Trichomonadidae            Trichomonas           Trichomonas_vaginalis            k__Eukaryota|p__Parabasalia|NA|o__Trichomonadida|f__Trichomonadidae|g__Trichomonas|s__Trichomonas_vaginalis                                                                   0.6564723             0.5321252              0.1328118              1.7745074          11
NMR     Eukaryota   Ascomycota      Pichiaceae                 Pichia                Pichia_kudriavzevii              k__Eukaryota|p__Ascomycota|c__Saccharomycetes|o__Saccharomycetales|f__Pichiaceae|g__Pichia|s__Pichia_kudriavzevii                                                             0.5206977             0.9007816              0.0000000              2.3987543           8
NMR     Eukaryota   Ascomycota      Aspergillaceae             Aspergillus           Aspergillus_puulaauensis         k__Eukaryota|p__Ascomycota|c__Eurotiomycetes|o__Eurotiales|f__Aspergillaceae|g__Aspergillus|s__Aspergillus_puulaauensis                                                       0.3496001             0.5269957              0.0000000              1.7242715           8
NMR     Eukaryota   Ascomycota      Debaryomycetaceae          Candida               Candida_albicans                 k__Eukaryota|p__Ascomycota|c__Saccharomycetes|o__Saccharomycetales|f__Debaryomycetaceae|g__Candida|s__Candida_albicans                                                        0.0329433             0.1092607              0.0000000              0.3623767           1
NMR     Eukaryota   Ascomycota      Sclerotiniaceae            Botrytis              Botrytis_cinerea                 k__Eukaryota|p__Ascomycota|c__Leotiomycetes|o__Helotiales|f__Sclerotiniaceae|g__Botrytis|s__Botrytis_cinerea                                                                  0.0263571             0.0604909              0.0000000              0.1781600           2
NMR     Eukaryota   Apicomplexa     Plasmodiidae               Plasmodium            Plasmodium_knowlesi              k__Eukaryota|p__Apicomplexa|c__Aconoidasida|o__Haemosporida|f__Plasmodiidae|g__Plasmodium|s__Plasmodium_knowlesi                                                              0.0245241             0.0603256              0.0000000              0.1924194           2
NMR     Eukaryota   Ascomycota      Aspergillaceae             Aspergillus           Aspergillus_chevalieri           k__Eukaryota|p__Ascomycota|c__Eurotiomycetes|o__Eurotiales|f__Aspergillaceae|g__Aspergillus|s__Aspergillus_chevalieri                                                         0.0173051             0.0397844              0.0000000              0.1280452           3
NMR     Eukaryota   Basidiomycota   Pucciniaceae               Puccinia              Puccinia_striiformis             k__Eukaryota|p__Basidiomycota|c__Pucciniomycetes|o__Pucciniales|f__Pucciniaceae|g__Puccinia|s__Puccinia_striiformis                                                           0.0106036             0.0117712              0.0000000              0.0285213           6
NMR     Eukaryota   Fornicata       Hexamitidae                Giardia               Giardia_intestinalis             k__Eukaryota|p__Fornicata|NA|o__Diplomonadida|f__Hexamitidae|g__Giardia|s__Giardia_intestinalis                                                                               0.0084750             0.0184802              0.0000000              0.0634814           7
NMR     Eukaryota   Chlorophyta     Bathycoccaceae             Ostreococcus          Ostreococcus_sp._'lucimarinus'   k__Eukaryota|p__Chlorophyta|c__Mamiellophyceae|o__Mamiellales|f__Bathycoccaceae|g__Ostreococcus|s__Ostreococcus_sp._lucimarinus                                               0.0068757             0.0228042              0.0000000              0.0756331           1
NMR     Eukaryota   Basidiomycota   Pucciniaceae               Puccinia              Puccinia_triticina               k__Eukaryota|p__Basidiomycota|c__Pucciniomycetes|o__Pucciniales|f__Pucciniaceae|g__Puccinia|s__Puccinia_triticina                                                             0.0043004             0.0066941              0.0000000              0.0180182           4
NMR     Eukaryota   Ascomycota      Schizosaccharomycetaceae   Schizosaccharomyces   Schizosaccharomyces_pombe        k__Eukaryota|p__Ascomycota|c__Schizosaccharomycetes|o__Schizosaccharomycetales|f__Schizosaccharomycetaceae|g__Schizosaccharomyces|s__Schizosaccharomyces_pombe                0.0034099             0.0060928              0.0000000              0.0169310           3
NMR     Eukaryota   Ascomycota      Mycosphaerellaceae         Zymoseptoria          Zymoseptoria_tritici             k__Eukaryota|p__Ascomycota|c__Dothideomycetes|o__Mycosphaerellales|f__Mycosphaerellaceae|g__Zymoseptoria|s__Zymoseptoria_tritici                                              0.0025170             0.0056850              0.0000000              0.0160337           2
NMR     Eukaryota   Apicomplexa     Sarcocystidae              Neospora              Neospora_caninum                 k__Eukaryota|p__Apicomplexa|c__Conoidasida|o__Eucoccidiorida|f__Sarcocystidae|g__Neospora|s__Neospora_caninum                                                                 0.0009016             0.0029902              0.0000000              0.0099174           1
NMR     Eukaryota   Ascomycota      Debaryomycetaceae          Debaryomyces          Debaryomyces_hansenii            k__Eukaryota|p__Ascomycota|c__Saccharomycetes|o__Saccharomycetales|f__Debaryomycetaceae|g__Debaryomyces|s__Debaryomyces_hansenii                                              0.0008571             0.0028426              0.0000000              0.0094279           1
NMR     Eukaryota   Ascomycota      Dipodascaceae              Yarrowia              Yarrowia_lipolytica              k__Eukaryota|p__Ascomycota|c__Saccharomycetes|o__Saccharomycetales|f__Dipodascaceae|g__Yarrowia|s__Yarrowia_lipolytica                                                        0.0003299             0.0010942              0.0000000              0.0036289           1



### 13.2 Archaea
How many archaea do we have 


``` r
mean_sd_relab.nmr_no_groups.asv%>%
  filter(Kingdom=="Archaea")%>%
  distinct(Species)%>%
  tally()
```

```
## # A tibble: 1 × 1
##       n
##   <int>
## 1    15
```

Total abundance of archaea in classified data


``` r
mean_sd_relab.nmr_no_groups.asv.nonbact%>%
  ungroup%>%
  filter(Kingdom=="Archaea")%>%
  distinct()%>%
  knitr::kable(format = "simple")
```



class   Kingdom   Phylum                        Family                     Genus                Species                              OTU                                                                                                                                                                         MeanRelativeAbundance   sdRelativeAbundance   minRelativeAbundance   maxRelativeAbundance   n_samples
------  --------  ----------------------------  -------------------------  -------------------  -----------------------------------  -------------------------------------------------------------------------------------------------------------------------------------------------------------------------  ----------------------  --------------------  ---------------------  ---------------------  ----------
NMR     Archaea   Euryarchaeota                 Methanobacteriaceae        Methanobrevibacter   Methanobrevibacter_sp._AbM4          k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanobrevibacter|s__Methanobrevibacter_sp._AbM4                                        1.2808634             1.6381718              0.0261897              5.8973226          11
NMR     Archaea   Euryarchaeota                 Methanobacteriaceae        Methanobrevibacter   Methanobrevibacter_millerae          k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanobrevibacter|s__Methanobrevibacter_millerae                                        0.2749037             0.3953089              0.0000000              1.3942289          10
NMR     Archaea   Euryarchaeota                 Methanobacteriaceae        Methanobrevibacter   Methanobrevibacter_ruminantium       k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanobrevibacter|s__Methanobrevibacter_ruminantium                                     0.2685322             0.3095590              0.0000000              1.1253462          10
NMR     Archaea   Euryarchaeota                 Methanobacteriaceae        Methanobrevibacter   Methanobrevibacter_olleyae           k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanobrevibacter|s__Methanobrevibacter_olleyae                                         0.2068656             0.3532260              0.0000000              1.1910172           7
NMR     Archaea   Euryarchaeota                 Methanobacteriaceae        Methanobrevibacter   Methanobrevibacter_smithii           k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanobrevibacter|s__Methanobrevibacter_smithii                                         0.1641641             0.2581924              0.0000000              0.8944981           9
NMR     Archaea   Euryarchaeota                 Methanobacteriaceae        Methanobrevibacter   Methanobrevibacter_arboriphilus      k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanobrevibacter|s__Methanobrevibacter_arboriphilus                                    0.1226819             0.1979728              0.0000000              0.6511967           6
NMR     Archaea   Euryarchaeota                 Methanobacteriaceae        Methanobrevibacter   Methanobrevibacter_sp._YE315         k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanobrevibacter|s__Methanobrevibacter_sp._YE315                                       0.0923706             0.1405313              0.0000000              0.4923608          10
NMR     Archaea   Euryarchaeota                 Methanobacteriaceae        Methanobrevibacter   Methanobrevibacter_sp._TLL-48-HuF1   k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanobrevibacter|s__Methanobrevibacter_sp._TLL-48-HuF1                                 0.0725664             0.0968016              0.0000000              0.3429399           9
NMR     Archaea   Euryarchaeota                 Methanobacteriaceae        Methanosphaera       Methanosphaera_stadtmanae            k__Archaea|p__Euryarchaeota|c__Methanobacteria|o__Methanobacteriales|f__Methanobacteriaceae|g__Methanosphaera|s__Methanosphaera_stadtmanae                                              0.0069179             0.0229441              0.0000000              0.0760971           1
NMR     Archaea   Nitrososphaerota              Nitrosopumilaceae          Nitrosomarinus       Candidatus_Nitrosomarinus_catalina   k__Archaea|p__Nitrososphaerota|c__Nitrososphaeria|o__Nitrosopumilales|f__Nitrosopumilaceae|g__Nitrosomarinus|s__Candidatus_Nitrosomarinus_catalina                                      0.0020272             0.0059048              0.0000000              0.0196714           2
NMR     Archaea   Euryarchaeota                 Haloarculaceae             Halomicroarcula      Halomicroarcula_marina               k__Archaea|p__Euryarchaeota|c__Halobacteria|o__Halobacteriales|f__Haloarculaceae|g__Halomicroarcula|s__Halomicroarcula_marina                                                           0.0006524             0.0011282              0.0000000              0.0037712           5
NMR     Archaea   Candidatus_Thermoplasmatota   Methanomassiliicoccaceae   Methanoplasma        Candidatus_Methanoplasma_termitum    k__Archaea|p__Candidatus_Thermoplasmatota|c__Thermoplasmata|o__Methanomassiliicoccales|f__Methanomassiliicoccaceae|g__Methanoplasma|s__Candidatus_Methanoplasma_termitum                0.0001795             0.0002953              0.0000000              0.0009550           5
NMR     Archaea   Thermoproteota                Desulfurococcaceae         Aeropyrum            Aeropyrum_camini                     k__Archaea|p__Thermoproteota|c__Thermoprotei|o__Desulfurococcales|f__Desulfurococcaceae|g__Aeropyrum|s__Aeropyrum_camini                                                                0.0001182             0.0003693              0.0000000              0.0012309           3
NMR     Archaea   Euryarchaeota                 Haloarculaceae             Halapricum           Halapricum_desulfuricans             k__Archaea|p__Euryarchaeota|c__Halobacteria|o__Halobacteriales|f__Haloarculaceae|g__Halapricum|s__Halapricum_desulfuricans                                                              0.0000637             0.0001479              0.0000000              0.0004444           2
NMR     Archaea   Euryarchaeota                 Methanomicrobiaceae        Methanofollis        Methanofollis_formosanus             k__Archaea|p__Euryarchaeota|c__Methanomicrobia|o__Methanomicrobiales|f__Methanomicrobiaceae|g__Methanofollis|s__Methanofollis_formosanus                                                0.0000167             0.0000553              0.0000000              0.0001834           1

Archaea phyla


``` r
mean_sd_relab.nmr_no_groups.asv.nonbact%>%
  filter(Kingdom=="Archaea")%>%
  group_by(Phylum)%>%
  distinct(Species,.keep_all = T)%>%
  tally()%>%
  knitr::kable(format = "simple")
```



Phylum                          n
----------------------------  ---
Candidatus_Thermoplasmatota     1
Euryarchaeota                  12
Nitrososphaerota                1
Thermoproteota                  1



### 13.3 Viruses
How many viruses do we have


``` r
mean_sd_relab.nmr_no_groups.asv.nonbact%>%
  filter(Kingdom=="Viruses")%>%
  distinct(Species)%>%
  tally()
```

```
## # A tibble: 1 × 1
##       n
##   <int>
## 1    12
```

Total abundance of viruses in classified data


``` r
mean_sd_relab.nmr_no_groups.asv.nonbact%>%
  ungroup%>%
  filter(Kingdom=="Viruses")%>%
  distinct()%>%
  knitr::kable(format = "simple")
```



class   Kingdom   Phylum        Family           Genus         Species                   OTU                                                                                                                     MeanRelativeAbundance   sdRelativeAbundance   minRelativeAbundance   maxRelativeAbundance   n_samples
------  --------  ------------  ---------------  ------------  ------------------------  ---------------------------------------------------------------------------------------------------------------------  ----------------------  --------------------  ---------------------  ---------------------  ----------
NMR     Viruses   Uroviricota   Suoliviridae     Buchavirus    Buchavirus_intestinalis   k__Viruses|p__Uroviricota|c__Caudoviricetes|o__Crassvirales|f__Suoliviridae|g__Buchavirus|s__Buchavirus_intestinalis                0.0075926             0.0144536                      0              0.0470114           6
NMR     Viruses   Uroviricota   Suoliviridae     Buchavirus    Buchavirus_oralis         k__Viruses|p__Uroviricota|c__Caudoviricetes|o__Crassvirales|f__Suoliviridae|g__Buchavirus|s__Buchavirus_oralis                      0.0007431             0.0015452                      0              0.0049456           3
NMR     Viruses   Uroviricota   Suoliviridae     Buchavirus    Buchavirus_splanchnicus   k__Viruses|p__Uroviricota|c__Caudoviricetes|o__Crassvirales|f__Suoliviridae|g__Buchavirus|s__Buchavirus_splanchnicus                0.0006591             0.0014414                      0              0.0046910           3
NMR     Viruses   Uroviricota   Suoliviridae     Buchavirus    Buchavirus_coli           k__Viruses|p__Uroviricota|c__Caudoviricetes|o__Crassvirales|f__Suoliviridae|g__Buchavirus|s__Buchavirus_coli                        0.0005628             0.0011747                      0              0.0037054           3
NMR     Viruses   Uroviricota   Suoliviridae     Buchavirus    Buchavirus_hiberniae      k__Viruses|p__Uroviricota|c__Caudoviricetes|o__Crassvirales|f__Suoliviridae|g__Buchavirus|s__Buchavirus_hiberniae                   0.0004616             0.0009761                      0              0.0031248           3
NMR     Viruses   Uroviricota   Salasmaviridae   Claudivirus   Claudivirus_stitch        k__Viruses|p__Uroviricota|c__Caudoviricetes|NA|f__Salasmaviridae|g__Claudivirus|s__Claudivirus_stitch                               0.0000660             0.0001842                      0              0.0006118           2
NMR     Viruses   Uroviricota   Suoliviridae     Canhaevirus   Canhaevirus_hiberniae     k__Viruses|p__Uroviricota|c__Caudoviricetes|o__Crassvirales|f__Suoliviridae|g__Canhaevirus|s__Canhaevirus_hiberniae                 0.0000635             0.0001627                      0              0.0005296           2
NMR     Viruses   Uroviricota   Suoliviridae     Buchavirus    Buchavirus_hominis        k__Viruses|p__Uroviricota|c__Caudoviricetes|o__Crassvirales|f__Suoliviridae|g__Buchavirus|s__Buchavirus_hominis                     0.0000448             0.0000880                      0              0.0002649           3
NMR     Viruses   Uroviricota   Suoliviridae     Cohcovirus    Cohcovirus_splanchnicus   k__Viruses|p__Uroviricota|c__Caudoviricetes|o__Crassvirales|f__Suoliviridae|g__Cohcovirus|s__Cohcovirus_splanchnicus                0.0000110             0.0000365                      0              0.0001211           1
NMR     Viruses   Uroviricota   Suoliviridae     Buchavirus    Buchavirus_faecalis       k__Viruses|p__Uroviricota|c__Caudoviricetes|o__Crassvirales|f__Suoliviridae|g__Buchavirus|s__Buchavirus_faecalis                    0.0000105             0.0000348                      0              0.0001155           1
NMR     Viruses   Uroviricota   Suoliviridae     Cohcovirus    Cohcovirus_hiberniae      k__Viruses|p__Uroviricota|c__Caudoviricetes|o__Crassvirales|f__Suoliviridae|g__Cohcovirus|s__Cohcovirus_hiberniae                   0.0000061             0.0000204                      0              0.0000676           1
NMR     Viruses   Uroviricota   Suoliviridae     Buchavirus    Buchavirus_copri          k__Viruses|p__Uroviricota|c__Caudoviricetes|o__Crassvirales|f__Suoliviridae|g__Buchavirus|s__Buchavirus_copri                       0.0000023             0.0000077                      0              0.0000255           1

Viral phyla 


``` r
mean_sd_relab.nmr_no_groups.asv.nonbact%>%
  filter(Kingdom=="Viruses")%>%
  group_by(Phylum)%>%
  distinct(Species,.keep_all = T)%>%
  tally()%>%
  knitr::kable(format = "simple")
```



Phylum          n
------------  ---
Uroviricota    12

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
##  [1] vegan_2.6-4     lattice_0.22-6  permute_0.9-8   ggtext_0.1.2   
##  [5] lubridate_1.9.4 forcats_1.0.0   stringr_1.5.1   dplyr_1.1.4    
##  [9] purrr_1.0.4     readr_2.1.5     tidyr_1.3.1     tibble_3.2.1   
## [13] ggplot2_4.0.0   tidyverse_2.0.0 phyloseq_1.50.0
## 
## loaded via a namespace (and not attached):
##  [1] ade4_1.7-23             tidyselect_1.2.1        viridisLite_0.4.3      
##  [4] farver_2.1.2            Biostrings_2.74.1       S7_0.2.0               
##  [7] fastmap_1.2.0           digest_0.6.37           timechange_0.3.0       
## [10] lifecycle_1.0.5         cluster_2.1.8.1         survival_3.8-3         
## [13] magrittr_2.0.3          compiler_4.4.3          rlang_1.1.5            
## [16] sass_0.4.10             tools_4.4.3             igraph_2.1.4           
## [19] yaml_2.3.12             data.table_1.17.8       knitr_1.51             
## [22] labeling_0.4.3          here_1.0.2              xml2_1.3.8             
## [25] plyr_1.8.9              RColorBrewer_1.1-3      withr_3.0.2            
## [28] BiocGenerics_0.52.0     grid_4.4.3              stats4_4.4.3           
## [31] multtest_2.62.0         biomformat_1.34.0       colorspace_2.1-2       
## [34] Rhdf5lib_1.28.0         scales_1.4.0            iterators_1.0.14       
## [37] MASS_7.3-65             cli_3.6.4               rmarkdown_2.30         
## [40] crayon_1.5.3            ragg_1.5.0              generics_0.1.4         
## [43] otel_0.2.0              rstudioapi_0.18.0       httr_1.4.7             
## [46] reshape2_1.4.4          tzdb_0.5.0              commonmark_2.0.0       
## [49] ape_5.8-1               cachem_1.1.0            rhdf5_2.50.2           
## [52] zlibbioc_1.52.0         splines_4.4.3           parallel_4.4.3         
## [55] XVector_0.46.0          vctrs_0.6.5             Matrix_1.7-4           
## [58] jsonlite_2.0.0          litedown_0.9            bookdown_0.46          
## [61] IRanges_2.40.1          hms_1.1.4               S4Vectors_0.44.0       
## [64] systemfonts_1.2.3       foreach_1.5.2           jquerylib_0.1.4        
## [67] glue_1.8.0              codetools_0.2-20        stringi_1.8.4          
## [70] gtable_0.3.6            GenomeInfoDb_1.42.3     UCSC.utils_1.2.0       
## [73] pillar_1.11.1           htmltools_0.5.8.1       rhdf5filters_1.18.1    
## [76] GenomeInfoDbData_1.2.13 R6_2.6.1                textshaping_1.0.3      
## [79] rprojroot_2.1.1         evaluate_1.0.5          Biobase_2.66.0         
## [82] markdown_2.0            gridtext_0.1.5          bslib_0.10.0           
## [85] Rcpp_1.0.14             nlme_3.1-167            mgcv_1.9-1             
## [88] xfun_0.56               pkgconfig_2.0.3
```

``` r
rm(list =setdiff(ls(all.names = TRUE), c("markdown.dir")))
gc()
```

```
##           used  (Mb) gc trigger  (Mb) max used  (Mb)
## Ncells 5634095 300.9   10688262 570.9 10688262 570.9
## Vcells 9797262  74.8   18426857 140.6 15287820 116.7
```

