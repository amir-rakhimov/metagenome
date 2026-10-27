---
output: 
  bookdown::html_document2:
     toc: true
---




# Differential microbial abundance tests with MaAsLin2




## Introduction
This script performs differential microbial abundance tests on 
different hosts. We will use MaAsLin2.




## Load necessary libraries.


``` r
# if(!requireNamespace("BiocManager")){
#   install.packages("BiocManager")
# }
# Current: Bioconductor version 3.20 (BiocManager 1.30.27), R 4.4.3 (2025-02-28 ucrt)
```

Maaslin2: 1.20.0


``` r
# BiocManager::install(c("Maaslin2","phyloseq"), version = "3.20")
# BiocManager::install("phyloseq", version = "3.17")
# install.packages(c("tidyverse","vegan"))
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
library(Maaslin2)
library(vegan)
```

```
## Loading required package: permute
## Loading required package: lattice
## This is vegan 2.6-4
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



## Specifying parameters and directory/file names. 


``` r
source(here::here("config/R/config.R"))# config file with global variables
source(here::here("config/R/themes.R"))# config file with themes
dir.create(diffabund.tables,recursive = TRUE)
```

```
## Warning in dir.create(diffabund.tables, recursive = TRUE):
## '\\wsl.localhost\Ubuntu-22.04\home\rakhimov\projects\metagenome\analysis\1-kraken2\20240409_17_32_40\02-differential_abundance\tables'
## already exists
```

``` r
dir.create(diffabund.rdafiles,recursive = TRUE)
```

```
## Warning in dir.create(diffabund.rdafiles, recursive = TRUE):
## '\\wsl.localhost\Ubuntu-22.04\home\rakhimov\projects\metagenome\analysis\1-kraken2\20240409_17_32_40\02-differential_abundance\rdafiles'
## already exists
```

``` r
dir.create(diffabund.figures,recursive = TRUE)
```

```
## Warning in dir.create(diffabund.figures, recursive = TRUE):
## '\\wsl.localhost\Ubuntu-22.04\home\rakhimov\projects\metagenome\analysis\1-kraken2\20240409_17_32_40\02-differential_abundance\figures'
## already exists
```

``` r
authorname <- "kraken2"
rare.status <- "unrare"
```

Import the phyloseq object


``` r
ps.q <- readRDS(file= ps.q.raw.fname)
ps.q.rel<-readRDS(file = ps.q.rel.raw.fname)


### Create reference levels. ####
prepare_data <- function(ps,
                         agglom.rank, comparison, 
                         ref.level, metadata,
                         relations = NULL){
  message("Preparing data for differential microbial abundance tests on ", comparison, " variable at ", agglom.rank, " level. ",
          "Reference: ", ref.level)
  comparison_meta <- list(
    age = list(custom.levels = "NMR",
               maaslin2.reference = paste("agegroup", ref.level, sep = ","),
               maaslin2.fixed_effects = "agegroup",
               maaslin2.random_effects = c("relation"),
               maaslin2.random_effects.file_suffix = "with-relations",
               ggplot.levels = c("agegroup0_10",
                                 "agegroup10_16" ),
               ggplot.labels = c("agegroup0_10" = "agegroup0_10",
                                 "agegroup10_16" = "agegroup10_16"),
               gg.labs.name = "Age group",
               gg.title.groups = "age groups",
               output.filename = paste("NMR",agglom.rank,comparison,
                                       "ref",ref.level,
                                       sep="-")
    ),
    sex = list(custom.levels = "NMR",
               maaslin2.reference = paste("sex", ref.level, sep = ","),
               maaslin2.fixed_effects = "sex" ,
               maaslin2.random_effects = c("relation"),
               maaslin2.random_effects.file_suffix = "with-relations",
               ggplot.levels = c("female" ,
                                 "male" ),
               ggplot.labels = c("female" = "Females",
                                 "male" = "Males"),
               gg.labs.name = "Host sex",
               gg.title.groups = "groups",
               output.filename = paste("NMR",agglom.rank,comparison,
                                       "ref",ref.level,
                                       sep="-")
    )
  )
  current.comparison_meta <- comparison_meta[[comparison]]
  ### Add relations metadata ####
  if (setequal(current.comparison_meta$custom.levels,"NMR")){
    metadata<-metadata %>%
      right_join(relations,by="Sample")
    rownames(metadata)<-metadata$Sample
  }
  # add metadata to the list
  current.comparison_meta$metadata <- metadata
  
  # add sample groups variable
  if(setequal(current.comparison_meta$custom.levels,"NMR")){
    if(comparison=="age"){
      current.comparison_meta$sample.groups <-c("agegroup0_10",
                                                "agegroup10_16")
    }else if(comparison=="sex"){
      current.comparison_meta$sample.groups<-c("female",
                                               "male")
    }
  }else{
    current.comparison_meta$sample.groups <- current.comparison_meta$custom.levels
  }
  
  #' Relevel the comparison vector for ANCOM-BC. The first level will be the reference 
  #' for custom leveling
  current.comparison_meta$ancombc.levels <- 
    c(ref.level, 
      current.comparison_meta$sample.groups[current.comparison_meta$sample.groups!=ref.level])
  
  # Create a phyloseq with relative abundances
  ps.rel <- transform_sample_counts(ps, function(x) 100*x/sum(x))
  ### Melt the phyloseq object. ####
  if(agglom.rank == "OTU"){
    ps.agg <- ps%>%
      psmelt()
    
    ps.rel.agg <- ps.rel%>%
      psmelt()
  }else{
    ps.agg <- ps%>%
      tax_glom(taxrank = agglom.rank)%>%
      psmelt()%>%
      dplyr::select(-OTU)
    ps.rel.agg <- ps.rel%>%
      tax_glom(taxrank = agglom.rank)%>%
      psmelt()%>%
      dplyr::select(-OTU)
  }
  
  ps.agg.wide <- ps.agg %>%
    dplyr::select(all_of(c("Sample", "Abundance", agglom.rank)))%>%
    pivot_wider(names_from = agglom.rank, # or OTU
                values_from = "Abundance",
                values_fill = 0)%>%
    as.data.frame()%>%
    column_to_rownames("Sample")
  
  current.comparison_meta$ps.agg <- ps.agg
  current.comparison_meta$ps.rel.agg <- ps.rel.agg
  current.comparison_meta$ps.agg.wide <- ps.agg.wide
  
  return(current.comparison_meta)
}
```

Function that performs a test with MaAsLin2:


``` r
perform_maaslin2_test <- function(ps.agg.wide, 
                                  output.filename,
                                  comparison,
                                  metadata,
                                  ref.level,
                                  maaslin2.random_effects,
                                  maaslin2.random_effects.file_suffix,
                                  maaslin2.fixed_effects,
                                  maaslin2.reference
){
  #' Maaslin2 test
  maaslin2.main.dir <- file.path(maaslin2.out,
                                 paste0(authorname,"-output"),
                                 rare.status)
  if(!dir.exists(maaslin2.main.dir)){
    dir.create(maaslin2.main.dir,recursive = TRUE)
  }
  maaslin2.output.path <- file.path(maaslin2.out,
                                    paste0(authorname,"-output"),
                                    rare.status,paste(
                                      output.filename, 
                                      maaslin2.random_effects.file_suffix,
                                      sep = "-"))
  fit.rds.fname <- file.path(diffabund.rdafiles,
                             paste(
                               paste("maaslin.fit_data",comparison, sep = "."),
                               paste0(output.filename,".rds"),sep = "-"))
  
  if(dir.exists(maaslin2.output.path) & file.exists(fit.rds.fname)){
    message(maaslin2.output.path, " exists. Skipping")
    maaslin.fit_data <- readRDS(fit.rds.fname)
  }else{
    set.seed(1)
    maaslin.fit_data = 
      Maaslin2(input_data = ps.agg.wide, 
               input_metadata = metadata, 
               min_prevalence = 0,
               normalization = "TSS",
               transform = "LOG",
               analysis_method = "LM",
               random_effects = maaslin2.random_effects, 
               standardize = FALSE,
               output = maaslin2.output.path, 
               fixed_effects = c(maaslin2.fixed_effects),
               reference = maaslin2.reference,
               max_significance = 0.05,
               save_models = TRUE)
    #' Save the fit data object as an rds file.
    if(!file.exists(fit.rds.fname)){
      saveRDS(maaslin.fit_data,
              file = fit.rds.fname)
      message("Saved maaslin.fit_data as ", fit.rds.fname)
    }
    
  }
  return(maaslin.fit_data)
}
```

Function to process the MaAsLin2 output:


``` r
process_maaslin2_output<-function(maaslin.fit_data,
                                  ps.agg,
                                  sample.groups,
                                  agglom.rank,
                                  output.filename){
  #' Downstream processing of Maaslin2 output.
  source(file.path(util.functions.r,"make_features_maaslin.R"))
  #' Extract features with qvalue < 0.05.
  if(min(maaslin.fit_data$results$qval)<0.05){
    maaslin2.signif.features<-maaslin.fit_data$results%>%
      filter(qval<0.05) # should be qval
  }else{
    maaslin2.signif.features<-maaslin.fit_data$results%>%
      arrange(qval)%>%
      head(n = 15) # if no significant results found
  }
  
  #' Make features pretty. We have to make this exchange because maaslin 
  #' output treats space and hyphen as the same thing.
  if(agglom.rank=="OTU"){
    maaslin2.signif.features$feature<-gsub("^X","",
                                           maaslin2.signif.features$feature)
  }else{
    foo<-ps.agg
    foo<-foo%>%mutate("maaslin"=get(agglom.rank))
    foo<-make_features_maaslin(foo,"maaslin")
    foo<-unique(foo[,c("maaslin",agglom.rank)])
    maaslin2.signif.features<-maaslin2.signif.features%>%
      left_join(foo[,c("maaslin",agglom.rank)],by=c("feature"="maaslin"))%>%
      distinct()
    rm(foo)
    maaslin2.signif.features$feature<-maaslin2.signif.features[,agglom.rank]
    maaslin2.signif.features<-subset(maaslin2.signif.features, 
                                     select=-get(agglom.rank))
  }
  
  #' Extract features that are downregulated in all hosts simultaneously.
  # Downreglated features have coef<0.
  # We also add association strength column
  maaslin2.signif.decreased<-maaslin2.signif.features%>%
    as_tibble()%>%
    filter(coef<0)%>%
    group_by(feature)%>%
    filter(n()==length(sample.groups)-1)%>% 
    arrange(pval,feature)%>%
    mutate(n=n())%>%
    mutate(assoc.str=-log(qval)*sign(coef))%>%
    ungroup()
  
  
  #' Check if all features (OTU/taxa/etc) are also present in the 
  #' ps.agg dataframe.
  table(maaslin2.signif.decreased$feature%in%ps.agg[,agglom.rank])%>%
    print()
  
  signif.features.fname <-file.path(diffabund.tables,paste(
    "maaslin2",output.filename,
    "signif.tsv",sep="-"))
  signif.decreased.fname <- file.path(diffabund.tables,
                                      paste("maaslin2.signif.decreased",
                                            output.filename,
                                            "signif.tsv",sep="-"))
  
  if(file.exists(signif.features.fname)){
    message(signif.features.fname, " exists. Skipping")
  }else{
    message("Saving maaslin2.signif.features as ", signif.features.fname)
    write.table(maaslin2.signif.features,
                file = signif.features.fname,
                row.names = F,sep = "\t")
    
  }
  if(file.exists(signif.decreased.fname)){
    message(signif.decreased.fname, " exists. Skipping")
  }else{
    message("Saving maaslin2.signif.decreased as ", signif.decreased.fname)
    write.table(maaslin2.signif.decreased,
                file = signif.decreased.fname,
                row.names = F,sep = "\t")
    
  }
  return(list(maaslin2.signif.features = maaslin2.signif.features,
              maaslin2.signif.decreased = maaslin2.signif.decreased))
  
}
```

Function that performs downstream processing of output from three tools 
(if you compare between hosts) or one tool (if you compare within host).


``` r
analyse_test_output <-function(agglom.rank,
                               comparison, 
                               comparison.col_name, # the column that corresponds to comparison
                               ref.level, 
                               output.filename,
                               mean_sd_relab,
                               maaslin2.signif.decreased,
                               sample.groups = NULL ,
                               aldex2.neg.effect = NULL,
                               ancombc.signif.decreased = NULL){
  
  # if( comparison =="host"){
  # function that filters a list to remove NULL vectors or empty vectors 
  drop_empty <- function(list_x ) list_x[sapply(list_x, function(x) !is.null(x) & !length(x) == 0)]
  # create a list
  if (!is.null(aldex2.neg.effect)){
    aldex2.neg.effect <- aldex2.neg.effect %>%
      group_by(Taxon)%>%
      add_count(Taxon)%>%
      filter(n ==length(sample.groups)-1)      
  }
  
  available.decreased <- list(maaslin2 = unique(maaslin2.signif.decreased$feature),
                              aldex2 = aldex2.neg.effect$Taxon,
                              ancombc = ancombc.signif.decreased$taxon_id)
  available.decreased <- drop_empty(available.decreased)
  
  #' Find common significantly decreased features between three tools.
  if(length(available.decreased) >0){
    # Find all intersections (venn diagram style, but in a vectorised form)
    venn.decreased <- venn(available.decreased, show.plot = FALSE)
    plot(venn.decreased)
    # Get all sets (both shared and unique to each tool)
    venn.decreased.all_combinations <- attr(venn.decreased, "intersections")
    # Merge the list into a dataframe 
    # The dataframe will have two columns (values, ind), where values are 
    # differentially abundant features, and ind are tool names
    venn.decreased.df <- stack(venn.decreased.all_combinations)
    names(venn.decreased.df) <- c(agglom.rank, "Tool")
    
    if(length(available.decreased) >=2){
      message("Tools that detected significantly decreased features: ", 
              paste(names(available.decreased), collapse = ", "))
      common.decreased<-Reduce(intersect, available.decreased)
      
      # If no common decreased features, unite them all:
      if(length(common.decreased) == 0){
        message("No shared features between all tools.")
        # Check pairwise combinations: look for colon (":") in the names of intersections
        # "intersections" can also have elements not shared between groups (single group)
        venn.decreased.pairs <- venn.decreased.all_combinations[grepl(":", names(venn.decreased.all_combinations))]
        message("Pairs of tools that have shared features: ", names(venn.decreased.pairs))
        print (venn.decreased.pairs)
        
        message("Combining all features into one vector:")
        common.decreased <- unlist(available.decreased, use.names = FALSE)
        common.decreased <- unique(common.decreased)
        print(common.decreased)
      }else{
        message("Features that are shared between all tools: ")
        print(common.decreased)
      }
    }else if(length(available.decreased) == 1){
      message("Only one tool detected significantly decreased features: ", names(available.decreased))
      common.decreased<- venn.decreased.df[[agglom.rank]]
      print(common.decreased)
    }
    
  }else if (length(available.decreased) ==0){
    message("No differentially abundant features!")
    return(NULL)
  }
  
  # If we identified differentially abundant features, save them
  if(nrow(venn.decreased.df) >=1){
    # Merge the list into a dataframe and filter features
    # The dataframe will have two columns (values, ind), where values are 
    # differentially abundant features, and ind are tool names
    # common.decreased.df <- stack(available.decreased)%>%
    #   filter(values %in% common.decreased )
    # names(common.decreased.df) <- c(agglom.rank, "Tool")
    # tool names for saving the df
    venn.decreased.df.tools <- unique (venn.decreased.df$Tool)
    
    venn.decreased.df<- mean_sd_relab%>% 
      filter(get(comparison.col_name) == ref.level, 
             get(agglom.rank)%in%venn.decreased.df[[agglom.rank]])%>%
      right_join(venn.decreased.df, by = agglom.rank)%>%
      relocate(Tool, .after = agglom.rank)%>%
      arrange(Tool, -MeanRelativeAbundance)%>%
      distinct()
    # Save the common significantly decreased features with their relative abundance
    message("Showing a few differentially abundant taxa")
    head(venn.decreased.df)%>%
      knitr::kable(format = "simple")%>%
      print()
    venn.decreased.df.across_tools.fname <- file.path(diffabund.tables,
                                                      paste("common.decreased",
                                                            paste(names(available.decreased),collapse = "-"),
                                                            output.filename,
                                                            "signif.tsv",sep="-"))
    if(file.exists(venn.decreased.df.across_tools.fname)){
      message(venn.decreased.df.across_tools.fname, " exists. Skipping")
    }else{
      write.table(venn.decreased.df,
                  file = venn.decreased.df.across_tools.fname,
                  row.names = F,sep = "\t")
      
    }
    return( list(venn.decreased = venn.decreased, # vectorised form that can be plotted
                 venn.decreased.df = venn.decreased.df)) # dataframe form with taxonomy
    
  }
  # }
}
```

Plot differentially abundant features.


``` r
plot_differentially_abundant_taxa<- function(agglom.rank,
                                             comparison,
                                             ref.level,
                                             output.filename,
                                             metadata,
                                             
                                             custom.levels,
                                             ps.rel.agg,
                                             ggplot.levels,
                                             ggplot.labels,
                                             gg.labs.name,
                                             gg.title.groups,
                                             features.to_plot,
                                             maaslin2.signif.features = NULL
){
  if(setequal(custom.levels,"NMR")){
    pretty.asv.names.df<-ps.rel.agg%>%
      ungroup()%>%
      filter(get(agglom.rank)%in%features.to_plot,class=="NMR")%>%
      distinct(get(agglom.rank),.keep_all = T)
    
    
    if("OTU" %in% colnames(ps.rel.agg)){
      pretty.asv.names.df<-pretty.asv.names.df%>%
        dplyr::select(OTU,Genus)%>%
        left_join(maaslin2.signif.features, by = join_by(OTU == "feature") )%>%
        mutate(Taxon=paste0("ASV from ","<i>",Genus,"</i> (p = ",
                            round(qval,digits=4),")"))
      pretty.asv.names<-c(pretty.asv.names.df$Taxon)
      names(pretty.asv.names)<-pretty.asv.names.df$OTU
      
    }else{
      pretty.asv.names.df<-pretty.asv.names.df %>%
        dplyr::select(Species,Genus)%>%
        left_join(maaslin2.signif.features, by = join_by(Species == "feature"))%>%
        mutate(Taxon =gsub("_"," ",Species),
               Taxon=paste0("<i>",Taxon,"</i>"," (p = ",round(qval,digits = 3),")"))
      
      pretty.asv.names<-c(pretty.asv.names.df$Taxon)
      names(pretty.asv.names)<-pretty.asv.names.df$Species
      
    }
    
    
    if(comparison =="age"){
      feature.plot<-ps.rel.agg%>%
        filter(get(agglom.rank)%in%features.to_plot,class=="NMR")%>%
        # left_join(metadata, by = "Sample")%>% # Maybe remove?
        group_by_at(c(agglom.rank,"agegroup"))%>%
        ggplot(aes(x=factor(agegroup,levels=names(ggplot.labels)),
                   y=Abundance,
                   fill=factor(agegroup)))
    }else if(comparison =="sex"){
      feature.plot<-ps.rel.agg%>%
        filter(get(agglom.rank)%in%features.to_plot,class=="NMR")%>%
        group_by_at(c(agglom.rank,"sex"))%>%
        ggplot(aes(x=factor(sex,levels=names(ggplot.labels)),
                   y=Abundance,
                   fill=factor(sex)))
    }
    feature.plot <- feature.plot +
      geom_boxplot(show.legend = FALSE)+
      facet_wrap(~Species,scales = "free_y",
                 ncol = 2,
                 labeller = as_labeller(pretty.asv.names))+
      theme_bw()+
      coord_flip()+
      labs(x="",
           y="Relative abundance (%)")+
      scale_x_discrete(labels=ggplot.labels,
                       limits=ggplot.levels)+ # rename boxplot labels (x axis)
      # scale_fill_manual(values = custom.fill)+
      scale_fill_viridis_d(direction = (-1),
                           end = 0.9,
                           alpha = 0.5)+
      # ggtitle(paste0("Relative abundance of differentially abundant ASVs \nin different naked mole-rat groups"))+
      diffabund.plot.theme
    
    if(length(features.to_plot < 10 )){ # increase the size if we have few features
      feature.plot <- feature.plot +
        theme(axis.text.x = element_text(size = 12),
              axis.text.y = ggtext::element_markdown(size = 12),
              axis.title.x = element_text(size = 12),
              strip.text.x = ggtext::element_markdown(size = 12))
    }
    
    for(image.format in image.formats){
      ggsave(paste0(paste(ref.level,"specific-bacteria",output.filename,
                          sep = "-"),".",image.format),
             plot = feature.plot,
             path = diffabund.figures,
             width = 11, height = 8,units="in",
             dpi = 300,device = image.format)
    }
  }else{
    feature.plot<-ps.rel.agg%>%
      filter(get(agglom.rank)%in%features.to_plot)%>%
      group_by_at(c("class",agglom.rank))%>%
      ggplot(aes(x=factor(class,levels=rev(ggplot.levels)),
                 y=Abundance,
                 fill=factor(class,levels=rev(ggplot.levels))))+
      geom_boxplot(show.legend = FALSE)+
      facet_wrap(~Genus,scales = "free_x",
                 ncol = 4)+
      theme_bw(base_size = 5)+
      coord_flip()+
      labs(x="",
           y="Relative abundance (%)")+
      scale_color_viridis_d(breaks = rev(unname(ggplot.labels)),
                            labels=rev(unname(ggplot.labels)))+
      scale_x_discrete(labels=rev(ggplot.labels),
                       limits=rev(ggplot.levels))+ # rename boxplot labels (x axis)
      scale_fill_viridis_d(option="C")+
      # ggtitle(paste0("Relative abundance of naked mole-rat-specific taxa"))+
      diffabund.plot.theme
    
    
    for(image.format in image.formats){
      ggsave(paste0(paste(output.filename,"specific-bacteria",
                          sep = "-"),".",image.format),
             plot=feature.plot,
             path = diffabund.figures,
             width=8, height=11,units="in",
             dpi=300,device = image.format)
    }
    
  }
  print(feature.plot)
  return(feature.plot)
}
```



## Compare NMR age groups (ASV level).
Choose what to compare:


``` r
comparison<-"age"
```

Choose the reference level:


``` r
ref.level<-"agegroup0_10" 
agglom.rank<-"Species"

nmr.relations<-read.table(file.path(metadatadir,"nmr-relations.tsv"),
                          header = T,
                          sep = "\t")
```

Prepare the data for testing.


``` r
custom.md<-readRDS(file.path(metadatadir,"custom.md.ages.rds"))%>%
  filter(sequencing_type == "Naked mole-rat whole metagenome sequencing")

mean_sd_relab.nmr_ages.asv.fname <-
  file.path(community.composition.tables,
            "mean-max-min-sd-relab-nmr_ages-Species.tsv")
mean_sd_relab.nmr_ages.asv <- read.table(mean_sd_relab.nmr_ages.asv.fname,
                                         header = TRUE,
                                         sep = "\t")

nmr.age.data.for_test<-prepare_data(ps = ps.q,
                                    agglom.rank = agglom.rank,
                                    comparison = comparison,
                                    ref.level = ref.level,
                                    metadata = custom.md,
                                    relations = nmr.relations)
```

```
## Preparing data for differential microbial abundance tests on age variable at Species level. Reference: agegroup0_10
```

```
## Warning: Using an external vector in selections was deprecated in tidyselect 1.1.0.
## ℹ Please use `all_of()` or `any_of()` instead.
##   # Was:
##   data %>% select(agglom.rank)
## 
##   # Now:
##   data %>% select(all_of(agglom.rank))
## 
## See <https://tidyselect.r-lib.org/reference/faq-external-vector.html>.
## This warning is displayed once per session.
## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
## generated.
```

Run MaAsLin2 test.


``` r
maaslin.fit_data.age<-
  perform_maaslin2_test(ps.agg.wide = nmr.age.data.for_test$ps.agg.wide,
                        output.filename = nmr.age.data.for_test$output.filename,
                        comparison = comparison,
                        metadata = nmr.age.data.for_test$metadata,
                        ref.level = ref.level,
                        maaslin2.random_effects = nmr.age.data.for_test$maaslin2.random_effects,
                        maaslin2.random_effects.file_suffix = nmr.age.data.for_test$maaslin2.random_effects.file_suffix,
                        maaslin2.fixed_effects = nmr.age.data.for_test$maaslin2.fixed_effects,
                        maaslin2.reference = nmr.age.data.for_test$maaslin2.reference)
```

```
## \\wsl.localhost/Ubuntu-22.04/home/rakhimov/projects/metagenome/analysis/3-maaslin2/kraken2-output/unrare/NMR-Species-age-ref-agegroup0_10-with-relations exists. Skipping
```

``` r
maaslin.processed_output.age <- 
  process_maaslin2_output(maaslin.fit_data = maaslin.fit_data.age,
                          ps.agg = nmr.age.data.for_test$ps.agg,
                          sample.groups = nmr.age.data.for_test$sample.groups,
                          agglom.rank = agglom.rank,
                          output.filename = nmr.age.data.for_test$output.filename)
```

```
## 
## TRUE 
##    1
```

```
## \\wsl.localhost/Ubuntu-22.04/home/rakhimov/projects/metagenome/analysis/1-kraken2/20240409_17_32_40/02-differential_abundance/tables/maaslin2-NMR-Species-age-ref-agegroup0_10-signif.tsv exists. Skipping
```

```
## \\wsl.localhost/Ubuntu-22.04/home/rakhimov/projects/metagenome/analysis/1-kraken2/20240409_17_32_40/02-differential_abundance/tables/maaslin2.signif.decreased-NMR-Species-age-ref-agegroup0_10-signif.tsv exists. Skipping
```

### Downstream analysis and plotting of differentially abundant features. 


``` r
diffabund.test.summary.age <- analyse_test_output(agglom.rank = agglom.rank,
                                                  comparison = comparison,
                                                  comparison.col_name = nmr.age.data.for_test$maaslin2.fixed_effects,
                                                  ref.level = ref.level,
                                                  output.filename = nmr.age.data.for_test$output.filename,
                                                  mean_sd_relab = mean_sd_relab.nmr_ages.asv,
                                                  maaslin2.signif.decreased = maaslin.processed_output.age$maaslin2.signif.decreased)
```

<img src="005-kraken2-maaslin2_files/figure-html/unnamed-chunk-18-1.png" alt="" width="1056" />

```
## Only one tool detected significantly decreased features: maaslin2
```

```
## [1] "Mageeibacillus_indolicus"
```

```
## Showing a few differentially abundant taxa
```

```
## 
## 
## agegroup       Kingdom    Phylum      Family             Genus            Species                    Tool       OTU                                                                                                                          MeanRelativeAbundance   sdRelativeAbundance   minRelativeAbundance   maxRelativeAbundance   n_samples
## -------------  ---------  ----------  -----------------  ---------------  -------------------------  ---------  --------------------------------------------------------------------------------------------------------------------------  ----------------------  --------------------  ---------------------  ---------------------  ----------
## agegroup0_10   Bacteria   Bacillota   Oscillospiraceae   Mageeibacillus   Mageeibacillus_indolicus   maaslin2   k__Bacteria|p__Bacillota|c__Clostridia|o__Eubacteriales|f__Oscillospiraceae|g__Mageeibacillus|s__Mageeibacillus_indolicus                 0.022118             0.0082389              0.0118936              0.0319761           4
```

```
## \\wsl.localhost/Ubuntu-22.04/home/rakhimov/projects/metagenome/analysis/1-kraken2/20240409_17_32_40/02-differential_abundance/tables/common.decreased-maaslin2-NMR-Species-age-ref-agegroup0_10-signif.tsv exists. Skipping
```

``` r
plot_taxa.age <- 
  plot_differentially_abundant_taxa(agglom.rank = agglom.rank,
                                    comparison = comparison,
                                    ref.level = ref.level,
                                    output.filename = nmr.age.data.for_test$output.filename,
                                    metadata = nmr.age.data.for_test$metadata,
                                    custom.levels = nmr.age.data.for_test$custom.levels,
                                    ps.rel.agg = nmr.age.data.for_test$ps.rel.agg,
                                    ggplot.levels = nmr.age.data.for_test$ggplot.levels,
                                    ggplot.labels = nmr.age.data.for_test$ggplot.labels,
                                    gg.labs.name = nmr.age.data.for_test$gg.labs.name,
                                    gg.title.groups = nmr.age.data.for_test$gg.title.groups,
                                    features.to_plot = unique(diffabund.test.summary.age$venn.decreased.df[[agglom.rank]]),
                                    maaslin2.signif.features = maaslin.processed_output.age$maaslin2.signif.features
  )
```

<img src="005-kraken2-maaslin2_files/figure-html/unnamed-chunk-18-2.png" alt="" width="1056" />

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
##  [1] gplots_3.3.0    vegan_2.6-4     lattice_0.22-6  permute_0.9-8  
##  [5] Maaslin2_1.20.0 phyloseq_1.50.0 lubridate_1.9.4 forcats_1.0.0  
##  [9] stringr_1.5.1   dplyr_1.1.4     purrr_1.0.4     readr_2.1.5    
## [13] tidyr_1.3.1     tibble_3.2.1    ggplot2_4.0.0   tidyverse_2.0.0
## 
## loaded via a namespace (and not attached):
##   [1] Rdpack_2.6.4            DBI_1.2.3               bitops_1.0-9           
##   [4] rlang_1.1.5             magrittr_2.0.3          ade4_1.7-23            
##   [7] otel_0.2.0              compiler_4.4.3          mgcv_1.9-1             
##  [10] systemfonts_1.2.3       vctrs_0.6.5             reshape2_1.4.4         
##  [13] pkgconfig_2.0.3         crayon_1.5.3            fastmap_1.2.0          
##  [16] XVector_0.46.0          labeling_0.4.3          caTools_1.18.3         
##  [19] rmarkdown_2.30          markdown_2.0            tzdb_0.5.0             
##  [22] nloptr_2.2.1            UCSC.utils_1.2.0        ragg_1.5.0             
##  [25] xfun_0.56               zlibbioc_1.52.0         cachem_1.1.0           
##  [28] litedown_0.9            GenomeInfoDb_1.42.3     jsonlite_2.0.0         
##  [31] biomformat_1.34.0       rhdf5filters_1.18.1     Rhdf5lib_1.28.0        
##  [34] parallel_4.4.3          cluster_2.1.8.1         R6_2.6.1               
##  [37] biglm_0.9-3             bslib_0.10.0            stringi_1.8.4          
##  [40] RColorBrewer_1.1-3      boot_1.3-31             jquerylib_0.1.4        
##  [43] Rcpp_1.0.14             bookdown_0.46           iterators_1.0.14       
##  [46] knitr_1.51              IRanges_2.40.1          Matrix_1.7-4           
##  [49] splines_4.4.3           igraph_2.1.4            timechange_0.3.0       
##  [52] tidyselect_1.2.1        rstudioapi_0.18.0       yaml_2.3.12            
##  [55] ggtext_0.1.2            codetools_0.2-20        plyr_1.8.9             
##  [58] Biobase_2.66.0          withr_3.0.2             S7_0.2.0               
##  [61] evaluate_1.0.5          survival_3.8-3          getopt_1.20.4          
##  [64] xml2_1.3.8              Biostrings_2.74.1       pillar_1.11.1          
##  [67] KernSmooth_2.23-26      foreach_1.5.2           stats4_4.4.3           
##  [70] reformulas_0.4.4        pcaPP_2.0-5             generics_0.1.4         
##  [73] rprojroot_2.1.1         S4Vectors_0.44.0        hms_1.1.4              
##  [76] commonmark_2.0.0        scales_1.4.0            minqa_1.2.8            
##  [79] gtools_3.9.5            glue_1.8.0              tools_4.4.3            
##  [82] robustbase_0.99-6       data.table_1.17.8       lme4_1.1-37            
##  [85] mvtnorm_1.3-3           rhdf5_2.50.2            grid_4.4.3             
##  [88] optparse_1.7.5          ape_5.8-1               rbibutils_2.3          
##  [91] colorspace_2.1-2        nlme_3.1-167            GenomeInfoDbData_1.2.13
##  [94] cli_3.6.4               textshaping_1.0.3       viridisLite_0.4.3      
##  [97] gtable_0.3.6            DEoptimR_1.1-4          sass_0.4.10            
## [100] digest_0.6.37           BiocGenerics_0.52.0     farver_2.1.2           
## [103] htmltools_0.5.8.1       multtest_2.62.0         lifecycle_1.0.5        
## [106] httr_1.4.7              here_1.0.2              gridtext_0.1.5         
## [109] MASS_7.3-65
```

``` r
rm(list =setdiff(ls(all.names = TRUE), c("markdown.dir")))
gc()
```

```
##            used  (Mb) gc trigger  (Mb) max used  (Mb)
## Ncells  5847689 312.4   17704840 945.6 13067375 697.9
## Vcells 10215232  78.0   17173757 131.1 17173438 131.1
```

