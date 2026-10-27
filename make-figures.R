#' ---
#' output: 
#'   bookdown::html_document2:
#'      toc: true
#' params:
#'   active.analysis: ''
#' ---
#' 
#' ```{r, setup make-figures.R, include=FALSE}
#' knitr::opts_knit$set(root.dir = '/home/rakhimov/projects/metagenome')
#' ```
#' ```{r, echo = FALSE}
#' # For showing images, tables, etc: Use global path
#' #knitr::spin("make-figures.R", knit = FALSE)
#' #file.rename("make-figures.Rmd", "markdown/make-figures.Rmd")
#' #rmarkdown::render('./markdown/make-figures.Rmd', 'html_document',
#' # knit_root_dir="/home/rakhimov/projects/metagenome/")
#' ```
#' 
#+ echo=FALSE
# Generating figures ####
#' # Generating figures
#' 
#+ echo=FALSE
## Introduction ####
#'
#' ## Introduction
#' This script generates all figures from WMS data:
#' 
#' 1) Taxonomic classification from Kraken2 by sample and by kingdom (combined)
#' 
#' 2) Taxonomic classification from GTDB-tk by sample
#' 
#' 3) Taxonomic classification from BLASTN by sample
#' 
#+ echo=FALSE
## 1. Load necessary libraries and scripts. ####
#'
#' ## Load necessary libraries and scripts.
library(tidyverse)
library(Polychrome)
library(ggtext)
library(cowplot)

source(here::here("config/R/config.R"))# config file with global variables
source(here::here("config/R/themes.R"))# config file with themes

#+ echo=FALSE
## 2. Specifying parameters and directory/file names. #### 
#'
#' ## Specifying parameters and directory/file names. 
agglom.rank<-"Species" # this is the taxonomic rank that was used for agglomeration

#+ echo=FALSE
## 3. Plots. ####
#'
#' ## Plots.

#+ echo=FALSE
### Barplot from Kraken2 with kingdom plot
#' ### Barplot from Kraken2 with kingdom plot
kingdom.plot.with_bact <- readRDS(file = file.path(kraken2.manuscript.figures,
                                                   "kraken2.kingdom.plot.with_bact.rds"))

kingdom.plot.with_bact <- kingdom.plot.with_bact +
  coord_cartesian(expand = c("bottom" = F,
                             "left" = F,
                             "right" = F,
                             "top" = T),
                  ylim = c(0,80))+
  geom_text(aes(x=Kingdom, 
                y=TotalAbundancePerTax, 
                label=round(TotalAbundancePerTax,digits = 3)),
            vjust=-0.5,size = 2)+ # add text with percentage above the bar
  theme(plot.margin=unit(c(0.2,0,0.1,0.1), 'cm'),
        axis.line = element_blank(), 
        axis.text.x = element_text(size=5),
        axis.text.y = element_text(size=5), # size of y axis ticks
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 7), # size of axis names
        plot.title = ggtext::element_markdown(size = 25), # the plot 
        # title will be recognised as a markdown object, so we can
        # add line breaks (cause host names are too long)
        plot.caption = element_text(size=23),# size of plot caption
        legend.text = element_text(size = 10),# size of legend text
        legend.title = element_text(size = 12), # size of legend title
        legend.position = "none",
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) # no legend 

#' Kraken2 taxonomy:
mainplot <- readRDS(file = file.path(kraken2.manuscript.figures,
                         "kraken2.mainplot.rds"))

mainplot <- mainplot +
  theme(plot.margin=unit(c(0.2,0,0.1,0.1), 'cm'),
        strip.text.x = ggtext::element_markdown(size = 10),# the name of
        # each facet will be recognised as a markdown object, so we can
        # add line breaks (cause host names are too long)
        panel.spacing = unit(0.8, "cm"), # increase distance between facets
        axis.text.x = element_text(angle=45,size=9,hjust=1),# rotate
        # the x-axis labels by 45 degrees and shift to the right
        axis.text.y = element_text(size=10), # size of y axis ticks
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12), # size of axis names
        plot.title = element_text(size = 12), # size of plot title
        plot.caption = element_text(size=12), # size of plot caption
        legend.box.spacing = unit(0.01, units="cm"),
        legend.key.size = unit(0.5, units= "cm"),
        legend.text = element_markdown(size = 9), # size of legend text
        legend.title = element_text(size = 9), # size of legend title
        legend.position = "right") # legend on the right

#+ echo=FALSE
### Basic taxonomy barplot (GTDB-tk) (**Figure 5**) ####
#' ### Basic taxonomy barplot (GTDB-tk) (**Figure 5**) ####
gtdbtk.coverm.drep.plot<- readRDS( file.path(mag.manuscript.figures, 
                                             "gtdbtk.coverm.drep.plot.rds"))

gtdbtk.coverm.drep.plot <- gtdbtk.coverm.drep.plot +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size=7),
        axis.text.y = element_text(size=10), # size of y axis ticks
        axis.title = element_text(size = 11), # size of axis names
        plot.caption = element_text(size=10), # size of plot caption
        legend.title = element_text(size = 12) # size of legend title
  )+
  theme(plot.title = element_text(size=10), 
        panel.spacing = unit(0.8, "cm"), # increase distance between facets
        legend.text = element_markdown(size = 10), # size of legend text
        legend.position = "right",
        legend.key.size = unit(0.4, units= "cm"),
        legend.box.spacing = unit(0.01, units="cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank())

#+ echo=FALSE
### Taxonomy barplot BLASTN (**Figure 5b**) ####
#' ### Taxonomy barplot BLASTN (**Figure 5b**) ####
blastn.coverm.megahit.plot <- readRDS(file.path(mag.manuscript.figures, 
                                                "blastn.coverm.megahit.plot.rds"))

blastn.coverm.megahit.plot<-blastn.coverm.megahit.plot +
  scale_fill_manual(values = blastn.col.vec)+ # custom fill that is based on our 
  theme(axis.text.x = element_text(size = 7),
        axis.text.y = element_text(size=10), # size of y axis ticks
        axis.title = element_text(size = 11), # size of axis names
        plot.title = element_text(size = 12), # size of plot title
        plot.caption = element_text(size=10), # size of plot caption
        legend.text = element_markdown(size = 10), # size of legend text
        legend.title = element_text(size = 12), # size of legend title
        legend.box.spacing = unit(0.01, units="cm"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        legend.position = "right") # legend on the right

#+ echo=FALSE
### The barplot of all genes for each annotation (4 bars) ####
#' ### The barplot of all genes for each annotation (4 bars) ####
uniq.gene.counts.bp <- readRDS(file = file.path(mag.manuscript.figures, 
                                                "uniq.gene.counts.bp.rds"))

uniq.gene.counts.all.bp <- uniq.gene.counts.bp +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size=20), # size of y axis ticks
        axis.title = element_text(size = 20), # size of axis names
        plot.title = element_text(size = 25), # size of plot title
        plot.caption = element_text(size=23), # size of plot caption
        legend.text = element_text(size = 20), # size of legend text
        legend.title = element_text(size = 25), # size of legend title
        legend.position = "inside", # legend inside the plot
        legend.position.inside=c(0.15,0.8),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) # legend on the left top


# Import metadata
custom.md<-readRDS(file.path(metadata.directory,"custom.md.ages.rds"))%>%
  filter(sequencing_type == "Naked mole-rat whole metagenome sequencing")

pretty.level.names<-c("NMR" = "*Heterocephalus glaber*") # better labels for facets)

custom.levels<-intersect(names(pretty.level.names),custom.md$class)
# Filter the phyloseq object to retain animal hosts from custom.levels
ps.q.agg<-ps.q.agg%>%
  filter(class%in%custom.levels,Abundance!=0)


#+ echo=FALSE
## Save figures ####
#' ## Save figures ####
figure4 <- ggdraw(mainplot)+
  draw_plot(kingdom.plot.with_bact, 
            x = 0.02,
            y = 0.62,
            width = 0.85,
            height = 0.5,
            scale=0.4,
            hjust = 0.17)#+
# draw_plot_label(
#   # c("b", "a"),
#   c(0.02, 0.12),
#   c(1, 0.98),
#   size = 20
# )

# save the plots: kingdom.plot.with_bact
for(image.format in image.formats){
  ggsave(filename = paste("figure4",image.format,sep = "."),
         path = kraken2.manuscript.figures,
         plot = figure4,
         device = image.format,
         dpi=300,
         width = 170,
         height = 130, #190
         units = "mm",
         scale = 1)
}

figure5<-plot_grid(gtdbtk.coverm.drep.plot,
                   blastn.coverm.megahit.plot, 
                   # uniq.gene.counts.all.bp,
                   labels = "auto",
                   label_size = 20,
                   ncol = 1)
# save the plots: kingdom.plot.with_bact
for(image.format in image.formats){
  ggsave(filename = paste("figure5",image.format,sep = "."),
         path = mag.manuscript.figures,
         plot = figure5,
         device = image.format,
         dpi = 300,
         width = 170,
         height = 180, # 190,
         units = "mm",
         scale = 1# 2.2
         )
}
sessionInfo()
rm(list = ls(all=TRUE))
gc()
