# Creating barplots ####

# After importing the kraken2 output files into R and processing them with phyloseq,
# it's time to explore the taxonomic composition of our data.
# We will use the Polychrome package to create a custom palette for the 
# barplots.
# install.packages(c("tidyverse","ggtext","Polychrome"))
library(tidyverse)
library(Polychrome)
library(ggtext)
# Import function that creates barplot legend
source("../amplicon_nmr/code/r-scripts/create_barplot_legend.R")
## Specifying parameters and directory/file names #### 
agglom.rank<-"Species" # this is the taxonomic rank that was used for agglomeration
barplot.directory<-"./images/barplots/" # set the path where barplots will
# be saved
image.formats<-c("png","tiff")
pipeline.name<-"kraken2"
pipeline.name<-"singlem"
rdafiles.directory<-"./output/rdafiles"
rtables.directory<-"./output/rtables"
metadata.directory<-"../amplicon_nmr/output/rdafiles" # path for custom.md.ages.rds

# Import abundance table from 001-phyloseq-qiime2.R
if(pipeline.name=="kraken2"){
  ps.q.agg.date_time<-"20241003_13_52_43"
}else if (pipeline.name=="singlem"){
  ps.q.agg.date_time<-"20240929_23_33_37"
}
ps.q.agg<-readRDS(file=file.path(
  rdafiles.directory,
  paste(ps.q.agg.date_time,"phyloseq",pipeline.name,agglom.rank,
        "table.rds",sep = "-")))
# Import metadata
custom.md<-readRDS(file.path(metadata.directory,"custom.md.ages.rds"))

# Pretty labels for barplot facets that correspond to animal hosts. Here,
# the left side of the vector (values) is taken from the metadata, while
# the right side (names) are the pretty labels that will be shown on the final
# barplot
pretty.level.names<-c("NMR" = "*Heterocephalus glaber*", # better labels for facets
                      "NMRwt"="Wild *Heterocephalus glaber*",
                      "DMR" = "*Fukomys Damarensis*",
                      "B6mouse" = "B6 mouse",
                      "MSMmouse" = "MSM/Ms mouse",
                      "FVBNmouse" = "FVB/N mouse",
                      "spalax" = "*Nannospalax leucodon*",
                      "pvo" = "*Pteromys volans orii*",
                      "hare" = "*Lepus europaeus*",
                      "rabbit" = "*Oryctolagus cuniculus*")

# Set custom levels for the barplot. These are the animal hosts that will be used
# in barplots.
# Use only the taxa that are present in the workspace
# (custom.md is metadata from the rdafile)
custom.levels<-intersect(names(pretty.level.names),custom.md$class)
# Filter the phyloseq object to retain animal hosts from custom.levels
ps.q.agg<-ps.q.agg%>%
  filter(class%in%custom.levels,Abundance!=0)

# Check the total number of taxa (not filtered by mean relative abundance)
ps.q.agg%>%
  ungroup()%>%
  select(matches(paste0("^",agglom.rank,"$")))%>%
  pull(.)%>%
  unique()%>%
  sort()%>%
  length()
# select(matches(paste0("^",agglom.rank,"$"))) will select variables that 
# match a pattern (regex).
# So, if we are agglomerating by Phylum, the command will select the Phylum column

# Creating a legend ####
### 1. Extract the unique taxa (corresponding to agglom.rank) from the dataset ####
# Result: taxa.list vector
taxa.list<-ps.q.agg%>%
  group_by_at(c("class",agglom.rank))%>%
  filter(MeanRelativeAbundance>=1)%>%
  ungroup()%>%
  select(matches(paste0("^",agglom.rank,"$")))%>%
  pull(.)%>%
  unique()
# select(matches(paste0("^",agglom.rank,"$"))) will select variables that 
# match a pattern (regex).
# So, if we are agglomerating by Phylum, the command will select the Phylum column

### 2. Get the taxonomic ranks for ordering ####
# The order of the legend is: "Remainder","Kingdom", "Phylum", "Class", 
# "Order", "Family". But that is when we agglomerate by Genus.
# If we agglomerate by a higher rank, we need to remove the unnecessary ranks.
# For example, if agglom.rank="Phylum", we need to match the agglom.rank
# with the custom.order vector and extract the ranks before agglom.rank
# (excluding the agglom.rank!). The resulting custom.order vector becomes shorter
taxa.list<-c(taxa.list,"Remainder")
custom_order <- c("Remainder","Kingdom", "Phylum", "Class", "Order", "Family","Genus","Species")
# If we agglomerate by higher level (Order,Class, etc), need to adjust the rank
if(agglom.rank%in%custom_order){
  agglom.rank.index<-match(agglom.rank,custom_order)
  custom_order<-custom_order[1:agglom.rank.index-1]
}

### 3. Find classified and unclassified taxa indices. ####
# The sub command removes everything before a space (e.g. Bacteria Kingdom becomes
# Bacteria) in the taxa.list vector. The match commmand matches the taxa.list
# to custom.order (taxonomic ranks before agglom.rank). This will find 
# unclassified taxa because they are from higher taxonomic rank (e.g. if we 
# agglomerate by Genus, taxa with "Family" or "Order" in their name are 
# unclassified ones). The rank vector has length equal to length of the
# taxa.list. Each value corresponds to an element from taxa.list that matches
# to the index of custom.order vector. So, if the 2nd element of taxa.list is
# "Kingdom Bacteria", the 2nd element of the rank vector is "1" because the 
# taxonomic rank "Kingdom" matches "Kingdom Bacteria" and is 1st in the custom.order
# vector.
# However, the classified taxa won't match the custom.order 
# vector (because they are classified). Thus, the classified elements will be NA.
# We substitute NA values with the length of custom.order +1. So, if we agglomerate
# by Genus, the length of custom.order is 6. But the classified taxa (which were
# not found in custom.order and thus have NA in the rank vector) are out of the 
# ranking. Their index is 6+1=7. So, in the rank vector all NAs will become 7.
rank <- match(sub(".* ", "", taxa.list),custom_order)
rank[is.na(rank)] <- length(custom_order) + 1

### 4. Extract classified taxa: their rank is the last in the custom_order vector ####
### 4.1. We match taxa.list to the rank vector and extract elements with rank values ####
# equal to length(custom.order)+1.
classified.taxa<-taxa.list[rank==length(custom_order) + 1]
### 4.2. Unclassified taxa are ones with indices found in the custom.order ####
unclassified.taxa<-taxa.list[rank!=length(custom_order) + 1]

### 4.3. In the ps.q.agg dataset, we find the index of the agglom.rank column and  ####
# preceding.rank column (e.g Genus and Family)
if(agglom.rank=="OTU"){
  agglom.rank.col<-which(colnames(ps.q.agg) =="Species")
}else{
  agglom.rank.col<-which(colnames(ps.q.agg) ==agglom.rank)
  preceding.rank.col<-which(colnames(ps.q.agg) ==agglom.rank)-1
  preceding.rank<-colnames(ps.q.agg)[preceding.rank.col]
}

### 4.4. Create a dataframe of classified taxa (agglom.rank + preceding rank) for the legend ####
# From the ps.q.agg dataset, we obtain classified taxa, select columns
# for agglom.rank and preceding.rank, then retain only unique rows
# To create the Taxon.bp column, we concatenate agglom.rank with preceding.rank
# The taxa.for_bp.df has columns for agglom.rank (e.g Genus), preceding.rank 
# (e.g. Family), and Taxon.bp (merged agglom.rank and preceding.rank)
taxa.for_bp.df<-ps.q.agg%>%
  ungroup()%>%
  filter(get(agglom.rank)%in%classified.taxa)%>%
  select(all_of(c(agglom.rank,preceding.rank,"Kingdom")))%>%
  distinct()%>%
  unite("Taxon.bp",agglom.rank:preceding.rank,sep = " (",remove = FALSE)%>%
  select(all_of(c("Kingdom",agglom.rank,preceding.rank,"Taxon.bp")))%>%
  mutate("Append"=")")%>%
  unite("Taxon.bp",Taxon.bp:Append,sep = "")

### 4.5. Order classified taxa by preceding.rank then agglom.rank ####
# These are agglom.rank (preceding.rank) format strings for barplot
taxa.for_bp.df<-taxa.for_bp.df%>%
  arrange(get(preceding.rank),get(agglom.rank))
### 4.6 We select the Taxon.bp column of the taxa.for_bp.df and store it as a  ####
# separate vector taxa.for_bp.list
taxa.for_bp.list<-taxa.for_bp.df$Taxon.bp

### 5. Now we order unclassified taxa: start by matching ####
### 5.1 Match the unclassified taxa to indices of custom.order ####
newrank <- match(sub(".* ", "", unclassified.taxa),custom_order)
# newrank[is.na(newrank)] <- length(custom_order) + 1

### 5.2 Order the unclassified taxa vector ####
unclassified.taxa<-unclassified.taxa[order(newrank)]

### 5.3 Sorting inside each rank: split the vector by taxonomic rank ####
# If agglom.rank="Genus", unclassified.taxa.split will be a list of 6 lists
# (Remainder, Kingdom, Phylum, Class, "Order", "Family").
unclassified.taxa.split <- split(unclassified.taxa, newrank[order(newrank)])
# Sort inside each rank
unclassified.taxa.sorted <- unlist(lapply(unclassified.taxa.split, sort))
# The elements in unclassified.taxa.sorted are named by their rank+position inside
# rank (e.g. the 19th taxon in Family will be 619: 6 for the Family list, 19 for the 
# index). Get rid of these names!
unclassified.taxa.sorted<-unname(unclassified.taxa.sorted)

### 6. Add sorted unclassified taxa to the sorted classified taxa (create taxa.for_bp.list) ####
taxa.for_bp.list<-c(unclassified.taxa.sorted,taxa.for_bp.list)
# Remainder is the first element. We will update the name
taxa.for_bp.list[1]<-"Remainder (Mean relative abundance < 1%)"

### 7. Add barplot taxa to the main dataframe ####
# Join the ps.q.agg dataset with classified taxa dataset (taxa.for_bp.df)
ps.q.agg.for_bp<-ps.q.agg%>%
  left_join(subset(taxa.for_bp.df,select = -Kingdom),by=agglom.rank)%>%
  ungroup()%>%
  select(-paste0(preceding.rank,".y"))%>% # remove the preceding.rank column from taxa.for_bp.df
  rename(!!preceding.rank:=paste0(preceding.rank,".x")) # rename the preceding.rank column
# !!preceding.rank:= will evaluate the variable

### 7.1 Remove NAs: unclassified taxa with MeanRelativeAbundance>=1 ####
# They are NA because they were not found in taxa.for_bp.df (which has only
# classified data)
# We remove NA by simply assigning agglom.rank to Taxon.bp. So, if Genus is 
# "Bacteria Kingdom", the Taxon.bp will also become "Bacteria Kingdom".
ps.q.agg.for_bp[which(ps.q.agg.for_bp$MeanRelativeAbundance>=1&is.na(ps.q.agg.for_bp$Taxon.bp)),"Taxon.bp"]<-
  ps.q.agg.for_bp[which(ps.q.agg.for_bp$MeanRelativeAbundance>=1&is.na(ps.q.agg.for_bp$Taxon.bp)),agglom.rank]

# Taxa with MeanRelativeAbundance<1% become "Remainder"
ps.q.agg.for_bp[which(ps.q.agg.for_bp$MeanRelativeAbundance<1),"Taxon.bp"]<-
  "Remainder (Mean relative abundance < 1%)"

# find class and agglom.rank columns, just in case
classcol<-which(colnames(ps.q.agg.for_bp) =="class")
if(agglom.rank=="OTU"){
  agglom.rank.col<-which(colnames(ps.q.agg.for_bp) =="Species")
}else{
  agglom.rank.col<-which(colnames(ps.q.agg.for_bp) ==agglom.rank)
}

# Relative abundance per kingdom (including bacteria)
if(pipeline.name=="kraken2"){
  unclassified_reads<-read.table("./output/kraken2_pipeline/20240409_17_32_40_unclassified_reads.tsv",
                                 header = T)
  colnames(unclassified_reads)[which(colnames(unclassified_reads)=="Unclassified")]<-"Abundance"
  unclassified_reads$Kingdom<-"Unclassified"
}

# Kingdom barplot ####
if(pipeline.name=="kraken2"){
  kingdom.plot.with_bact<-ps.q.agg.for_bp%>%
    group_by(Sample,Kingdom)%>%
    summarise(Abundance=sum(Abundance))%>%
    bind_rows(unclassified_reads)
}else{
  kingdom.plot.with_bact<-ps.q.agg.for_bp%>%
    mutate(Kingdom=ifelse(grepl("Kingdom|Phylum|Class|Order|Family|Genus",Species)&Kingdom=="Bacteria",
                          "Unclassified Bacteria",Kingdom))%>%
    mutate(Kingdom=ifelse(grepl("Kingdom|Phylum|Class|Order|Family|Genus",Species)&Kingdom=="Archaea",
                          "Unclassified Archaea",Kingdom))%>%
    group_by(Sample,Kingdom)%>%
    summarise(Abundance=sum(Abundance))
}
kingdom.plot.with_bact<-kingdom.plot.with_bact%>%
  group_by(Kingdom)%>%
  summarise(TotalAbundance=sum(Abundance))%>%
  mutate(TotalAbundancePerTax=TotalAbundance/sum(TotalAbundance)*100)%>%
  arrange(-TotalAbundancePerTax)%>%
  ggplot(aes(x=reorder(Kingdom,-TotalAbundancePerTax),
             y=TotalAbundancePerTax,
             fill=Kingdom))+
  geom_bar(stat = "identity")+
  ylab("Relative abundance from all samples (%)")+
  xlab("")+
  theme_bw()+
  ggtitle("Whole metagenome sequencing profile of naked mole-rat gut microbiota
         (Kingdom level)")+
  geom_text(aes(x=Kingdom, 
                y=TotalAbundancePerTax, 
                label=round(TotalAbundancePerTax,digits = 3)),
            vjust=-0.5,size=8)+ # add text with percentage above the bar
  theme(axis.line = element_blank(), 
        axis.text.x = element_text(size=20),
        axis.text.y = element_text(size=20), # size of y axis ticks
        axis.title = element_text(size = 20), # size of axis names
        plot.title = ggtext::element_markdown(size = 25), # the plot 
        # title will be recognised as a markdown object, so we can
        # add line breaks (cause host names are too long)
        plot.caption = element_text(size=23),# size of plot caption
        legend.text = element_text(size = 20),# size of legend text
        legend.title = element_text(size = 25), # size of legend title
        legend.position = "none") # legend on the right

# save the plots: kingdom.plot.with_bact
# for(image.format in image.formats){
#   ggsave(paste0(barplot.directory,
#                 paste(paste(format(Sys.time(),format="%Y%m%d"),
#                             format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                       "NMR","kingdom.plot.with_unclas",
#                       sep = "-"),".",image.format),
#          plot=kingdom.plot.with_bact,
#          width = 5000,height = 3000,
#          units = "px",dpi=300,device = image.format)
# }

# Check library size per sample
ps.q.agg.for_bp%>%
  group_by(Sample)%>%
  summarise(TotalSample=sum(Abundance))%>%
  ggplot(aes(x=Sample,y=TotalSample,fill=TotalSample))+
  geom_bar(stat = "identity")



ps.q.legend<-create_barplot_legend(tax.df = ps.q.agg.for_bp,
                                   inside.host = TRUE,
                                   is.metagenome = TRUE,
                                   split.by.comparison=FALSE,
                                   comparison.var.name = "agegroup",
                                   tax.rank = "Species",
                                   left.group.name = "agegroup0_10",
                                   right.group.name = "agegroup10_16",
                                   taxa.to.plot =taxa.for_bp.list,
                                   metadata.table=custom.md # TODO fix because it's only when we compare inside nmr
)

## Plot the barplots ####
# We need to choose colors for the taxa in our barplot. They should be 
# distinguishable, so we can't choose similar colors. Or at least we shouldn't 
# put them next to each other.
# 
# We will use `createPalette` function from the `Polychrome` package.
# The function takes the number of colors for the palette and seed colors that 
# will be used for palette generation. In this case, we use the number of
# rows in our legend (taxa) and rainbow colors (to avoid having similar colors 
#                                               next to each other). 
# We also need to set the random seed because the output
# is a bit random. The output is a vector of colors.
set.seed(1)
plot.cols<-createPalette(nrow(ps.q.legend)-1,
                         seedcolors =rainbow(7))# input: number of rows

# in our legend and the seed colors that we decide to be rainbow

# The vector of colors should be named according to our legend
# because we will use this mapping in the barplot. So, each color in the vector
# correspond to each color in the legend. Remember, remainder taxa are 
# all merged into a single entry ("Remainder"), so there's just one color for 
# remainder portion.

# Decrease alpha for unclassified
plot.cols<-c("#C1CDCD",plot.cols)

col.vec<-setNames(plot.cols,ps.q.legend$Taxon.bp)
# Create the barplot with ggplot2. First, we take the agglomerated
# dataset that we obtained in the `001-phyloseq-qiime2.R` and merge it with
# the dataset of total abundances, so we can know how many reads were in
# each sample. We will concatenate the sample name on the x-axis with the 
# number of reads in that sample. It will look like "`Sample 1 (n = 25000)`".
# Actually, this kind of string will be stored in a new column 
# called `NewSample`.
# 
# Then, we will convert the `class` column (host names) into factors that we 
# defined in `custom.levels`. This will order our bars according to the vector of
# levels. It must be factor because this allows us ordering `custom.levels` as
# we want, otherwise it would be alphabetic. And in our vector, the first level 
# is naked mole-rats. So, the first facet will also be naked mole-rats.
# Facets are panels that correspond to animal host.
# 
# The `ggplot` command needs an aesthetic: the x axis will correspond to
# the `NewSample` column (sample names with sample size), while the y-axis
# will be the column of relative abundances. We also need a `fill` value
# which is basically the vector that will be used for coloring the barplot.
# We use taxa from the `Taxon.bp` column because each section of each bar
# is a taxon that must be colored. But we must convert the `Taxon.bp` into
# factor, so it can map to the vector of color. **The order of factorised
# `Taxon.bp` is based on the `Taxon` column from the legend**.
dominant.species.bact<-c("Phascolarctobacterium_faecium",
                         "Bacteroides_xylanisolvens",
                         "Enterococcus_faecium",
                         "Anaerostipes_caccae",
                         "Bacteroides_thetaiotaomicron",
                         "Parabacteroides_distasonis")
dominant.species.nonbact<-c("Ipomoea_triloba",
                            "Cryptomeria_japonica",
                            "Macadamia_integrifolia")
pretty.level.names<-names(table(custom.md$old_agegroup))
names(pretty.level.names)<-names(table(custom.md$agegroup))
custom.levels<-names(pretty.level.names)


sample.levels<-custom.md%>%
  ungroup()%>%
  select(Sample,age)%>%
  arrange(age)%>%
  distinct()%>%
  mutate(NewSample=paste0(Sample," (",age,")"))%>%
  mutate(Sample=factor(Sample,levels=Sample),
         NewSample=factor(NewSample,levels=NewSample))

mainplot<-ps.q.agg.for_bp%>%
  left_join(sample.levels)%>%
  # filter(Species%in%dominant.species.bact)%>%
  ungroup()%>%
  ggplot(aes(x=NewSample, y=RelativeAbundance,
             fill=factor(Taxon.bp, levels=ps.q.legend$Taxon.bp)))+
  geom_bar(stat = "identity")+ # barplot
  # facet_grid(~agegroup, # separate animal hosts
  #            scales="free",  # each species will have its own bars inside
  #            # facet (instead of all bars)
  #            space = "free", # bars will have same widths
  #            labeller = labeller(agegroup=pretty.level.names))+ # labeller will
  # change facet labels to custom
  # guides(fill=guide_legend(ncol=1))+ # legend as one column
  coord_cartesian(expand=FALSE) +
  scale_fill_manual(values = col.vec,
                    labels=ps.q.legend$new.colors)+ # custom fill that is based on our
  # custom palette
  # scale_fill_manual(values = col.vec[grepl(paste(dominant.species.bact,collapse = "|"),names(col.vec))],
  #                   labels=ps.q.legend$new.colors[grepl(paste(dominant.species.bact,collapse = "|"),ps.q.legend$new.colors)])+ # custom fill that is based on our
  # custom palette
  xlab("") +
  ylab("Relative Abundance (%)")+
  labs(fill="Taxon")+
  theme_bw()+
  guides(fill=guide_legend(ncol=1))+ # legend as one column
  ggtitle(paste(agglom.rank,"level gut microbiota profiles of fecal samples from naked mole-rats"))+
  theme(plot.margin=unit(c(1,1,1,1.5), 'cm'),
        axis.line = element_blank(), #TODO: what does it do?
        strip.text.x = ggtext::element_markdown(size = 20),# the name of
        # each facet will be recognised as a markdown object, so we can
        # add line breaks (cause host names are too long)
        panel.spacing = unit(0.8, "cm"), # increase distance between facets
        axis.text.x = element_text(angle=45,size=20,hjust=1),# rotate
        # the x-axis labels by 45 degrees and shift to the right
        axis.text.y = element_text(size=20), # size of y axis ticks
        axis.title = element_text(size = 20), # size of axis names
        plot.title = element_text(size = 25), # size of plot title
        plot.caption = element_text(size=23), # size of plot caption
        legend.text = element_markdown(size = 20), # size of legend text
        legend.title = element_text(size = 25), # size of legend title
        legend.position = "right") # legend under the plot

for(image.format in image.formats){
  ggsave(paste0(barplot.directory,
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                      pipeline.name,"barplot","NMR",
                      agglom.rank,
                      sep = "-"),".",image.format),
         plot=mainplot,
         width = 6000,height = 3200,
         units = "px",dpi=300,device = image.format)
}



akkermansia.plot<-ps.q.agg.for_bp%>%
  left_join(custom.md,by="Sample")%>%
  group_by(Sample)%>%
  filter(Species=="Akkermansia_muciniphila")%>%
  ggplot(aes(x=Sample,y=RelativeAbundance))+
  ylab("Relative Abundance (%)")+
  geom_bar(stat="identity")+
  ggtitle("Akkermansia muciniphila relative abundance")+
  facet_grid(~agegroup,scales = "free_x")+
  theme_bw()+
  theme(axis.line = element_blank(), 
        strip.text.x = ggtext::element_markdown(size = 20),# the name of 
        # each facet will be recognised as a markdown object, so we can
        # add line breaks (cause host names are too long)
        axis.text.x = element_text(angle=45,size=20,hjust=1),# rotate 
        # the x-axis labels by 45 degrees and shift to the right
        axis.text.y = element_text(size=20), # size of y axis ticks
        axis.title = element_text(size = 20), # size of axis names
        plot.title = ggtext::element_markdown(size = 25), # the plot 
        # title will be recognised as a markdown object, so we can
        # add line breaks (cause host names are too long)
        plot.caption = element_text(size=23),# size of plot caption
        legend.text = element_text(size = 20),# size of legend text
        legend.title = element_text(size = 25), # size of legend title
        legend.position = "right")
ggsave(paste0("./images/barplots/",
              paste(Sys.Date(),"barplot","NMR-akkermansia.plot",
                    agglom.rank,sep = "-"),".png"),
       plot=akkermansia.plot,
       width = 6000,height = 3000,
       units = "px",dpi=300,device = "png")
ggsave(paste0("./images/barplots/",
              paste(Sys.Date(),"barplot","NMR-akkermansia.plot",
                    agglom.rank,sep = "-"),".tiff"),
       plot=akkermansia.plot,
       width = 6000,height = 3000,
       units = "px",dpi=300,device = "tiff")


campylobacter_jejuni.plot<-ps.q.agg.for_bp%>%
  left_join(custom.md,by="Sample")%>%
  group_by(Sample)%>%
  filter(Species=="Campylobacter_jejuni")%>%
  ggplot(aes(x=Sample,y=RelativeAbundance))+
  ylab("Relative Abundance (%)")+
  geom_bar(stat="identity")+
  ggtitle("Campylobacter jejuni relative abundance")+
  facet_grid(~agegroup,scales = "free_x")+
  theme_bw()+
  theme(axis.line = element_blank(), 
        strip.text.x = ggtext::element_markdown(size = 20),# the name of 
        # each facet will be recognised as a markdown object, so we can
        # add line breaks (cause host names are too long)
        axis.text.x = element_text(angle=45,size=20,hjust=1),# rotate 
        # the x-axis labels by 45 degrees and shift to the right
        axis.text.y = element_text(size=20), # size of y axis ticks
        axis.title = element_text(size = 20), # size of axis names
        plot.title = ggtext::element_markdown(size = 25), # the plot 
        # title will be recognised as a markdown object, so we can
        # add line breaks (cause host names are too long)
        plot.caption = element_text(size=23),# size of plot caption
        legend.text = element_text(size = 20),# size of legend text
        legend.title = element_text(size = 25), # size of legend title
        legend.position = "right")
ggsave(paste0("./images/barplots/",
              paste(Sys.Date(),"barplot","NMR-campylobacter_jejuni.plot",
                    agglom.rank,sep = "-"),".png"),
       plot=campylobacter_jejuni.plot,
       width = 6000,height = 3000,
       units = "px",dpi=300,device = "png")
ggsave(paste0("./images/barplots/",
              paste(Sys.Date(),"barplot","NMR-campylobacter_jejuni.plot",
                    agglom.rank,sep = "-"),".tiff"),
       plot=campylobacter_jejuni.plot,
       width = 6000,height = 3000,
       units = "px",dpi=300,device = "tiff")


taxa.comparison.plot<-ps.q.agg.for_bp%>%
  left_join(custom.md,by="Sample")%>%
  group_by(Sample)%>%
  filter(Species=="Campylobacter_jejuni"| Species=="Akkermansia_muciniphila")%>%
  ggplot(aes(x=Sample,y=RelativeAbundance))+
  ylab("Relative Abundance (%)")+
  geom_bar(stat="identity")+
  ggtitle("Akkermansia muciniphila vs Campylobacter jejuni relative abundance")+
  facet_grid(~Species,scales = "free_x")+
  theme_bw()+
  theme(axis.line = element_blank(), 
        strip.text.x = ggtext::element_markdown(size = 20),# the name of 
        # each facet will be recognised as a markdown object, so we can
        # add line breaks (cause host names are too long)
        axis.text.x = element_text(angle=45,size=20,hjust=1),# rotate 
        # the x-axis labels by 45 degrees and shift to the right
        axis.text.y = element_text(size=20), # size of y axis ticks
        axis.title = element_text(size = 20), # size of axis names
        plot.title = ggtext::element_markdown(size = 25), # the plot 
        # title will be recognised as a markdown object, so we can
        # add line breaks (cause host names are too long)
        plot.caption = element_text(size=23),# size of plot caption
        legend.text = element_text(size = 20),# size of legend text
        legend.title = element_text(size = 25), # size of legend title
        legend.position = "right")
ggsave(paste0("./images/barplots/",
              paste(Sys.Date(),"barplot","NMR-taxa.comparison.plot",
                    agglom.rank,sep = "-"),".png"),
       plot=taxa.comparison.plot,
       width = 6000,height = 3000,
       units = "px",dpi=300,device = "png")
ggsave(paste0("./images/barplots/",
              paste(Sys.Date(),"barplot","NMR-taxa.comparison.plot",
                    agglom.rank,sep = "-"),".tiff"),
       plot=taxa.comparison.plot,
       width = 6000,height = 3000,
       units = "px",dpi=300,device = "tiff")

ps.q.agg.for_bp%>%
  left_join(custom.md,by="Sample")%>%
  group_by(Sample)%>%
  filter(Species=="Campylobacter_jejuni"| Species=="Akkermansia_muciniphila")%>%
  ggplot(aes(x=Sample,y=RelativeAbundance))+
  ylab("Relative Abundance (%)")+
  geom_bar(stat="identity")+
  ggtitle("Akkermansia muciniphila vs Campylobacter jejuni relative abundance")+
  facet_grid(~Species,scales = "free_x")+
  theme_bw()
