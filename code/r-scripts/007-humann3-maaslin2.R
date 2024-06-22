library(tidyverse)
library(Maaslin2)
merged.date_time<-"20240618_11_21_44"
metadatadir<-paste0("../amplicon_nmr/data/metadata/pooled-metadata/") # directory with metadata

comparison<-"age"
host<-"NMR"
ref.level<-"agegroup0_10" # choose the reference level

# Pathways ####
# pathabundance.filename<-file.path("output/humann3_pipeline",
#           paste("humann_out",merged.date_time,"merged",sep = "_"),
#           paste(merged.date_time,"humann_out_pathabundance-cpm-filtered-total-renorm.tsv",sep = "_"))
pathabundance.filename<-file.path("output/humann3_pipeline",
                                  paste("humann_out",merged.date_time,"merged",sep = "_"),
                                  paste(merged.date_time,"humann_out_pathabundance-cpm-total-filtered-renorm.tsv",sep = "_"))

pathabundance<-read.table(pathabundance.filename, 
                            header = T, 
                            sep = "\t",
                            fill=TRUE,
                            quote = "", 
                            comment.char = "@",
                            na.strings = c("","NA"))
# pathab.long<-pathabundance%>%pivot_longer(names_to = "Sample",cols = 2:12)
# pathab.long%>%
#   filter(grepl("ARGSYNBSUB-PWY: L-arginine biosynthesis II (acetyl cycle)",X..Pathway,fixed = T))%>%
#   group_by(Sample)%>%
#   summarise(sumab=sum(value))

dim(pathabundance)
head(pathabundance)
# pathabundance[is.na(pathabundance)]<-"NA"
colnames(pathabundance)<-gsub("X..Pathway","Pathway",colnames(pathabundance))
colnames(pathabundance)<-gsub("X2D10","2D10",colnames(pathabundance))
colnames(pathabundance)<-gsub("X2D14","2D14",colnames(pathabundance))
colnames(pathabundance)<-
  gsub("_trim_merged_Abundance","",colnames(pathabundance))
colnames(pathabundance)<-
  gsub(".CPM","",colnames(pathabundance))

# Check how many pathways are unmapped or unintegrated ####
nrow(pathabundance[grepl("UNINTEGRATED|UNMAPPED",pathabundance$Pathway),])
pathabundance[grepl("UNINTEGRATED|UNMAPPED",pathabundance$Pathway),1]

nrow(pathabundance[!grepl("UNINTEGRATED|UNMAPPED",pathabundance$Pathway),])

# pathabundance.col_sums.total<-colSums(pathabundance[,-1])
# pathabundance.col_sums.total
# pathabundance.col_sums.unmapped<-colSums(pathabundance[
#   grepl("UNINTEGRATED$|UNMAPPED$",pathabundance$Pathway),-1])
# pathabundance.col_sums.unmapped/pathabundance.col_sums.total*100
# pathabundance.col_sums.total-pathabundance.col_sums.unmapped

# Remove unmapped and unintegrated pathways ####
# pathabundance<-pathabundance[!grepl("UNINTEGRATED|UNMAPPED",pathabundance$Pathway),]

# Split the Pathway column ####
# Transform into a long format but don't touch the Pathway column
pathabundance<-pathabundance%>%
  pivot_longer(!Pathway,
               names_to = "Sample")
# Split the Pathway column by ":"
# For example, ARO-PWY: chorismate biosynthesis I|g__Phascolarctobacterium.s__Phascolarctobacterium_faecium
# will be split into two columns: pathabundance$Pathway$X1 and pathabundance$Pathway$X2
# pathabundance$Pathway$X1 will have ARO-PWY and pathabundance$Pathway$X2 will
# have chorismate biosynthesis I|g__Phascolarctobacterium.s__Phascolarctobacterium_faecium
# But those that don't have a ":" will be duplicated (both columns will have the 
# same contents)
pathabundance<-within(pathabundance, 
            Pathway<-data.frame(do.call('rbind', 
                                             strsplit(as.character(Pathway), 
                                                      ':', fixed=TRUE))))
# copy pathabundance to another dataframe
# pathabundance_full<-pathabundance
# pathway_strat has the pathway name with stratification
# pathabundance_full$pathway_strat<-pathabundance$Pathway$X2
# Pathway_ID has the Metacyc ID or EC number of the pathway
# pathabundance_full$Pathway_ID<-pathabundance$Pathway$X1
# pathabundance_full<-pathabundance_full%>%
#   select(-Pathway)
# Split the pathway_strat column by "|": all stratifications will be put
# into a separate column. 
# So, from chorismate biosynthesis I|g__Phascolarctobacterium.s__Phascolarctobacterium_faecium
# pathabundance_full$pathway_strat$X1 will have chorismate biosynthesis I
# and pathabundance_full$pathway_strat$X2 will have g__Phascolarctobacterium.s__Phascolarctobacterium_faecium
pathabundance_full<-within(pathabundance_full, 
                      pathway_strat<-data.frame(do.call('rbind', 
                                                  strsplit(as.character(pathway_strat), 
                                                           '|', fixed=TRUE))))
# pathway column has the name of the pathway
# pathabundance_full$pathway<-pathabundance_full$pathway_strat$X1
# # strat column has the stratification
# pathabundance_full$strat<-pathabundance_full$pathway_strat$X2
# pathabundance_full<-pathabundance_full%>%
#   select(-pathway_strat)
# 
# pathabundance$Pathway_ID<-pathabundance$Pathway$X1
# 
# pathabundance_full%>%
#   filter(pathway==strat|strat=="unclassified")%>%
#   group_by(Sample,Pathway_ID)%>%
#   summarise(Abundance=sum(value))%>%head
# 
# pathabundance_full%>%
#   filter(pathway!=strat&strat!="unclassified")%>%
#   group_by(Sample,Pathway_ID)%>%
#   summarise(Abundance=sum(value))%>%head
# 
# 
# pathabundance<-pathabundance%>%
#   group_by(Sample,Pathway_ID)%>%
#   summarise(Abundance=sum(value))

# pathabundance.wide<-pathabundance%>%
#   pivot_wider(names_from = Sample,
#               values_from = Abundance)
# colSums(pathabundance.wide[-1])
# pathabundance.clr<-clr(pathabundance.wide[,-1])


############### VVVVVVVVVV
pathabundance.wide<-pathabundance%>%
  pivot_wider(names_from = Pathway,
              values_from = value,
              values_fill = 0)%>%
  column_to_rownames("Sample")
########## ^^^^^

# pathabundance.wide<-pathabundance.wide%>%column_to_rownames("Pathway_ID")
# Check the number of unique pathways ####
# nrow(pathabundance.wide)


# Gene families ####
gene_family.filename<-
  file.path("output/humann3_pipeline",paste("humann_out",
                                            merged.date_time,"merged",
                                            sep = "_"),
            paste(merged.date_time,
                  "humann_out_genefamilies-cpm-rxn-named-filtered-total-renorm.tsv",sep = "_"))

gene_family<-read.table(gene_family.filename, 
                          header = T, 
                          sep = "\t",
                          fill=TRUE,
                          quote = "", 
                          comment.char = "@",
                          na.strings = c("","NA"))
# gene_family.long<-gene_family%>%
#   pivot_longer(!X..Gene.Family,
#                names_to = "Sample")
# gene_family.long%>%
#   filter(grepl("ASPARTATEKIN-RXN: (expasy) Aspartate kinase [2.7.2.4]",
#                X..Gene.Family,fixed = T),Sample=="X2D10_trim_merged_Abundance.CPM")%>%
#   group_by(Sample)%>%
#   summarise(sumab=sum(value))


dim(gene_family)
head(gene_family)
# gene_family[is.na(gene_family)]<-"NA"
colnames(gene_family)<-gsub("X..Gene.Family","Gene.Family",colnames(gene_family))
colnames(gene_family)<-gsub("X2D10","2D10",colnames(gene_family))
colnames(gene_family)<-gsub("X2D14","2D14",colnames(gene_family))
colnames(gene_family)<-
  gsub("_trim_merged_Abundance","",colnames(gene_family))
colnames(gene_family)<-
  gsub(".CPM","",colnames(gene_family))

# Check how many gene families are unmapped or ungrouped ####
nrow(gene_family[grepl("UNGROUPED|UNMAPPED",gene_family$Gene.Family),])
nrow(gene_family[!grepl("UNGROUPED|UNMAPPED",gene_family$Gene.Family),])

# gene_family.col_sums.total<-colSums(gene_family[,-1])
# gene_family.col_sums.unmapped<-colSums(gene_family[
#   grepl("UNGROUPED|UNMAPPED",gene_family$Gene.Family),-1])
# gene_family.col_sums.unmapped/gene_family.col_sums.total*100
# gene_family.col_sums.total-gene_family.col_sums.unmapped
# 
# # Remove unmapped and ungrouped gene families ####
# gene_family<-gene_family[!grepl("UNGROUPED|UNMAPPED",gene_family$Gene.Family),]

# Split Gene.Family column ####
gene_family<-gene_family%>%
  pivot_longer(!Gene.Family,
               names_to = "Sample")
# Split the Gene.Family column by ":"
# For example, LINOLENOYL-RXN: (expasy) Long-chain-fatty-acid--CoA ligase [6.2.1.3]|g__Bacteroides.s__Bacteroides_xylanisolvens
# will be split into two columns: gene_family$Gene.Family$X1 and gene_family$Gene.Family$X2
# gene_family$Gene.Family$X1 will have LINOLENOYL-RXN and gene_family$Gene.Family$X2 will
# have (expasy) Long-chain-fatty-acid--CoA ligase [6.2.1.3]|g__Bacteroides.s__Bacteroides_xylanisolvens
# But those that don't have a ":" will be duplicated (both columns will have the 
# same contents)
# gene_family<-within(gene_family, 
#                     Gene.Family<-data.frame(do.call('rbind', 
#                                                   strsplit(as.character(Gene.Family), 
#                                                            ':', fixed=TRUE))))
# copy gene_family to another dataframe
# gene_family_full<-gene_family
# enzyme_strat has the enzyme name with stratification
# gene_family$enzyme_strat<-gene_family$Gene.Family$X2
# Gene.Family_ID has the Metacyc ID or EC number of the Gene Family 
# gene_family$Gene.Family_ID<-gene_family$Gene.Family$X1

# gene_family<-gene_family%>%
#   select(-Gene.Family)

# Split the enzyme_strat column by "|": all stratifications will be put
# into a separate column. 
# So, from (expasy) Long-chain-fatty-acid--CoA ligase [6.2.1.3]|g__Bacteroides.s__Bacteroides_xylanisolvens
# gene_family$enzyme_strat$X1 will have (expasy) Long-chain-fatty-acid--CoA ligase [6.2.1.3]
# and gene_family$enzyme_strat$X2 will have g__Bacteroides.s__Bacteroides_xylanisolvens
# gene_family<-within(gene_family, 
#                          enzyme_strat<-data.frame(do.call('rbind', 
#                                                             strsplit(as.character(enzyme_strat), 
#                                                                      '|', fixed=TRUE))))
# gene_family column has the name of the enzyme_strat
# gene_family$gene_family<-gene_family$enzyme_strat$X1
# strat column has the stratification
# gene_family$strat<-gene_family$enzyme_strat$X2
# gene_family<-gene_family%>%
#   select(-enzyme_strat)

# Verify that sums of gene family stratifications are the same as 
# total gene family abundance
# Extract totals
# gene_family.total<-gene_family%>%
#   filter(gene_family ==strat)%>%
#   group_by(Sample,Gene.Family_ID)%>%
#   summarise(Abundance=sum(value))
# Extract stratifications and sum by gene family
# gene_family.sum_strat<-gene_family%>%
#   filter(gene_family!=strat)%>%
#   group_by(Sample,Gene.Family_ID)%>%
#   summarise(Abundance=sum(value))
# check the differences
# gene_family.total%>%
#   left_join(gene_family.sum_strat,by=c("Sample","Gene.Family_ID"))%>%
#   mutate(checksum=abs(Abundance.x-Abundance.y))%>%
#   arrange(-checksum)
# 
# gene_family%>%
#   filter(gene_family !=strat)%>%View
# gene_family%>%
#   filter(gene_family ==strat)%>%View

# gene_family%>%
#   filter(gene_family ==strat)%>%
#   filter(Sample=="2D10")

# Duplicates???
# gene_family%>%
#   filter(gene_family ==strat)%>%
#   select(-strat,-gene_family)%>%
#   summarise(n = dplyr::n(), .by = c(Sample, Gene.Family_ID))%>%
#   arrange(-n)%>%View


# gene_family.wide<-gene_family%>%
#   filter(gene_family ==strat,value!=0)%>%
#   group_by(Sample,Gene.Family_ID)%>%
#   summarise(Abundance=sum(value))%>%
#   pivot_wider(names_from = Gene.Family_ID,
#               values_from = Abundance)

# Rows are samples and columns are features
gene_family.wide<-gene_family%>%
  pivot_wider(names_from = Gene.Family,
              values_from = value,
              values_fill = 0)%>%
  column_to_rownames("Sample")


# gene_family.sums.wide<-gene_family%>%
#   filter(gene_family!=strat,value!=0)%>%
#   group_by(Sample,Gene.Family_ID)%>%
#   summarise(Abundance=sum(value))%>%
#   pivot_wider(names_from = Gene.Family_ID,
#               values_from = Abundance)

# gene_family.clr<-clr(gene_family.wide[,-1])

rowSums(gene_family.wide)



# Create metadata ####
custom.md<-readRDS("../amplicon_nmr/output/rdafiles/custom.md.rds")
# metadata.filename<-paste0(metadatadir,
#                           paste("filenames-single-pooled-raw-supercomp.tsv", sep = "-"))
# custom.md<-read.table(metadata.filename, header = T)
# colnames(custom.md)[1]<-"Sample" # set the first column name as Sample
custom.md$Sample<-rownames(custom.md)


custom.md<-custom.md%>%
  filter(Sample%in%pathabundance$Sample)
custom.md<-custom.md%>%
  filter(Sample%in%gene_family$Sample)#%>%
  # column_to_rownames(var = "Sample")%>%
  # select(-absolute.filepath)
# custom.md<-custom.md%>%
#   filter(Sample%in%colnames(gene_family)[-1])%>%
#   column_to_rownames(var = "Sample")%>%
#   select(-absolute.filepath)
# 
# custom.md$class<-as.factor(custom.md$class)
# custom.md$animal<-as.factor(custom.md$animal)
# custom.md$sex<-as.factor(custom.md$sex)
# custom.md$birthday<-as.Date(custom.md$birthday)
custom.md$age<-year(as.period(interval(custom.md$birthday,as.Date("2023-11-16"))))


# age.breaks<-c(min(custom.md$age),10,max(custom.md$age))
age.breaks<-c(0,10,16)
custom.md<-custom.md%>%
  mutate(agegroup=cut(age, breaks = age.breaks, 
                       right=FALSE))%>%
  mutate(agegroup=as.factor(agegroup))
# substr(gsub("\\[|\\]","",as.character(custom.md$agegroup)[1]),2,2)

unique_levels <- custom.md %>%
  ungroup()%>%
  distinct(agegroup)%>%
  arrange(agegroup) %>%
  mutate(new_agegroup = paste0("agegroup", agegroup))%>%
  mutate(new_agegroup = gsub("\\(|\\)|\\[|\\]","",new_agegroup))%>%
  mutate(new_agegroup = gsub("\\,","_",new_agegroup))

custom.md <- custom.md %>%
  left_join(unique_levels, by = "agegroup")
colnames(custom.md)[which(colnames(custom.md)=="agegroup")]<-"old_agegroup"
colnames(custom.md)[which(colnames(custom.md)=="new_agegroup")]<-"agegroup"
rownames(custom.md)<-custom.md$Sample
custom.levels<-names(table(custom.md$agegroup))

maaslin.reference<-paste("agegroup",ref.level,sep = ",")
maaslin.comparison<-"agegroup"


# Add relations
relations<-read.table("../amplicon_nmr/data/metadata/pooled-metadata/nmr-relations.tsv",
                      header = T,
                      sep = "\t")
custom.md<-custom.md%>%
  left_join(relations,by="Sample")
rownames(custom.md)<-custom.md$Sample


# Maaslin2 on pathways ####
set.seed(1)
maaslin.fit_data =
  Maaslin2(input_data = pathabundance.wide,
           input_metadata = custom.md,
           min_prevalence = 0,
           analysis_method = "LM",
           normalization = "TSS",
           transform = "LOG",
           random_effects = c("relation"), 
           standardize = FALSE,
           output = file.path("./output/maaslin2","humann3-output",
                              paste(
                                paste(format(Sys.time(),format="%Y%m%d"),
                                      format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                "pathway",host,comparison,
                                paste(custom.levels,collapse = '-'),
                                "ref",ref.level,sep = "-")),
           fixed_effects = maaslin.comparison,
           reference = maaslin.reference,
           max_significance = 0.5)
maaslin.fit_data$results%>%head

save.image(file.path("./output/rdafiles",paste(
  paste(format(Sys.time(),format="%Y%m%d"),
        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
  "maaslin2-humann3-pathabundance",host,
  comparison,paste(custom.levels,collapse = '-'),"ref",
  ref.level,"workspace.RData",sep="-")))

# Maaslin2 on gene families ####
set.seed(1)
maaslin.fit_data = 
  Maaslin2(input_data = gene_family.wide, 
           input_metadata = custom.md, 
           min_prevalence = 0,
           analysis_method = "LM",
           normalization = "CLR",
           transform = "LOG",
           random_effects = c("relation"), 
           standardize = FALSE,
           output = file.path("./output/maaslin2","humann3-output",
                              paste(
                                paste(format(Sys.time(),format="%Y%m%d"),
                                      format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                "gene_family",host,comparison,
                                paste(custom.levels,collapse = '-'),
                                "ref",ref.level,sep = "-")),
           fixed_effects = maaslin.comparison,
           reference = maaslin.reference,
           max_significance = 0.5)
maaslin.fit_data$results%>%head

save.image(file.path("./output/rdafiles",paste(
  paste(format(Sys.time(),format="%Y%m%d"),
        format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
  "maaslin2-humann3-gene_family",host,
  comparison,paste(custom.levels,collapse = '-'),"ref",
  ref.level,"workspace.RData",sep="-")))

# Downstream analysis on gene families ####
maaslin.workspace.gene_family.datetime<-"20240621_18_59_28"
comparison<-"age"
host<-"NMR"
ref.level<-"agegroup0_10" # choose the reference level
custom.levels<-c("agegroup0_10","agegroup10_16")
load(file.path("./output/rdafiles",
            paste(maaslin.workspace.gene_family.datetime,"maaslin2-humann3-gene_family",host,
                  comparison,paste(custom.levels,collapse = '-'),"ref",
                  ref.level,
                  "workspace.RData",sep="-")))
source("../amplicon_nmr/code/r-scripts/make_features_maaslin.R")

if(min(maaslin.fit_data$results$qval)<0.05){
  maaslin.signif.features<-maaslin.fit_data$results%>%
    filter(qval<0.05) # should be qval
}else{
  maaslin.signif.features<-maaslin.fit_data$results%>%
    arrange(qval)%>%
    head(n = 10) # if no significant results found
}


foo<-gene_family
foo$maaslin<-foo$Gene.Family
foo<-make_features_maaslin(foo,"maaslin")
foo<-unique(foo[,c("maaslin","Gene.Family")])
maaslin.signif.features<-maaslin.signif.features%>%
  left_join(foo[,c("maaslin","Gene.Family")],by=c("feature"="maaslin"))%>%
  distinct()
rm(foo)
maaslin.signif.features$feature<-maaslin.signif.features$Gene.Family
maaslin.signif.features<-subset(maaslin.signif.features, select=-Gene.Family)

table(maaslin.signif.features$feature%in%gene_family$Gene.Family)

maaslin.signif.features<-maaslin.signif.features%>%
  left_join(unique(gene_family[,c("Gene.Family","enzyme_strat")]),
            by=c("feature"="Gene.Family"))


maaslin.signif.decreased<-maaslin.signif.features%>%
  as_tibble()%>%
  filter(coef<0)%>%
  arrange(feature)%>%
  # mutate(n=n())%>%
  mutate(assoc.str=-log(qval)*sign(coef))#%>%
# select(feature,assoc.str,name)

maaslin.signif.increased<-maaslin.signif.features%>%
  as_tibble()%>%
  filter(coef>0)%>%
  arrange(feature)%>%
  # mutate(n=n())%>%
  mutate(assoc.str=-log(qval)*sign(coef))#%>%
# select(feature,assoc.str,name)

table(maaslin.signif.decreased$feature%in%gene_family$Gene.Family)
table(maaslin.signif.increased$feature%in%gene_family$Gene.Family)

maaslin.signif.decreased$feature[!maaslin.signif.decreased$feature%in%gene_family$Gene.Family]
maaslin.signif.increased$feature[!maaslin.signif.increased$feature%in%gene_family$Gene.Family]


maaslin.features.list<-maaslin.signif.increased%>%
  bind_rows(maaslin.signif.decreased)%>%
  filter(abs(assoc.str)>1)%>%
  select(feature,coef,assoc.str,metadata,value,enzyme_strat)
# write.table(maaslin.features.list,
#             file = paste0("./output/rtables/",
#                           paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                       format(Sys.time(),format = "%H_%M_%S"),
#                                       sep = "_"),
#                           "maaslin2-humann3-signif-gene_family",
#                           paste(custom.levels,collapse = '-'),
#                           "ref",ref.level,sep = "-"),
#                           ".tsv"),
#             row.names = F,
#             col.names = F,
#             quote = F,
#             sep = "\t")
# 
# write.table(maaslin.features.list$feature,
#             file = paste0("./output/rtables/",
#                           paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                 format(Sys.time(),format = "%H_%M_%S"),
#                                 sep = "_"),
#                           "maaslin2-humann3-signif-gene_family-list",
#                           paste(custom.levels,collapse = '-'),
#                           "ref",ref.level,sep = "-"),".tsv"),
#             row.names = F,
#             col.names = F,
#             quote = F,
#             sep = "\t")


tot.pathw<-maaslin.features.list%>%
  group_by(enzyme_strat)%>%
  mutate(tot=n())%>%
  arrange(-tot)
# write.table(tot.pathw,
#             file = paste0("./output/rtables/",
#                           paste(paste(format(Sys.time(),format="%Y%m%d"),
#                                 format(Sys.time(),format = "%H_%M_%S"),
#                                 sep = "_"),
#                           "maaslin2-humann3-signif-gene_family-tot.pathw",
#                           paste(custom.levels,collapse = '-'),"ref",
#                           ref.level,sep = "-"),".tsv"),
#             row.names = F,
#             quote = F,
#             sep = "\t")
# Barplot for gene families
level.names<-names(sort(table(maaslin.features.list$enzyme_strat),decreasing = F))
level.names<-gsub("^ \\(expasy\\)|^ \\(metacyc\\)","",level.names)
level.names<-gsub("^ ", "",level.names)
level.names<-str_to_sentence(level.names)

# Frequency plot for gene families that belong to the same pathway
maaslin.gene_family.freq.plot<-maaslin.features.list%>%
  mutate(new_enzyme_strat=gsub("^ \\(expasy\\)|^ \\(metacyc\\)","",enzyme_strat))%>%
  mutate(new_enzyme_strat=gsub("^ ","",new_enzyme_strat))%>%
  mutate(new_enzyme_strat=str_to_sentence(new_enzyme_strat))%>%
  
  mutate(new_enzyme_strat=factor(new_enzyme_strat,levels=level.names))%>%
  ggplot(aes(y=new_enzyme_strat))+
  geom_bar()+
  theme_bw()+
  labs(x="",
       y="",
       title = "Gene family abundance from MaAsLin2")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = ggtext::element_markdown(hjust=1,size=18),
        axis.text.y = element_text(size=20),
        axis.title = element_text(size = 20),
        strip.text.x = element_text(size=20),
        plot.title = element_text(size = 27),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        legend.position = "none")
ggsave(paste0("./images/barplots/",
              paste( paste(format(Sys.time(),format="%Y%m%d"),
                                 format(Sys.time(),format = "%H_%M_%S"),
                           sep = "_"),"barplot",
                    "humann3-maaslin2-gene_family.freq",sep = "-"),
              ".png"),
       plot = maaslin.gene_family.freq.plot,
       width = 6000,height = 3000,
       units = "px",dpi=300,device = "png")



# Downstream analysis on pathabundance ####
maaslin.workspace.pathabundance.datetime<-"20240524_17_30_26"
comparison<-"age"
host<-"NMR"
ref.level<-"agegroup0_10" # choose the reference level
# this is for file names
if(host=="NMR"){
  host.labels<-c("NMR" = "*Heterocephalus glaber*")
}else{
  host.labels<-
    c("B6mouse" = "B6 mouse",
      "MSMmouse" = "MSM/Ms mouse",
      "FVBNmouse" = "FVB/N mouse")
}
custom.levels<-c("agegroup0_10","agegroup10_16")
load(file.path("./output/rdafiles",
               paste(maaslin.workspace.pathabundance.datetime,"maaslin2-humann3-pathabundance",host,
                     comparison,paste(custom.levels,collapse = '-'),"ref",
                     ref.level,
                     "workspace.RData",sep="-")))

source("../amplicon_nmr/code/r-scripts/make_features_maaslin.R")
if(min(maaslin.fit_data$results$qval)<0.05){
  maaslin.signif.features<-maaslin.fit_data$results%>%
    filter(qval<0.05) # should be qval
}else{
  maaslin.signif.features<-maaslin.fit_data$results%>%
    arrange(qval)%>%
    head(n = 10) # if no significant results found
}

foo<-pathabundance
foo$maaslin<-foo$Pathway_ID
foo<-make_features_maaslin(foo,"maaslin")
foo<-unique(foo[,c("maaslin","Pathway_ID")])
maaslin.signif.features<-maaslin.signif.features%>%
  left_join(foo[,c("maaslin","Pathway_ID")],by=c("feature"="maaslin"))%>%
  distinct()
rm(foo)
maaslin.signif.features$feature<-maaslin.signif.features$Pathway_ID
maaslin.signif.features<-subset(maaslin.signif.features, select=-Pathway_ID)


maaslin.signif.features<-maaslin.signif.features%>%
  left_join(unique(pathabundance_full[,c("Pathway_ID","pathway_strat")]),by=c("feature"="Pathway_ID"))


table(maaslin.signif.features$feature%in%pathabundance$Pathway_ID)

maaslin.signif.decreased<-maaslin.signif.features%>%
  as_tibble()%>%
  filter(coef<0)%>%
  arrange(feature)%>%
  # mutate(n=n())%>%
  mutate(assoc.str=-log(qval)*sign(coef))#%>%
# select(feature,assoc.str,name)

maaslin.signif.increased<-maaslin.signif.features%>%
  as_tibble()%>%
  filter(coef>0)%>%
  arrange(feature)%>%
  # mutate(n=n())%>%
  mutate(assoc.str=-log(qval)*sign(coef))#%>%
# select(feature,assoc.str,name)

table(maaslin.signif.decreased$feature%in%pathabundance$Pathway_ID)
table(maaslin.signif.increased$feature%in%pathabundance$Pathway_ID)

maaslin.signif.decreased$feature[!maaslin.signif.decreased$feature%in%pathabundance$Pathway_ID]
maaslin.signif.increased$feature[!maaslin.signif.increased$feature%in%pathabundance$Pathway_ID]
maaslin.features.list<-maaslin.signif.increased%>%
  bind_rows(maaslin.signif.decreased)%>%
  filter(abs(assoc.str)>1)


write.table(maaslin.features.list,
            file = paste0("./output/rtables/",
                          paste(paste(format(Sys.time(),format="%Y%m%d"),
                                      format(Sys.time(),format = "%H_%M_%S"),
                                      sep = "_"),
                          "maaslin2-humann3-signif-pathabundance",
                          paste(custom.levels,collapse = '-'),
                          "ref",ref.level,sep = "-"),
                          ".tsv"),
            row.names = F,
            col.names = F,
            quote = F,
            sep = "\t")

write.table(maaslin.features.list$feature,
            file = paste0("./output/rtables/",
                          paste(paste(format(Sys.time(),format="%Y%m%d"),
                                format(Sys.time(),format = "%H_%M_%S"),
                                sep = "_"),
                          "maaslin2-humann3-signif-pathabundance-list",
                          paste(custom.levels,collapse = '-'),
                          "ref",ref.level,sep = "-"),".tsv"),
            row.names = F,
            col.names = F,
            quote = F,
            sep = "\t")
            

# Gene families plot ####
merged.date_time<-"20240217"
image.formats<-c("png","tiff")
host<-"NMR"
comparison<-"age"
custom.levels<-c("between10and15","between2and10")
metadatadir<-paste0("../amplicon_nmr/data/metadata/pooled-metadata/") # directory with metadata
metadata.filename<-paste0(metadatadir,
                          paste("filenames-single-pooled-raw-supercomp.tsv", sep = "-"))
custom.md<-read.table(metadata.filename, header = T)
colnames(custom.md)[1]<-"Sample" # set the first column name as Sample
custom.md<-custom.md%>%
  filter(class=="NMR")%>%
  column_to_rownames(var = "Sample")%>%
  select(-absolute.filepath)
custom.md$Sample<-rownames(custom.md)
custom.md$class<-as.factor(custom.md$class)
custom.md$animal<-as.factor(custom.md$animal)
custom.md$sex<-as.factor(custom.md$sex)
custom.md$birthday<-as.Date(custom.md$birthday)
custom.md$age<-year(as.period(interval(custom.md$birthday,as.Date("2023-11-16"))))
age.breaks<-c(min(custom.md$age),10,max(custom.md$age))
custom.md<-custom.md%>%
  mutate(agegroup=cut(age, breaks = age.breaks, 
                       include.lowest = TRUE))%>%
  mutate(agegroup=as.factor(agegroup))

gene_family.filename<-
  file.path("output/humann3_pipeline",paste("humann_out",
                                            merged.date_time,"merged",
                                            sep = "_"),
            paste(merged.date_time,
                  "humann_out_genefamilies-cpm-rxn-named.tsv",sep = "_"))
gene_family<-read.table(gene_family.filename, 
                        header = T, 
                        sep = "\t",
                        fill=TRUE,
                        quote = "", 
                        comment.char = "@",
                        na.strings = c("","NA"))
dim(gene_family)
head(gene_family)
# gene_family[is.na(gene_family)]<-"NA"
colnames(gene_family)<-gsub("X..Gene.Family","Gene.Family",colnames(gene_family))
colnames(gene_family)<-gsub("X2D10","2D10",colnames(gene_family))
colnames(gene_family)<-gsub("X2D14","2D14",colnames(gene_family))
colnames(gene_family)<-
  gsub("_trim_merged_Abundance","",colnames(gene_family))
colnames(gene_family)<-
  gsub(".CPM","",colnames(gene_family))
gene_family<-gene_family%>%
  pivot_longer(!Gene.Family,
               names_to = "Sample")
gene_family<-within(gene_family, 
                    Gene.Family<-data.frame(do.call('rbind', 
                                                    strsplit(as.character(Gene.Family), 
                                                             ':', fixed=TRUE))))
gene_family_full<-gene_family
gene_family_full$enzyme_strat<-gene_family$Gene.Family$X2
gene_family_full$Gene.Family<-gene_family_full$Gene.Family$X1
gene_family_full<-within(gene_family_full, 
                         enzyme_strat<-data.frame(do.call('rbind', 
                                                              strsplit(as.character(enzyme_strat), 
                                                                       '|', fixed=TRUE))))



foo<-gene_family_full
gene_family_full$pathway_strat<-gene_family_full$enzyme_strat$X1
gene_family_full$strat<-gene_family_full$enzyme_strat$X2
gene_family_full<-gene_family_full%>%
  select(-enzyme_strat)


sample.factors<-gene_family_full%>%
  left_join(custom.md,by="Sample")%>%
  filter(Gene.Family=="LINOLENOYL-RXN")%>%
  filter(!strat%in%" (expasy) Long-chain-fatty-acid--CoA ligase [6.2.1.3]")%>%
  group_by(Sample)%>%
  mutate(total.ab=sum(value))%>%
  distinct(total.ab)%>%
  arrange(-total.ab)%>%
  pull(Sample)
sample.factors


species.groups<-read.table("./output/rtables/species_groups.tsv", 
                           header = T, 
                           sep = "\t",
                           fill=TRUE,
                           quote = "", 
                           comment.char = "@",
                           na.strings = c("","NA"))
library(Polychrome)
library(RColorBrewer)
library(ggtext)
group.sum<-
  species.groups%>%
  group_by(Group)%>%
  arrange(Group)%>%
  summarise(n=n())
cols<-c("Bacteroides xylanisolvens"="#0000FF",
        "Hungatella hathewayi"="#8DEEEE",
        "Bacteroides thetaiotaomicron"="#556B2F",
        "Bacteroides uniformis"="#A2CD5A",
        "Intestinimonas butyriciproducens"="#CD6600",
        "Clostridium clostridioforme"="#FFB90F",
        "Clostridium symbiosum"="#EEDC82",
        "Weissella cibaria"="#8B0000",
        "Alistipes shahii" = "#4A1486",
        "Alistipes timonensis"= "#8B2252",
        "Bacteroides faecichinchillae" = "#9F79EE",
        "Eisenbergiella tayi" = "#FF00FF",
        "Faecalicatena contorta" = "#CD6090",
        "Clostridium citroniae" ="#BF3EFF",
        "Parabacteroides distasonis" ="#68228B",
        "Parasutterella excrementihominis" ="#FECAB1",
        "unclassified"="#696969")
names(cols)<-paste0("<i>",names(cols),"</i>")

humann3.plot<-gene_family_full%>%
  left_join(custom.md,by="Sample")%>%
  filter(Gene.Family=="LINOLENOYL-RXN")%>%
  filter(!strat%in%" (expasy) Long-chain-fatty-acid--CoA ligase [6.2.1.3]")%>%
  mutate(strat=gsub("g__.*\\.s__","",strat))%>%
  mutate(strat=gsub("_"," ",strat))%>%
  mutate(istrati=paste0("<i>",strat,"</i>"))%>%
  left_join(species.groups,by=c("strat"="Species"))%>%
  mutate(Sample=factor(Sample,levels=sample.factors))%>%
  ggplot(aes(x=Sample,y=value,fill=factor(istrati,levels=names(cols))))+
  geom_bar(stat="identity")+
  theme_bw()+
  facet_grid(~agegroup,
             scales = "free_x")+
  scale_fill_manual(values = cols,
                    breaks = factor(names(cols),levels=names(cols)))+
  guides(fill=guide_legend(ncol=1,title = "<b>Contributing species:</b>"))+
  labs(x="",y="Abundance (unspecified units)")+
  theme(strip.text.x = ggtext::element_markdown(size = 20),# the name of
        # each facet will be recognised as a markdown object, so we can 
        # add line breaks (cause host names are too long)
        panel.spacing = unit(0.8, "cm"), # increase distance between facets
        axis.text.x = element_text(size=20,color = "black"),# rotate
        # the x-axis labels by 45 degrees and shift to the right
        axis.text.y = element_text(size=20,color="black"), # size of y axis ticks
        axis.title = element_text(size = 20), # size of axis names
        plot.title = element_text(size = 25), # size of plot title
        plot.caption = element_text(size=23), # size of plot caption
        legend.text = element_markdown(size = 20), # size of legend text
        legend.title = element_markdown(size = 25), # size of legend title
        legend.position = "bottom")
for (image.format in image.formats){
  ggsave(paste0("./images/barplots/",
                paste(paste(format(Sys.time(),format="%Y%m%d"),
                            format(Sys.time(),format = "%H_%M_%S"),
                            sep = "_"),"humann3.plot-linolenoyl-rxn",
                      host,comparison, paste(custom.levels,collapse = '-'),
                      ref.level,sep = "-"),".",image.format),
         plot = humann3.plot,
         width = 6000,height = 3000,
         units = "px",dpi=300,device = image.format)
  
}
