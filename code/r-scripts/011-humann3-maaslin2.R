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
# pathabundance<-within(pathabundance, 
#             Pathway<-data.frame(do.call('rbind', 
#                                              strsplit(as.character(Pathway), 
#                                                       ':', fixed=TRUE))))
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
# pathabundance_full<-within(pathabundance_full, 
#                       pathway_strat<-data.frame(do.call('rbind', 
#                                                   strsplit(as.character(pathway_strat), 
#                                                            '|', fixed=TRUE))))
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
custom.md<-readRDS("./output/rdafiles/custom.md.ages.rds")
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
# custom.md$age<-year(as.period(interval(custom.md$birthday,as.Date("2023-11-16"))))


# age.breaks<-c(min(custom.md$age),10,max(custom.md$age))
# age.breaks<-c(0,10,16)
# custom.md<-custom.md%>%
#   mutate(agegroup=cut(age, breaks = age.breaks, 
#                        right=FALSE))%>%
#   mutate(agegroup=as.factor(agegroup))
# substr(gsub("\\[|\\]","",as.character(custom.md$agegroup)[1]),2,2)

# unique_levels <- custom.md %>%
#   ungroup()%>%
#   distinct(agegroup)%>%
#   arrange(agegroup) %>%
#   mutate(new_agegroup = paste0("agegroup", agegroup))%>%
#   mutate(new_agegroup = gsub("\\(|\\)|\\[|\\]","",new_agegroup))%>%
#   mutate(new_agegroup = gsub("\\,","_",new_agegroup))

# custom.md <- custom.md %>%
#   left_join(unique_levels, by = "agegroup")
# colnames(custom.md)[which(colnames(custom.md)=="agegroup")]<-"old_agegroup"
# colnames(custom.md)[which(colnames(custom.md)=="new_agegroup")]<-"agegroup"
# rownames(custom.md)<-custom.md$Sample
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

# save.image(file.path("./output/rdafiles",paste(
#   paste(format(Sys.time(),format="%Y%m%d"),
#         format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#   "maaslin2-humann3-gene_family",host,
#   comparison,paste(custom.levels,collapse = '-'),"ref",
#   ref.level,"workspace.RData",sep="-")))

# Downstream analysis on gene families ####
library(tidyverse)
maaslin.workspace.gene_family.datetime<-"20240621_18_53_46"
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


# Split the feature column by ":"
# For example, LINOLENOYL-RXN: (expasy) Long-chain-fatty-acid--CoA ligase [6.2.1.3]|g__Bacteroides.s__Bacteroides_xylanisolvens
# will be split into two columns: gene_family$Gene.Family$X1 and gene_family$Gene.Family$X2
# gene_family$Gene.Family$X1 will have LINOLENOYL-RXN and gene_family$Gene.Family$X2 will
# have (expasy) Long-chain-fatty-acid--CoA ligase [6.2.1.3]|g__Bacteroides.s__Bacteroides_xylanisolvens
# But those that don't have a ":" will be duplicated (both columns will have the 
# same contents)
maaslin.signif.features<-within(maaslin.signif.features,
                    feature<-data.frame(do.call('rbind',
                                                  strsplit(as.character(feature),
                                                           ':', fixed=TRUE))))



write.table(maaslin.signif.features$feature$X1,
            file=file.path("./output/rtables",
                           paste(
                             paste(format(Sys.time(),format="%Y%m%d"),
                                   format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                             "maaslin2-singif-gene_family",host,comparison,
                             paste(custom.levels,collapse = '-'),
                             "ref",ref.level,".tsv",sep = "-")),
            row.names = F,col.names = F,quote = F,sep = "\t")

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
maaslin.workspace.pathabundance.datetime<-"20240621_19_16_03"
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
custom.levels<-c("between10and16","between0and10")
metadatadir<-paste0("../amplicon_nmr/data/metadata/pooled-metadata/") # directory with metadata
metadata.filename<-paste0(metadatadir,
                          paste("filenames-single-pooled-raw-supercomp.tsv", sep = "-"))
custom.md<-readRDS("./output/rdafiles/custom.md.ages.rds")

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


# Dominant gene families and pathways? ####
gene_family<-gene_family%>%
  left_join(custom.md,by="Sample")
# Gene families are already in CPM units
gene_family<-gene_family%>%
  rename("AbundanceCPM"=value)%>%
  group_by(class,Sample)%>%
  mutate(TotalSample=sum(AbundanceCPM))%>%
  group_by_at(c("class","Sample","Gene.Family"))%>%
  mutate(RelativeAbundanceCPM=AbundanceCPM/TotalSample*100)%>%
  group_by(class)%>% # group by class (animal host),
  mutate(TotalClass=sum(AbundanceCPM))%>%
  group_by_at(c("class","Gene.Family"))%>%
  mutate(TotalGene.Family=sum(AbundanceCPM))%>%
  mutate(MeanRelativeAbundanceCPM=TotalGene.Family/TotalClass*100)#%>%
  # select(-TotalClass,-TotalSample)

  
# gene_family<-gene_family%>%
#   group_by(agegroup)%>% # group by class (animal host),
#   mutate(TotalAgegroup=sum(AbundanceCPM))%>%
#   group_by(agegroup,Species)%>%
#   mutate(TotalGene.FamilyAge=sum(AbundanceCPM))%>%
#   mutate(MeanRelativeAbundanceCPMAgegroup=TotalGene.FamilyAge/TotalAgegroup*100)


dominant.gene_families<-gene_family%>%
  group_by(Gene.Family)%>%
  arrange(-MeanRelativeAbundanceCPM)%>%
  distinct(Gene.Family,.keep_all = T)%>%
  select(-TotalClass,-TotalSample,-TotalGene.Family,-class,-animal)%>%
  head(10)%>%
  pull(Gene.Family)

pretty.level.names<-names(table(custom.md$old_agegroup))
names(pretty.level.names)<-names(table(custom.md$agegroup))
custom.fill<-c("#EE2C2C","#5CACEE")
names(custom.fill)<-custom.levels

sample.levels<-custom.md%>%
  select(Sample,agegroup)%>%
  arrange(agegroup)%>%
  distinct()
sample.levels$Sample<-factor(sample.levels$Sample,
                             levels=sample.levels$Sample)

gene_family%>%
  filter(Gene.Family%in%dominant.gene_families)%>%
  mutate(Gene.Family=factor(Gene.Family,levels=dominant.gene_families))%>%
  ggplot(aes(x=Sample,y=RelativeAbundanceCPM,fill=factor(agegroup)))+
  geom_bar(stat="identity")+
  facet_wrap(~Gene.Family,
             scales = "free",
             ncol=2)+
  theme_bw()+
  scale_color_manual(breaks = pretty.level.names,
                     labels=unname(pretty.level.names))+
  scale_x_discrete(labels=pretty.level.names,
                   limits=sample.levels$Sample)+
  scale_fill_manual(values = custom.fill,
                    labels=pretty.level.names)



gene_family%>%
  group_by(Sample)%>%
  summary()



# Sulfur metabolism ####
gene_family.sulfur<-gene_family[grepl("sulf|thio",gene_family$Gene.Family),]
gene_family.sulfur.wide<-gene_family.sulfur%>%
  ungroup%>%
  select(Sample,Gene.Family,AbundanceCPM)%>%
  pivot_wider(names_from = "Gene.Family",
              values_from = "AbundanceCPM")%>%
  column_to_rownames("Sample")

gene_family.sulfur%>%
  group_by(Gene.Family)%>%
  arrange(-MeanRelativeAbundanceCPM)%>%
  distinct(Gene.Family,.keep_all = T)%>%
  select(-TotalClass,-TotalSample,-TotalGene.Family,-class,-animal)%>%
  head(10)%>%
  select(Gene.Family,MeanRelativeAbundanceCPM)

gene_family.sulfur.split<-
  within(gene_family.sulfur,Gene.Family<-
           data.frame(do.call('rbind',strsplit(as.character(Gene.Family),
                                               '\\(expasy\\)|\\(metacyc\\)'))))
gene_family.sulfur.split$fullGeneFamily<-gene_family.sulfur.split$Gene.Family$X2
gene_family.sulfur.split$IDGeneFamily<-gene_family.sulfur.split$Gene.Family$X1
gene_family.sulfur.split<-gene_family.sulfur.split%>%
  relocate(IDGeneFamily,.before = everything())%>%
  relocate(fullGeneFamily,.after=IDGeneFamily)
gene_family.sulfur.split<-gene_family.sulfur.split%>%
  ungroup()%>%
  select(-Gene.Family)
gene_family.sulfur.split%>%
  group_by(fullGeneFamily,agegroup)%>%
  distinct(IDGeneFamily,.keep_all = T)%>%
  select(IDGeneFamily,fullGeneFamily)%>%
  arrange(fullGeneFamily)%>%
  tally%>%
  arrange(-n)

# Total CPM % taken by sulfur metabolism
gene_family.sulfur%>%
  group_by(Sample)%>%
  summarise(sumab=sum(RelativeAbundanceCPM))

gene_family.sulfur.split%>%
  group_by(fullGeneFamily,agegroup)%>%
  summarise(sumab=sum(RelativeAbundanceCPM))%>%
  arrange(-sumab)


gene_family.sulfur.split_wide<-gene_family.sulfur.split%>%
  ungroup%>%
  select(Sample,fullGeneFamily,AbundanceCPM)%>%
  group_by(Sample,fullGeneFamily)%>%
  summarise(AbundanceCPM=sum(AbundanceCPM))%>%
  pivot_wider(names_from = "fullGeneFamily",
              values_from = "AbundanceCPM")%>%
  column_to_rownames("Sample")

young.ECs.sulfur<-gene_family.sulfur.split%>%
  filter(agegroup=="agegroup0_10",AbundanceCPM!=0)%>%
  distinct(fullGeneFamily)%>%
  pull

old.ECs.sulfur<-gene_family.sulfur.split%>%
  filter(agegroup=="agegroup10_16",AbundanceCPM!=0)%>%
  distinct(fullGeneFamily)%>%
  pull

young.IDs.sulfur<-gene_family.sulfur.split%>%
  filter(agegroup=="agegroup0_10",AbundanceCPM!=0)%>%
  distinct(IDGeneFamily)%>%
  pull

old.IDs.sulfur<-gene_family.sulfur.split%>%
  filter(agegroup=="agegroup10_16",AbundanceCPM!=0)%>%
  distinct(IDGeneFamily)%>%
  pull

setdiff(young.ECs.sulfur,old.ECs.sulfur)
setdiff(old.ECs.sulfur,young.ECs.sulfur)

gene_family.sulfur.uniq.vector<-c(setdiff(young.ECs.sulfur,old.ECs.sulfur),
                                  setdiff(old.ECs.sulfur,young.ECs.sulfur))
gene_family.sulfur.uniq.vector<-gene_family.sulfur.split%>%
  arrange(agegroup,-MeanRelativeAbundanceCPM)%>%
  filter(fullGeneFamily%in%gene_family.sulfur.uniq.vector,AbundanceCPM!=0)%>%
  select(fullGeneFamily)%>%
  distinct()%>%
  pull()

gene_family.sulfur.split%>%
  filter(fullGeneFamily%in%gene_family.sulfur.uniq.vector)%>%
  mutate(fullGeneFamily=factor(fullGeneFamily,levels=gene_family.sulfur.uniq.vector))%>%
  ggplot(aes(x=Sample,y=AbundanceCPM,fill=factor(agegroup)))+
  geom_bar(stat="identity")+
  facet_wrap(~fullGeneFamily,
             scales = "free",
             ncol=2)+
  theme_bw()+
  scale_color_manual(breaks = pretty.level.names,
                     labels=unname(pretty.level.names))+
  scale_x_discrete(labels=pretty.level.names,
                   limits=sample.levels$Sample)+
  scale_fill_manual(values = custom.fill,
                    labels=pretty.level.names)+
  labs(y="Relative abundance (Copies per million)",
       fill="Age group")+
  theme(axis.title.y = element_text(size = 25),
        axis.title = element_text(size = 20),
        axis.text.y = ggtext::element_markdown(size=18),
        axis.text.x = element_text(size=20),
        strip.text.x = ggtext::element_markdown(size=20),
        plot.title = element_text(size = 27),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        legend.position = "right")



ggsave(paste0("./images/barplots/",
              paste(paste(format(Sys.time(),format="%Y%m%d"),
                          format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                    "humann3-sulfur-gene-families",
                    sep = "-"),".png"),
       plot=last_plot(),
       width = 6000,height = 4000,
       units = "px",dpi=300,device = "png")


# CSE and CBE
gene_family.sulfur.split%>%
  filter(grepl("Cystathionine gamma-lyase \\[4\\.4\\.1\\.1]|Cystathionine beta-synthase \\[4\\.2\\.1\\.22\\]",
               fullGeneFamily,))%>%
  group_by(fullGeneFamily,agegroup)%>%
  summarise(sumab=sum(AbundanceCPM))%>%
  arrange(-sumab)

gene_family.sulfur.split%>%
  filter(grepl("Cystathionine gamma-lyase \\[4\\.4\\.1\\.1]|Cystathionine beta-synthase \\[4\\.2\\.1\\.22\\]|3-mercaptopyruvate sulfurtransferase \\[2\\.8\\.1\\.2\\]",
               fullGeneFamily))%>%
  ggplot(aes(x=Sample,y=AbundanceCPM,fill=factor(agegroup)))+
  geom_bar(stat="identity")+
  facet_wrap(~fullGeneFamily,
             scales = "free",
             ncol=1)+
  theme_bw()+
  labs(y="Relative abundance (Copies per million)",
       fill="Age group")+
  scale_color_manual(breaks = pretty.level.names,
                     labels=unname(pretty.level.names))+
  scale_x_discrete(labels=pretty.level.names,
                   limits=sample.levels$Sample)+
  scale_fill_manual(values = custom.fill,
                    labels=pretty.level.names)+
  theme(axis.title.y = element_text(size = 25),
        axis.title = element_text(size = 20),
        axis.text.y = ggtext::element_markdown(size=18),
        axis.text.x = element_text(size=20),
        strip.text.x = ggtext::element_markdown(size=20),
        plot.title = element_text(size = 27),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        legend.position = "right")
ggsave(paste0("./images/barplots/",
              paste(paste(format(Sys.time(),format="%Y%m%d"),
                          format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                    "humann3-sulfur-gene-families-cse-cbe",
                    sep = "-"),".png"),
       plot=last_plot(),
       width =4000,height = 3000,
       units = "px",dpi=300,device = "png")


# Can we do alpha diversity on gene families? ####
library(vegan)
library(Polychrome)
gg.labs.name<-"Age group"
gg.title.groups<-"age groups"

pretty.level.names<-c("Young group","Old group")
names(pretty.level.names)<-names(table(custom.md.ages$agegroup))
custom.levels<-names(pretty.level.names)
host<-"NMR"
host.class<-c("NMR"="naked mole-rat")
custom.levels<-c("agegroup0_10","agegroup10_16")

metric.labs=c('sobs'="Richness \n(Observed species)",
              'shannon' = "Shannon",
              # 'simpson' = "Simpson",
              'invsimpson' = "Inverse \nSimpson")
plot.metrics<-c("sobs","shannon", # "simpson",
                "invsimpson") # metrics to plot
div.indices<-c("sobs","shannon",# "simpson",
               "invsimpson")

set.seed(1)
custom.fill<-createPalette(length(custom.levels),
                           seedcolors = c("#EE2C2C","#5CACEE","#00CD66",
                                          "#FF8C00","#BF3EFF", "#00FFFF",
                                          "#FF6EB4","#00EE00","#EEC900"))
names(custom.fill)<-custom.levels
swatch(custom.fill)

gene_family.df<-gene_family%>%
  select(Sample,Gene.Family,AbundanceCPM,class,sex,agegroup,old_agegroup)%>%
  filter(AbundanceCPM!=0)

all.div<-gene_family.df%>%
  group_by(Sample)%>%
  reframe(sobs=specnumber(AbundanceCPM), # richness (num of species)
          shannon=diversity(AbundanceCPM,index = "shannon"),
          # simpson=diversity(Abundance, index="simpson"),
          invsimpson=diversity(AbundanceCPM, index="invsimpson"), # inverse simpson
          tot=sum(AbundanceCPM),
          agegroup=agegroup)%>%
  group_by(Sample)%>%
  pivot_longer(cols=c(sobs,shannon,invsimpson),
               names_to="metric")%>%
  distinct()


kt.results<-data.frame(matrix(nrow = 2,ncol=length(div.indices)))
colnames(kt.results)<-div.indices
rownames(kt.results)<-c("statistic", "pvalue")

combinations<-combn(custom.levels,2) # all unique pairwise combinations
w.results<-data.frame(matrix(nrow = ncol(combinations),ncol=length(div.indices))) # ncol(combinations) pairwise comparisons
colnames(w.results)<-div.indices
w.results<-array(dim = c(length(table(combinations[1,])), # nrows
                         length(table(combinations[2,])), # ncols
                         length(div.indices)), # num of 2D arrays (stacking) 
                 dimnames = list(NULL, NULL, div.indices))



for (div.metric in div.indices) {
  metric.ds<-all.div%>%
    filter(metric==div.metric)%>%
    distinct()
  # perform kruskal-wallis test
  if(comparison=="age"){
    kt<-kruskal.test(value~agegroup,data=metric.ds)
    
  }else if (comparison=="sex"){
    kt<-kruskal.test(value~sex,data=metric.ds)
    
  }else if(comparison=="strain"){
    kt<-kruskal.test(value~class,data=metric.ds)
    
  }
  
  kt.results["statistic",div.metric]<- kt$statistic
  kt.results["pvalue",div.metric]<- kt$p.value
  # low pvalue indicates statistical difference between some groups
  # But which groups are different from others?
  # Perform pairwise wilcoxon test
  
  # NB: You must not run pairwise wilcoxon test unless kruskal wallis results 
  # are significant
  if(kt$p.value<0.05){
    if(comparison=="age"){
      w.test<-pairwise.wilcox.test(metric.ds$value,
                                   metric.ds$agegroup,
                                   p.adjust.method = "BH",
                                   exact=FALSE)
    }else if (comparison=="sex"){
      w.test<-pairwise.wilcox.test(metric.ds$value,
                                   metric.ds$sex,
                                   p.adjust.method = "BH",
                                   exact=FALSE)
      
    }else if(comparison=="strain"){
      w.test<-pairwise.wilcox.test(metric.ds$value,
                                   metric.ds$class,
                                   p.adjust.method = "BH",
                                   exact=FALSE)
    }
    
    w.results[,,div.metric]<-w.test$p.value
    
  }else(
    w.results[,,div.metric]<-matrix(data = "n.s.",
                                    nrow = nrow(w.test$p.value),
                                    ncol = ncol(w.test$p.value))
    
  )
  dimnames(w.results)[[1]]<-dimnames(w.test$p.value)[[1]] # change rownames of w.results
  dimnames(w.results)[[2]]<-dimnames(w.test$p.value)[[2]] # change colnames of w.results
}

kt.results
w.results
stopifnot(all(kt.results[2,]<0.05))


max.values<-all.div%>%
  group_by(metric)%>%
  summarise(max_val=max(value))

# Prepare data for plotting
if(comparison=="age"){
  all.div$agegroup<-factor(all.div$agegroup,levels=custom.levels)
}else if (comparison=="sex"){
  all.div$sex<-factor(all.div$sex,levels=custom.levels)
  
}else if(comparison=="strain"){
  all.div$agegroup<-factor(all.div$class,levels=custom.levels)
}

## Plot alpha diversity metrics ####
if(comparison=="age"){
  div.plot<-ggplot(all.div[all.div$metric %in%
                             plot.metrics,],
                   # aes(x=reorder(class,-value),y=value,fill=class))+
                   aes(x=factor(all.div$agegroup,
                                level=custom.levels),y=value,fill=factor(agegroup)))
}else if (comparison=="sex"){
  div.plot<-ggplot(all.div[all.div$metric %in%
                             plot.metrics,],
                   # aes(x=reorder(class,-value),y=value,fill=class))+
                   aes(x=factor(all.div$sex,
                                level=custom.levels),y=value,fill=factor(sex)))
}else if(comparison=="strain"){
  div.plot<-ggplot(all.div[all.div$metric %in%
                             plot.metrics,],
                   # aes(x=reorder(class,-value),y=value,fill=class))+
                   aes(x=factor(all.div$class,
                                level=custom.levels),y=value,fill=factor(class)))
}

div.plot<-div.plot+
  geom_boxplot(show.legend = FALSE,alpha=0.8)+
  facet_wrap(~factor(metric, # reorder facets
                     levels=plot.metrics),
             ncol=length(plot.metrics),
             scales="free_y", # free y axis range
             labeller = as_labeller(metric.labs))+ # rename facets
  theme_bw()+ 
  labs(color=gg.labs.name)+
  scale_color_manual(breaks = unname(pretty.level.names),
                     labels=unname(pretty.level.names))+
  scale_x_discrete(labels=pretty.level.names,
                   limits=custom.levels)+ # rename boxplot labels (x axis)
  scale_fill_manual(values = custom.fill)+
  theme(plot.margin=unit(c(1,1,1,2), 'cm'),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = ggtext::element_markdown(angle=45,hjust=1,size=18),
        axis.text.y = element_text(size=20),
        axis.title = element_text(size = 20),
        strip.text.x = element_text(size=20),
        plot.title = element_text(size = 27),
        legend.text = element_text(size = 20),
        legend.title = element_text(size = 25),
        legend.position = "none")+
  ggtitle(paste("Alpha diversity of gene families in",as.character(host.class[host]),gg.title.groups))

div.plot.with.dots<-div.plot+
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.8)
# 
# for(image.format in image.formats){
#   ggsave(paste0("./images/diversity/alpha/",
#                 paste(paste(format(Sys.time(),format="%Y%m%d"),
#                             format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                       "alpha",paste(plot.metrics,collapse = "-"),host,
#                       comparison,agglom.rank,
#                       sep = "-"),".",image.format),
#          plot=div.plot,
#          width = 5000,height = 3000,
#          units = "px",dpi=300,device = image.format)
#   ggsave(paste0("./images/diversity/alpha/",
#                 paste(paste(format(Sys.time(),format="%Y%m%d"),
#                             format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
#                       "alpha-dots",paste(plot.metrics,collapse = "-"),host,
#                       comparison,agglom.rank,
#                       sep = "-"),".",image.format),
#          plot=div.plot.with.dots,
#          width = 5000,height = 3000,
#          units = "px",dpi=300,device = image.format)
# }


