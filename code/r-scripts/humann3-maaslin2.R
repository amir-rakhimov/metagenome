library(tidyverse)
library(Maaslin2)
merged.date<-"20240217"
pathabundance.filename<-paste0("output/humann3_pipeline/humann_merged_",
                                 merged.date,"/",merged.date,
                                 "_humann_out_pathabundance-cpm.tsv")
metadatadir<-paste0("../amplicon_nmr/data/metadata/pooled-metadata/") # directory with metadata

comparison<-"age"
host<-"NMR"
ref.level<-"between2and10" # choose the reference level
# this is for file names
if(host=="NMR"){
  host.labels<-c("NMR" = "*Heterocephalus glaber*")
}else{
  host.labels<-
    c("B6mouse" = "B6 mouse",
      "MSMmouse" = "MSM/Ms mouse",
      "FVBNmouse" = "FVB/N mouse")
}
custom.levels<-c("NMR"
)


pathabundance<-read.table(pathabundance.filename, 
                            header = T, 
                            sep = "\t",
                            fill=TRUE,
                            quote = "", 
                            comment.char = "@",
                            na.strings = c("","NA"))
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

pathabundance<-pathabundance[!grepl("UNMAPPED",pathabundance$Pathway),]
pathabundance<-pathabundance[!grepl("UNINTEGRATED",pathabundance$Pathway),]
pathabundance<-pathabundance%>%
  pivot_longer(!Pathway,
               names_to = "Sample")
pathabundance<-within(pathabundance, 
            Pathway<-data.frame(do.call('rbind', 
                                             strsplit(as.character(Pathway), 
                                                      ':', fixed=TRUE))))
pathabundance_full<-pathabundance
pathabundance_full$full_pathway<-pathabundance$Pathway$X2
pathabundance_full$Pathway<-pathabundance$Pathway$X1
pathabundance_full<-within(pathabundance_full, 
                      full_pathway<-data.frame(do.call('rbind', 
                                                  strsplit(as.character(full_pathway), 
                                                           '|', fixed=TRUE))))
pathabundance_full$full_pathway<-pathabundance_full$full_pathway$X1

pathabundance$Pathway<-pathabundance$Pathway$X1
pathabundance<-pathabundance%>%
  group_by(Sample,Pathway)%>%
  summarise(Abundance=sum(value))

pathabundance.wide<-pathabundance%>%
  pivot_wider(names_from = Sample,
              values_from = Abundance)

# pathabundance.clr<-clr(pathabundance.wide[,-1])

pathabundance.wide<-pathabundance.wide%>%column_to_rownames("Pathway")



#################
gene_family.filename<-paste0("output/humann3_pipeline/humann_merged_",
                               merged.date,"/",merged.date,
                               "_humann_out_genefamilies-cpm-rxn-named.tsv")
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

gene_family<-gene_family[!grepl("UNMAPPED",gene_family$Gene.Family),]
gene_family<-gene_family[!grepl("UNGROUPED",gene_family$Gene.Family),]

gene_family<-gene_family%>%
  pivot_longer(!Gene.Family,
               names_to = "Sample")
gene_family<-within(gene_family, 
                    Gene.Family<-data.frame(do.call('rbind', 
                                                  strsplit(as.character(Gene.Family), 
                                                           ':', fixed=TRUE))))
gene_family_full<-gene_family
gene_family_full$full_gene_family<-gene_family$Gene.Family$X2
gene_family_full$Gene.Family<-gene_family_full$Gene.Family$X1
gene_family_full<-within(gene_family_full, 
                         full_gene_family<-data.frame(do.call('rbind', 
                                                            strsplit(as.character(full_gene_family), 
                                                                     '|', fixed=TRUE))))
gene_family_full$full_gene_family<-gene_family_full$full_gene_family$X1


gene_family$Gene.Family<-gene_family$Gene.Family$X1
gene_family<-gene_family%>%
  group_by(Sample,Gene.Family)%>%
  summarise(Abundance=sum(value))

gene_family.wide<-gene_family%>%
  pivot_wider(names_from = Sample,
              values_from = Abundance)

# gene_family.clr<-clr(gene_family.wide[,-1])

gene_family.wide<-gene_family.wide%>%column_to_rownames("Gene.Family")














# Create metadata
metadata.filename<-paste0(metadatadir,
                          paste("filenames-single-pooled-raw-supercomp.tsv", sep = "-"))
custom.md<-read.table(metadata.filename, header = T)
colnames(custom.md)[1]<-"Sample" # set the first column name as Sample
custom.md<-custom.md%>%
  filter(Sample%in%gene_family$Sample)%>%
  column_to_rownames(var = "Sample")%>%
  select(-absolute.filepath)
custom.md$class<-as.factor(custom.md$class)
custom.md$animal<-as.factor(custom.md$animal)
custom.md$sex<-as.factor(custom.md$sex)
custom.md$birthday<-as.Date(custom.md$birthday)
custom.md$age<-year(as.period(interval(custom.md$birthday,now())))


age.breaks<-c(min(custom.md$age),10,max(custom.md$age))
custom.md<-custom.md%>%
  mutate(age_group=cut(age, breaks = age.breaks, 
                       include.lowest = TRUE))%>%
  mutate(age_group=as.factor(age_group))
substr(gsub("\\[|\\]","",as.character(custom.md$age_group)[1]),2,2)

unique_levels <- custom.md %>%
  distinct(age_group) %>%
  arrange(age_group) %>%
  mutate(new_age_group=gsub("\\[|\\]|\\(","",age_group))%>%
  separate(new_age_group,c("leftbound","rightbound"))%>%
  mutate(new_age_group = paste0("between", leftbound,"and",rightbound))%>%
  select(-leftbound,-rightbound)
custom.md$Sample<-rownames(custom.md)

custom.md <- custom.md %>%
  left_join(unique_levels, by = "age_group")
colnames(custom.md)[which(colnames(custom.md)=="age_group")]<-"old_age_group"
colnames(custom.md)[which(colnames(custom.md)=="new_age_group")]<-"age_group"
rownames(custom.md)<-custom.md$Sample

maaslin.reference<-paste("age_group",ref.level,sep = ",")
maaslin.comparison<-"age_group"

set.seed(1)
maaslin.fit_data =
  Maaslin2(input_data = pathabundance.wide,
           input_metadata = custom.md,
           min_prevalence = 0,
           normalization = "NONE",
           output = paste0("./output/maaslin2/","humann3","-output/",
                           paste("pathway",host,
                                                 comparison,
                                                 ref.level,sep="-")),
           fixed_effects = maaslin.comparison,
           reference = maaslin.reference,
           max_significance = 0.5)
save.image(paste0("./output/rdafiles/",
                  paste("maaslin-humann3-pathabundance",host,
                        comparison,ref.level,
                        "workspace.RData",sep="-")))

set.seed(1)
maaslin.fit_data = 
  Maaslin2(input_data = gene_family.wide, 
           input_metadata = custom.md, 
           min_prevalence = 0,
           normalization = "NONE",
           output = paste0("./output/maaslin2/","humann3","-output/",
                           paste("gene_family",host,
                                 comparison,
                                 ref.level,sep="-")), 
           fixed_effects = maaslin.comparison,
           reference = maaslin.reference,
           max_significance = 0.5)
save.image(paste0("./output/rdafiles/",
                  paste("maaslin-humann3-gene_family",host,
                        comparison,ref.level,
                        "workspace.RData",sep="-")))


# Downstream analysis
load(paste0("./output/rdafiles/",
            paste("maaslin-humann3-gene_family",host,
                  comparison,ref.level,
                  "workspace.RData",sep="-")))
source("../amplicon_nmr/code/r-scripts/make_features_maaslin.R")

# maaslin.signif.features<-maaslin.fit_data$results%>%
#   filter(qval<0.05)
maaslin.signif.features<-maaslin.fit_data$results%>%
  filter(qval<0.05) # should be qval





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
  left_join(unique(gene_family_full[,c("Gene.Family","full_gene_family")]),by=c("feature"="Gene.Family"))


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
  select(feature,coef,assoc.str,metadata,value,full_gene_family)
# write.table(maaslin.features.list,
#             file = paste0("./output/rtables/","maaslin2-humann3-","signif-gene_family-",
#                           paste(custom.levels,collapse = '-'),".tsv"),
#             row.names = F,
#             col.names = F,
#             quote = F,
#             sep = "\t")
# 
# write.table(maaslin.features.list$feature,
#             file = paste0("./output/rtables/","maaslin2-humann3-","signif-gene_family-list-",
#                           paste(custom.levels,collapse = '-'),".tsv"),
#             row.names = F,
#             col.names = F,
#             quote = F,
#             sep = "\t")


tot.pathw<-maaslin.features.list%>%group_by(full_gene_family)%>%mutate(tot=n())%>%arrange(-tot)
# write.table(tot.pathw,
#             file = paste0("./output/rtables/","maaslin2-humann3-","signif-gene_family-tot.pathw-",
#                           paste(custom.levels,collapse = '-'),".tsv"),
#             row.names = F,
#             quote = F,
#             sep = "\t")
# Barplot for gene families
level.names<-names(sort(table(maaslin.features.list$full_gene_family),decreasing = F))
level.names<-gsub("^ \\(expasy\\)|^ \\(metacyc\\)","",level.names)
level.names<-gsub("^ ", "",level.names)
level.names<-str_to_sentence(level.names)


maaslin.gene_family.freq.plot<-maaslin.features.list%>%
  mutate(new_full_gene_family=gsub("^ \\(expasy\\)|^ \\(metacyc\\)","",full_gene_family))%>%
  mutate(new_full_gene_family=gsub("^ ","",new_full_gene_family))%>%
  mutate(new_full_gene_family=str_to_sentence(new_full_gene_family))%>%
  
  mutate(new_full_gene_family=factor(new_full_gene_family,levels=level.names))%>%
  ggplot(aes(y=new_full_gene_family))+
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
              paste(Sys.Date(),"barplot",
                    "humann3-maaslin2-gene_family.freq",sep = "-"),
              ".png"),
       plot = maaslin.gene_family.freq.plot,
       width = 6000,height = 3000,
       units = "px",dpi=300,device = "png")











load(paste0("./output/rdafiles/",
            paste("maaslin-humann3-pathabundance",host,
                  comparison,ref.level,
                  "workspace.RData",sep="-")))
source("../amplicon_nmr/code/r-scripts/make_features_maaslin.R")
maaslin.signif.features<-maaslin.fit_data$results%>%
  filter(qval<0.1) # should be qval
foo<-pathabundance
foo$maaslin<-foo$Pathway
foo<-make_features_maaslin(foo,"maaslin")
foo<-unique(foo[,c("maaslin","Pathway")])
maaslin.signif.features<-maaslin.signif.features%>%
  left_join(foo[,c("maaslin","Pathway")],by=c("feature"="maaslin"))%>%
  distinct()
rm(foo)
maaslin.signif.features$feature<-maaslin.signif.features$Pathway
maaslin.signif.features<-subset(maaslin.signif.features, select=-Pathway)


maaslin.signif.features<-maaslin.signif.features%>%
  left_join(unique(pathabundance_full[,c("Pathway","full_pathway")]),by=c("feature"="Pathway"))


table(maaslin.signif.features$feature%in%pathabundance$Pathway)

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

table(maaslin.signif.decreased$feature%in%pathabundance$Pathway)
table(maaslin.signif.increased$feature%in%pathabundance$Pathway)

maaslin.signif.decreased$feature[!maaslin.signif.decreased$feature%in%pathabundance$Pathway]
maaslin.signif.increased$feature[!maaslin.signif.increased$feature%in%pathabundance$Pathway]
maaslin.features.list<-maaslin.signif.increased%>%
  bind_rows(maaslin.signif.decreased)%>%
  filter(abs(assoc.str)>1)


write.table(maaslin.features.list,
            file = paste0("./output/rtables/","maaslin2-humann3-","signif-pathabundance-",
                          paste(custom.levels,collapse = '-'),".tsv"),
            row.names = F,
            col.names = F,
            quote = F,
            sep = "\t")
write.table(maaslin.features.list$feature,
            file = paste0("./output/rtables/","maaslin2-humann3-","signif-pathabundance-list-",
                          paste(custom.levels,collapse = '-'),".tsv"),
            row.names = F,
            col.names = F,
            quote = F,
            sep = "\t")


################
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
gene_family_full$full_gene_family<-gene_family$Gene.Family$X2
gene_family_full$Gene.Family<-gene_family_full$Gene.Family$X1
gene_family_full<-within(gene_family_full, 
                         full_gene_family<-data.frame(do.call('rbind', 
                                                              strsplit(as.character(full_gene_family), 
                                                                       '|', fixed=TRUE))))



foo<-gene_family_full
gene_family_full$full_pathway<-gene_family_full$full_gene_family$X1
gene_family_full$strat<-gene_family_full$full_gene_family$X2
gene_family_full<-gene_family_full%>%
  select(-full_gene_family)


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
  arrange(Group)
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
  facet_grid(~age_group,
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
ggsave(paste0("./images/barplots/",paste(Sys.Date(),"humann3.plot-linolenoyl-rxn",
                                   host,comparison,ref.level,sep = "-"),
              ".png"),
       plot = humann3.plot,
       width = 6000,height = 3000,
       units = "px",dpi=300,device = "png")


ggsave(paste0("./images/barplots/",paste(Sys.Date(),"humann3.plot-linolenoyl-rxn",
                                         host,comparison,ref.level,sep = "-"),
              ".tiff"),
       plot = humann3.plot,
       width = 6000,height = 3000,
       units = "px",dpi=300,device = "tiff")
