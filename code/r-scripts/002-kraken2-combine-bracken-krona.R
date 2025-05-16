library(tidyverse)
# install.packages(
#   "microViz",
#   repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
# )
metadatadir<-paste0("../amplicon_nmr/data/metadata/pooled-metadata/") # directory with metadata

date_time="20240409_17_32_40_"
bracken_krona_txt.dir<-"output/kraken2_pipeline/bracken_krona_txt"
bracken_krona_txt.files<-
  list.files(bracken_krona_txt.dir)[grep("trim_no_minimizer_data_table.tsv",list.files(bracken_krona_txt.dir))]


list.tables<-list()


for(filename in bracken_krona_txt.files){
  sample.base.name<-gsub("_trim_no_minimizer_data_table.tsv","",filename)
  sample.base.name<-gsub(date_time,"",sample.base.name)
  bracken_krona_table<-
    read.table(file.path(bracken_krona_txt.dir,filename), 
               header = T, 
               sep = "\t",
               fill=TRUE,
               quote = "", 
               comment.char = "@",
               na.strings = c("","NA"))
  bracken_krona_table<-bracken_krona_table%>%
    filter(Abundance!=0)%>%
    relocate(Kingdom:Species)
  colnames(bracken_krona_table)[which(colnames(bracken_krona_table)=="Abundance")]<-sample.base.name
  
  list.tables[[sample.base.name]]<-bracken_krona_table
  
}

red.merged<-list.tables %>%
  Reduce(function(dtf1,dtf2) full_join(dtf1,
                                       dtf2,
                                       by=c("Kingdom","Phylum",
                                            "Class","Order",
                                            "Family","Genus","Species")), .)
# red.merged<-red.merged%>%
#   replace(is.na(.),0)

write.table(red.merged,
            file=file.path("./output/rtables",
                           paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "combined_report.tsv",sep = "_")),
            sep ="\t",
            quote = F,
            row.names = F)


### Singlem
date_time="20240928_17_00_00_"
bracken_krona_txt.dir<-"output/singlem_pipeline/singlem_output/singlem_profiles"
bracken_krona_txt.files<-
  list.files(bracken_krona_txt.dir)[grep("wms_profile_clean.tsv",list.files(bracken_krona_txt.dir))]

tax_levels=c("domain"="Kingdom",
             "phylum"="Phylum",
             "class"="Class",
             "order"="Order",
             "family"="Family",
             "genus"="Genus",
             "species"="Species")
list.tables<-list()

for(filename in bracken_krona_txt.files){
  sample.base.name<-gsub("_wms_profile_clean.tsv","",filename)
  sample.base.name<-gsub(date_time,"",sample.base.name)
  bracken_krona_table<-
    read.table(file=file.path(bracken_krona_txt.dir,filename),
               sep="\t",
               header = T)
  bracken_krona_table=subset(bracken_krona_table,Classification_level!="root")
  new_bracken_krona_table=subset(bracken_krona_table,Classification_level=="species")
  for (tax_rank in rev(tax_levels)[-1]){
    # current_classif_level is stored in the column "Classification_level"
    # tax_rank is column names or tax_levels (they're same)
    current_classif_level=names(tax_levels[tax_levels==tax_rank])
    current_level_position=which(tax_rank == colnames(bracken_krona_table))
    temp_bracken_krona_table=subset(bracken_krona_table,Classification_level==current_classif_level)
    
    for (taxon_name in temp_bracken_krona_table[,tax_rank]){
      if (taxon_name%in%new_bracken_krona_table[,tax_rank]){
        # If, for example, a genus is found in the dataframe with species names
        # (search the genus column), we want to know how much % in the genus is unclassified
        
        # new_bracken_krona_table_taxon_rows are the rows with species names that
        # belong to a genus.
        new_bracken_krona_table_taxon_rows=new_bracken_krona_table[which(new_bracken_krona_table[[tax_rank]]==taxon_name),]
        # Sum their abundances
        sum_abund_taxon_rows=sum(new_bracken_krona_table_taxon_rows$Relative_abundance)
        # Get the relative abundance of the whole genus from the higher rank bracken_krona_table
        taxon_name_total_rel_abundance=temp_bracken_krona_table[which(temp_bracken_krona_table[[tax_rank]]==taxon_name),"Relative_abundance"]
        # Find the difference between the sum of individual species abundance and 
        # the whole abundance of the genus
        unclassified_abundance=taxon_name_total_rel_abundance-sum_abund_taxon_rows
        if (unclassified_abundance>0){
          # If the difference is not 0, some reads are unclassified.
          # So, we create a row with unclassified %
          newrow=temp_bracken_krona_table[which(temp_bracken_krona_table[[tax_rank]]==taxon_name),]
          
          # Unclassified is the lower rank. so, we put the name of the genus into the 
          # species column and add 'Genus' to show it's unclassified. It'd look 
          # like "g__Mogibacterium Genus".
          newrow[current_level_position+1]=paste(taxon_name, tax_levels[tax_levels==tax_rank])
          newrow$Relative_abundance=unclassified_abundance
          new_bracken_krona_table=rbind(new_bracken_krona_table,newrow)
        }
      }else{
        # If there are no entries in the lower rank (no species that belong to a genus),
        # we add an "unclassified" row
        newrow=temp_bracken_krona_table[which(temp_bracken_krona_table[[tax_rank]]==taxon_name),]
        
        # Unclassified is the lower rank. so, we put the name of the genus into the 
        # species column and add 'Genus' to show it's unclassified. It'd look 
        # like "g__Mogibacterium Genus".
        newrow[current_level_position+1]=paste(taxon_name, tax_levels[tax_levels==tax_rank])
        newrow$Relative_abundance=temp_bracken_krona_table[which(temp_bracken_krona_table[[tax_rank]]==taxon_name),"Relative_abundance"]
        new_bracken_krona_table=rbind(new_bracken_krona_table,newrow)
      }
    }
  }
  new_bracken_krona_table<-new_bracken_krona_table%>%
    mutate(Phylum=ifelse(Phylum=="",Kingdom,Phylum),
           Class=ifelse(Class=="",Phylum,Class),
           Order=ifelse(Order=="",Class,Order),
           Family=ifelse(Family=="",Order,Family),
           Genus=ifelse(Genus=="",Family,Genus),
           Species=ifelse(Species=="",Genus,Species))
  stopifnot(sum(new_bracken_krona_table$Relative_abundance)==100)
  print(sum(new_bracken_krona_table$Relative_abundance))
  print(table(subset(new_bracken_krona_table,Classification_level=="species")[,"Species"]%in%bracken_krona_table$Species))
  print(table(subset(new_bracken_krona_table,Classification_level=="genus")[,"Genus"]%in%bracken_krona_table$Genus))
  print(table(subset(new_bracken_krona_table,Classification_level=="family")[,"Family"]%in%bracken_krona_table$Family))
  print(table(subset(new_bracken_krona_table,Classification_level=="order")[,"Order"]%in%bracken_krona_table$Order))
  print(table(subset(new_bracken_krona_table,Classification_level=="class")[,"Class"]%in%bracken_krona_table$Class))
  print(table(subset(new_bracken_krona_table,Classification_level=="phylum")[,"Phylum"]%in%bracken_krona_table$Phylum))
  rm(new_bracken_krona_table_taxon_rows)
  rm(temp_bracken_krona_table)
  new_bracken_krona_table<-new_bracken_krona_table%>%
    select(-Classification_level)
  colnames(new_bracken_krona_table)[which(colnames(new_bracken_krona_table)=="Relative_abundance")]<-sample.base.name
  list.tables[[sample.base.name]]<-new_bracken_krona_table
}
red.merged<-list.tables %>%
  Reduce(function(dtf1,dtf2) full_join(dtf1,
                                       dtf2,
                                       by=c("Kingdom","Phylum",
                                            "Class","Order",
                                            "Family","Genus","Species")), .)
red.merged<-red.merged%>%
  replace(is.na(.),0)%>%
  relocate(Kingdom:Species)
table(colSums(red.merged[,8:18])==100)
write.table(red.merged,
            file=file.path("./output/rtables",
                           paste(paste(format(Sys.time(),format="%Y%m%d"),
                                       format(Sys.time(),format = "%H_%M_%S"),sep = "_"),
                                 "singlem_combined_profile.tsv",sep = "_")),
            sep ="\t",
            quote = F,
            row.names = F)
