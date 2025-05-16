library(tidyverse)
fname="./output/singlem_pipeline/singlem_output/singlem_profiles/20240928_17_00_00_2D10_wms_profile_clean.tsv"
df<-read.table(file=fname,sep="\t",header = T)
df=subset(df,Classification_level!="root")
tax_levels=c("domain"="Kingdom",
             "phylum"="Phylum",
             "class"="Class",
             "order"="Order",
             "family"="Family",
             "genus"="Genus",
             "species"="Species")

new_df=subset(df,Classification_level=="species")
# temp_df=subset(df,Classification_level=="genus")
# for (taxonomic_level in tax_levels){
#   for (genus_name in temp_df$Genus){
#     if (genus_name%in%new_df$Genus){
#       new_df_genus_rows=subset(new_df,Genus==genus_name)
#       sum_abund_genus_rows=sum(new_df_genus_rows$Relative_abundance)
#       genus_name_total_rel_abundance=subset(temp_df,Genus==genus_name)[,"Relative_abundance"]
#       unclassified_abundance=genus_name_total_rel_abundance-sum_abund_genus_rows
#       if (unclassified_abundance>0){
#         newrow=subset(temp_df,Genus==genus_name)
#         newrow$Species=paste("Unclassified",genus_name,sep="_")
#         newrow$Relative_abundance=unclassified_abundance
#         new_df=rbind(new_df,newrow)
#       }
#     }
#   }
# }


for (tax_rank in rev(tax_levels)[-1]){
  # current_classif_level is stored in the column "Classification_level"
  # tax_rank is column names or tax_levels (they're same)
  current_classif_level=names(tax_levels[tax_levels==tax_rank])
  current_level_position=which(tax_rank == colnames(df))
  temp_df=subset(df,Classification_level==current_classif_level)
  
  for (taxon_name in temp_df[,tax_rank]){
    if (taxon_name%in%new_df[,tax_rank]){
      # If, for example, a genus is found in the dataframe with species names
      # (search the genus column), we want to know how much % in the genus is unclassified
      
      # new_df_taxon_rows are the rows with species names that
      # belong to a genus.
      new_df_taxon_rows=new_df[which(new_df[[tax_rank]]==taxon_name),]
      # Sum their abundances
      sum_abund_taxon_rows=sum(new_df_taxon_rows$Relative_abundance)
      # Get the relative abundance of the whole genus from the higher rank df
      taxon_name_total_rel_abundance=temp_df[which(temp_df[[tax_rank]]==taxon_name),"Relative_abundance"]
      # Find the difference between the sum of individual species abundance and 
      # the whole abundance of the genus
      unclassified_abundance=taxon_name_total_rel_abundance-sum_abund_taxon_rows
      if (unclassified_abundance>0){
        # If the difference is not 0, some reads are unclassified.
        # So, we create a row with unclassified %
        newrow=temp_df[which(temp_df[[tax_rank]]==taxon_name),]
        
        # Unclassified is the lower rank. so, we put the 'unclassified' into the 
        # species column.
        newrow[current_level_position+1]=paste(taxon_name, tax_levels[tax_levels==tax_rank])
        newrow$Relative_abundance=unclassified_abundance
        new_df=rbind(new_df,newrow)
      }
    }else{
      # If there are no entries in the lower rank (no species that belong to a genus),
      # we add an "unclassified" row
      newrow=temp_df[which(temp_df[[tax_rank]]==taxon_name),]
      
      # Unclassified is the lower rank. so, we put the 'unclassified' into the 
      # species column.
      newrow[current_level_position+1]=paste(taxon_name,tax_levels[tax_levels==tax_rank])
      newrow$Relative_abundance=temp_df[which(temp_df[[tax_rank]]==taxon_name),"Relative_abundance"]
      new_df=rbind(new_df,newrow)
    }
  }
}

sum(new_df$Relative_abundance)
table(subset(new_df,Classification_level=="species")[,"Species"]%in%df$Species)
table(subset(new_df,Classification_level=="genus")[,"Genus"]%in%df$Genus)
table(subset(new_df,Classification_level=="family")[,"Family"]%in%df$Family)
table(subset(new_df,Classification_level=="order")[,"Order"]%in%df$Order)
table(subset(new_df,Classification_level=="class")[,"Class"]%in%df$Class)
table(subset(new_df,Classification_level=="phylum")[,"Phylum"]%in%df$Phylum)

# fill NA columns
new_df<-new_df%>%
  mutate(Phylum=ifelse(Phylum=="",Kingdom,Phylum),
         Class=ifelse(Class=="",Phylum,Class),
         Order=ifelse(Order=="",Class,Order),
         Family=ifelse(Family=="",Order,Family),
         Genus=ifelse(Genus=="",Family,Genus),
         Species=ifelse(Species=="",Genus,Species))
new_df%>%
  mutate(Sample="2D10")%>%
  mutate(Species=ifelse(Relative_abundance<1,"Remainder",Species))%>%
  ggplot(aes(x=Sample,y=Relative_abundance,fill=Species))+
  geom_bar(stat = "identity")
