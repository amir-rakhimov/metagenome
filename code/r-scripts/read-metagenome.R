library(tidyverse)
breport.colnames<-c("Abundance","Kingdom","Phylum","Class","Order","Family",
                    "Genus","Species")
# h4.breport.filename<-"./output/bracken_krona_txt/20240122_H4_wms.b.krona.txt"
h4.breport.filename<-"./output/H4_wms_table.tsv"

h4.breport<-read.table(h4.breport.filename, 
                       header = T, 
                       sep = "\t",
                       fill=TRUE,
                       quote = "")

h3.breport.filename<-"./output/H3_wms_table.tsv"

h3.breport<-read.table(h3.breport.filename, 
                       header = T, 
                       sep = "\t",
                       fill=TRUE,
                       quote = "")


nrow(h4.breport)
nrow(h3.breport)

h4.breport%>%
  group_by(Phylum)%>%
  ggplot(aes(x=Phylum,y=Abundance))+
  geom_bar(stat='identity')+
  coord_flip()
sum(h4.breport$Abundance)

h4.phyla<-h4.breport%>%
  select(Abundance,Phylum)%>%
  group_by(Phylum)%>%
  summarise(Abundance=sum(Abundance))
h3.phyla<-h3.breport%>%
  select(Abundance,Phylum)%>%
  group_by(Phylum)%>%
  summarise(Abundance=sum(Abundance))

df<-h3.phyla%>%
  full_join(h4.phyla,by="Phylum",suffix = c("H3","H4"),)%>%
  replace(is.na(.),0)%>%
  column_to_rownames("Phylum")


