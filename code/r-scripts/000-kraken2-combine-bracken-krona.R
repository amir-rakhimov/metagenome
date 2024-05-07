library(tidyverse)
# install.packages(
#   "microViz",
#   repos = c(davidbarnett = "https://david-barnett.r-universe.dev", getOption("repos"))
# )
agglom.rank<-"Species"
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
