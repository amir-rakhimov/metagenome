library(tidyverse)
dbcan.output.dir<-"./output/mag_assembly/dbcan_output"
dbcan.date_time<-"20250401_16_20_06"
# Step 1: Load and inspect the data
dbcan.output.file_list<-list.files(dbcan.output.dir,pattern = "overview.txt$", 
                                   recursive = TRUE)

# Helper function to process one file
process_overview <- function(fname,file_dir=NULL,file.date_time=NULL, ...) {
  sample_name <- str_replace(fname,paste0(file.date_time,"_"), "")%>%
    str_replace("_dbcan_proteins\\/overview.txt", "")
  # Read overview file
  fname.fullpath<-file.path(file_dir,paste(file.date_time,sample_name,
                                           "dbcan_proteins/overview.txt",sep="_"))
  df <- read.delim(fname.fullpath, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  df<-df%>%
    rename("Num_of_tools"="X.ofTools")

  # Filter for confident CAZyme hits (≥2 tools agree): count how many tools
  # gave a hit per gene
  df.filtered <- df %>%
    filter(Num_of_tools >= 2)%>%
    mutate(Sample = sample_name)
  
  # Pivot to long format, drop empty entries
  df.long <- df.filtered %>%
    pivot_longer(cols = -c("Gene.ID","EC.","Num_of_tools","Sample"), 
                 names_to = "Tool",
                 values_to = "Family") %>%
    filter(Family != "-")%>%
    distinct(`Gene.ID`, Family,Sample)  # Avoid duplicates
  
    
  
  # Split multi-domain strings (e.g. "GH5(1-100)+CBM50(110-180)")
  df.split <- df.long %>%
    separate_rows(Family, sep = "\\+") %>%
    mutate(Family = str_extract(Family, "^[^\\(]+"))  # remove coordinates
  
  return(df.split)
}

# Apply to all files
dbcan_data.all <-map_dfr(dbcan.output.file_list, process_overview,
        file_dir = dbcan.output.dir, 
        file.date_time = dbcan.date_time)

# Preview
head(dbcan_data.all)



# Summarize counts: how many times a family is observed in a sample
summary_counts <- dbcan_data.all %>%
  count(Sample, Family, sort = TRUE)

# View top 10 families per sample
summary_counts %>%
  group_by(Sample) %>%
  slice_max(n, n = 10) %>%
  ungroup()

# Plot top families across all samples
top_families <- summary_counts %>%
  group_by(Family) %>%
  summarise(Total = sum(n)) %>%
  slice_max(Total, n = 10)

ggplot(summary_counts %>% filter(Family %in% top_families$Family),
       aes(x = reorder(Family, n), y = n, fill = Sample)) +
  geom_col(position = "dodge") +
  coord_flip() +
  labs(title = "Top CAZy Families Across Samples",
       x = "CAZy Family", y = "Count") +
  theme_minimal()

# Count the CAZymes in each group (GT, GH, etc)
dbcan_data.all%>%
  mutate(Family=gsub("_[0-9]+$|_e[0-9]+$","", Family))%>%
  distinct(Family)%>%View

# TODO: Average relative abundance of top 10 CAZymes



# TODO: heatmap of cazymes: skewed distribution (one Family is too high so nothing looks ok)
heat.mat<-summary_counts%>%
  filter(Family%in%top_families$Family)%>%
  pivot_wider(names_from = "Family",values_from = "n")%>%
  column_to_rownames("Sample")%>%
  as.matrix.data.frame()
heatmap(heat.mat)



# dbcan.output.2d10.path<-file.path(dbcan.output.dir,"20250401_16_20_06_2D10_dbcan_proteins/overview.txt")
# dbcan.output.2d10<-read.delim(dbcan.output.2d10.path,header=TRUE,
#                               sep = "\t",stringsAsFactors = FALSE)
# dbcan.output.2d10<-dbcan.output.2d10%>%
#   rename("Num_of_tools"="X.ofTools")
# head(dbcan.output.2d10)
# 
# # Step 2: Filter for confident CAZyme hits (≥2 tools agree)
# # Count how many tools gave a hit per gene
# dbcan.output.2d10.filtered <- dbcan.output.2d10 %>%
#   filter(Num_of_tools >= 2)
# 
# 
# # Step 3: Summarize and visualize CAZyme families
# # Let’s count which CAZy families were most commonly found:
# # Combine family assignments into one column
# cazy_families <- dbcan.output.2d10.filtered %>%
#   pivot_longer(cols = -c("Gene.ID","EC.","Num_of_tools"), names_to = "Tool", values_to = "Family") %>%
#   filter(Family != "-") %>%
#   distinct(`Gene.ID`, Family)  # Avoid duplicates
# 
# # TODO: split GH43_12(57-354)+CBM91(387-564)
# # Count occurrences of each CAZy family
# family_counts <- cazy_families %>%
#   count(Family, sort = TRUE)
# 
# # View top 10
# head(family_counts, 10)
# 
# # Step 4: Plot CAZy family abundance
# # Plot top 20 CAZy families
# top_families <- family_counts %>%
#   slice_max(n, n = 20)
# 
# ggplot(top_families, aes(x = reorder(Family, n), y = n)) +
#   geom_col(fill = "#2E86AB") +
#   coord_flip() +
#   labs(title = "Top CAZy Families in Sample S1",
#        x = "CAZy Family", y = "Count") +
#   theme_bw()#+coord_cartesian(expand = FALSE)


dbcan_data.all<-dbcan_data.all%>%
  mutate(contig_id=gsub("_[0-9]+$","",Gene.ID))%>%
  mutate(contig_id=paste(Sample,contig_id,sep = "_"))
dbcan_data.contigs<-dbcan_data.all%>%
  distinct(contig_id)%>%
  pull

head(combined.mag_stats)

combined.mag_stats$contig_id_old<-str_replace(combined.mag_stats$contig_id,
                                              paste0("_",combined.mag_stats$Bin,"_k"),
                                              "_k")

combined.mag_stats.contigs<-combined.mag_stats%>%
  distinct(contig_id_old)%>%
  pull

length(dbcan_data.contigs)
length(combined.mag_stats.contigs)

setdiff(combined.mag_stats.contigs,dbcan_data.contigs)
setdiff(dbcan_data.contigs,combined.mag_stats.contigs)
intersect(dbcan_data.contigs,combined.mag_stats.contigs)
