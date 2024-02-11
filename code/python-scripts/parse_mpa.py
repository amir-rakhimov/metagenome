import pandas as pd
import re

# Read the tab-separated table into a DataFrame
file_path = "output/kraken2_pipeline/kraken2_reports/20240122_combined_mpa.tsv"  # Replace with your actual file path
file_out_path="output/kraken2_pipeline/kraken2_reports/20240122_combined_mpa_clean.tsv"

taxonomy_data = {"Kingdom": [],
                 "Phylum": [],
                 "Class": [], 
                 "Order": [], 
                 "Family": [], 
                 "Genus": [],
                 "Species": []}
with open(file_path) as f:
    content = f.read().splitlines()


for line in content:
    line=line.split("\t")
    if line[0]=="#Classification":
        continue
    line=line[0]
    line=line.split("|")

    kingdom=re.compile("k__")
    kingdom=list(filter(kingdom.match, line))
    
    phylum=re.compile("p__")
    phylum=list(filter(phylum.match, line))
    
    taxonomic_class=re.compile("c__")
    taxonomic_class=list(filter(taxonomic_class.match, line))
    
    order=re.compile("o__")
    order=list(filter(order.match, line))

    family=re.compile("f__")
    family=list(filter(family.match, line))

    genus=re.compile("g__")
    genus=list(filter(genus.match, line))

    species=re.compile("s__")
    species=list(filter(species.match, line))

    taxonomy_data["Kingdom"].append(''.join(kingdom))
    taxonomy_data["Phylum"].append(''.join(phylum))
    taxonomy_data["Class"].append(''.join(taxonomic_class))
    taxonomy_data["Order"].append(''.join(order))
    taxonomy_data["Family"].append(''.join(family))
    taxonomy_data["Genus"].append(''.join(genus))
    taxonomy_data["Species"].append(''.join(species))

taxonomy_split = pd.DataFrame(taxonomy_data)

abundance_df = pd.read_csv(file_path, sep='\t')
abundance_df=abundance_df.drop(['#Classification'],axis=1)
# frames=[abundance_df,taxonomy_split]
# df=pd.concat(frames)
df=abundance_df.join(taxonomy_split)

df.to_csv(file_out_path, index=False,sep="\t")
