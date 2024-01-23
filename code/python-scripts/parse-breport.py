import pandas as pd
import re

# "./output/bracken_krona_txt/20240122_H4_wms.b.krona.txt"
# Initialize an empty dictionary to store taxonomic information
taxonomy_data = {"Abundance": [],
                 "Kingdom": [],
                 "Phylum": [],
                 "Class": [], 
                 "Order": [], 
                 "Family": [], 
                 "Genus": [],
                 "Species": []}

filename="./output/bracken_krona_txt/20240122_H4_wms.b.krona.txt"
# Open the file and read its content.
with open(filename) as f:
    content = f.read().splitlines()

# Display the file's content line by line.
for line in content:
    line=line.split("\t")
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

    taxonomy_data["Abundance"].append(line[0])   
    taxonomy_data["Kingdom"].append(''.join(kingdom))
    taxonomy_data["Phylum"].append(''.join(phylum))
    taxonomy_data["Class"].append(''.join(taxonomic_class))
    taxonomy_data["Order"].append(''.join(order))
    taxonomy_data["Family"].append(''.join(family))
    taxonomy_data["Genus"].append(''.join(genus))
    taxonomy_data["Species"].append(''.join(species))


# Create a DataFrame from the dictionary
df = pd.DataFrame(taxonomy_data)



# Save the DataFrame to a new CSV file with missing values as NA
df.to_csv("H4_wms_table.tsv", index=False,sep="\t")

