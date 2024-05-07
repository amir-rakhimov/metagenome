import pandas as pd
import re
import argparse
import sys

# for FILE in ${bracken_krona_txt_dir}/${date_time}_*_trim_no_minimizer_data.b.krona.txt
# do 
#  SAMPLE=$(echo ${FILE} | sed "s/\.b.krona.txt//"|sed "s/${date_time}_//"|sed "s/_no_minimizer_data//")
#  base_name=$(basename "$SAMPLE" )
#  python code/python-scripts/parse-breport.py \
#  	--input_file ${bracken_krona_txt_dir}/${date_time}_${base_name}_no_minimizer_data.b.krona.txt \
# 	--output_file ${bracken_krona_txt_dir}/${date_time}_${base_name}_no_minimizer_data_table.tsv;
# done  
# File names have this format:
# "./output/bracken_krona_txt/20240409_17_32_40_H4_wms.b.krona.txt"
def breport_to_table(args):
    input_filename, output_filename=args
    input_filename=input_filename[0]
    output_filename=output_filename[0]
    # Initialize an empty dictionary to store taxonomic information
    taxonomy_data = {"Abundance": [],
                    "Kingdom": [],
                    "Phylum": [],
                    "Class": [], 
                    "Order": [], 
                    "Family": [], 
                    "Genus": [],
                    "Species": []}

    # sample_name="2D10"
    # date_time="20240409_17_32_40"
    # input_filename=filename_dir+date_time+"_"+sample_name+"_trim_no_minimizer_data.b.krona.txt"
    # output_filename=filename_dir+date_time+"_"+sample_name+"_trim_no_minimizer_data_table.tsv"
    # Open the file and read its content.
    with open(input_filename) as f:
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
    df.to_csv(output_filename, index=False,sep="\t")



if __name__ == "__main__":
    if len(sys.argv) <2:
        print("Usage: python3 parse-breport.py <input_file> <output_file>")
        sys.exit(1)
    parser = argparse.ArgumentParser(description='Process some Krona TXT files from Bracken.')
    parser.add_argument('--input_file', nargs=1, required=True, help='the TXT file to process')
    parser.add_argument('--output_file', nargs=1, required=True, help='the output file name')
    args = parser.parse_args()
    breport_to_table([args.input_file,args.output_file])
    
