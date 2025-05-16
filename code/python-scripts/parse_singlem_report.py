import pandas as pd
import re
import argparse
import sys
#TODO:change comments
# bracken_krona_txt_dir=output/kraken2_pipeline/bracken_krona_txt
# date_time=20240409_17_32_40
# for FILE in ${bracken_krona_txt_dir}/${date_time}_*_trim_no_minimizer_data.b.krona.txt
# do 
#  SAMPLE=$(echo ${FILE} | sed "s/\.b.krona.txt//"|sed "s/${date_time}_//"|sed "s/_no_minimizer_data//")
#  base_name=$(basename "$SAMPLE" )
#  python code/python-scripts/parse_singlem_report.py \
#  	--input_file ${bracken_krona_txt_dir}/${date_time}_${base_name}_no_minimizer_data.b.krona.txt \
# 	--output_file ${bracken_krona_txt_dir}/${date_time}_${base_name}_no_minimizer_data_table.tsv;
# done  
# File names have this format:
# "./output/bracken_krona_txt/20240409_17_32_40_H4_wms.b.krona.txt"


# input_filename="output/singlem_pipeline/singlem_output/2D10_wms_profile.with_extras.tsv"
# with open(input_filename) as f:
#         content = f.read().splitlines()
# lines=content[1:]
# line=lines[3]
# orig_line=line.split("\t")
# taxa_line=orig_line[5]
# taxa_line=taxa_line.split("; ")
# # skip the "Root" string which is in all rows
# taxa_line.remove("Root")


def singlem_report_to_table(args):
    input_filename, output_filename=args
    input_filename=input_filename[0]
    output_filename=output_filename[0]
    # Initialize an empty dictionary to store taxonomic information
    taxonomy_data = {"Classification_level": [],
                    "Relative_abundance": [],
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
    # Skip the header
    for line in content[1:]:
        orig_line=line.split("\t")
        # the columns are: sample, coverage, full_coverage, relative_abundance, level, taxonomy
        taxa_line=orig_line[5]
        taxa_line=taxa_line.split("; ")
        # skip the "Root" string which is in all rows
        taxa_line.remove("Root")
        
        kingdom=re.compile("d__")
        kingdom=list(filter(kingdom.match, taxa_line))
        
        phylum=re.compile("p__")
        phylum=list(filter(phylum.match, taxa_line))
        
        taxonomic_class=re.compile("c__")
        taxonomic_class=list(filter(taxonomic_class.match, taxa_line))
        
        order=re.compile("o__")
        order=list(filter(order.match, taxa_line))

        family=re.compile("f__")
        family=list(filter(family.match, taxa_line))

        genus=re.compile("g__")
        genus=list(filter(genus.match, taxa_line))

        species=re.compile("s__")
        species=list(filter(species.match, taxa_line))

        # this is the lowest level in the string
        taxonomy_data["Classification_level"].append(orig_line[4])   
        taxonomy_data["Relative_abundance"].append(orig_line[3])   
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
        print("Usage: python3 parse_singlem_report.py <input_file> <output_file>")
        sys.exit(1)
    parser = argparse.ArgumentParser(description='Process SingleM output in long format.')
    parser.add_argument('--input_file', nargs=1, required=True, help='the TXT file to process')
    parser.add_argument('--output_file', nargs=1, required=True, help='the output file name')
    args = parser.parse_args()
    singlem_report_to_table([args.input_file,args.output_file])
    
