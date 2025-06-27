"""
Script to extract 'Processed contigs' download links from a list of assembly JSON URLs 
and generate a TSV file mapping sample IDs to their corresponding contig URLs.

Usage:
    python script.py sample_list.txt -n output.tsv

Arguments:
    sample_list.txt : Text file containing a list of assembly JSON URLs.
    -n, --name      : Name of the output TSV file.
"""


# URL for the JSON file download
    url = str(asmbly)

    # Download the file
    response = requests.get(url)

    #Check if the download was succesful
    if response.status_code == 200:
        #Save file locally
        json_file = "tmp.jason"
        with open(json_file, 'w', encoding= 'utf-8') as file:
            file.write(response.text)
        #print(f"JSON file downloaded and saved as {json_file}")


    # open the file and load its contents
    with open(json_file, 'r', encoding = 'utf-8') as file:
        data = json.load(file)

    #Display the top-level keys in the JSON data
    #print(f"Top-level keys in the JSON: {list(data.keys())}")

    #Loop through the 'data' array
    if 'data' in data:
        data_array = data['data']
        for item in data_array:
            if 'attributes' in  item and 'description' in item['attributes']:
                key = item['attributes']['description'].get('label', [])
                if key == "Processed contigs":
                    #print(f"{key} found")
                    if 'links' in  item:
                        link = item['links'].get('self', [])
                        print(f"File found in: {link}")
        
                        return link
                else :
                    print(f"Link not found")
                    return None
            else :
                print(f"Key not found")
                return None
        else :
            print(f"Attributes not found")
            return None
    else :
            print(f"Data not found")
            return None

#Turns text file into a list
def lines_to_list(filename):
    """
    Reads a text file and creates a list where each element is a line from the file.

    Parameters:
        filename (str): The path to the text file.

    Returns:
        list: A list of strings, each representing a line from the file.
    """
    try:
        with open(filename, 'r') as file:
            lines = [line.strip() for line in file.readlines()]
        return lines
    except FileNotFoundError:
        print(f"Error: The file {filename} was not found.")
        return []
    except Exception as e:
        print(f"An error occurred: {e}")
        return []

#Function takes as input a txt file with assemblies and outputs a list with download urls 
def find_links(jlist):
    print("Retrieving URLs")
    urlist = []
    index = 0
    with open(jlist, 'r') as file:
        samples = [line.strip() for line in file]
        for sample in samples:
            url = json_to_url(sample)
            urlist.append(url)
            index = index + 1
            print(f"Retrieved {index} URL ")
    return urlist


#Function takes as inputs two lists and creates a .tsv file.
def create_table_and_save_tsv(ids, contigs):
    """
    Creates a table with IDs in the first column and corresponding values in the second column,
    and saves it as a TSV file without column names.

    Parameters:
        ids (list): A list of IDs.
        values (list): A list of values corresponding to the IDs.

    Returns:
        None
    """
    if len(ids) != len(contigs):
        raise ValueError("The lists 'ids' and 'values' must have the same length.")
    
    # Create a DataFrame with the provided lists
    table = pd.DataFrame({
        'ID': ids,
        'URL': contigs
    })

    return(table)

#MAIN
import argparse
import json
import requests
import os
import pandas as pd
from urllib.parse import urlparse

parser = argparse.ArgumentParser(description = "A script that takes a list of assembly json files and donwloads the processed contigs.")
parser.add_argument("sample", help="The list containing the url of assemblies as txt")
parser.add_argument("-n", "--name", type = str, help="Filename of tsv")
args = parser.parse_args()

print("Running")

in_sample = args.sample
filename = args.name

#Extracts links using sample ids
links = find_links(in_sample)

#Extract sample id from sample url
sample_list = []
for sample in lines_to_list(in_sample):
    sample_list.append(sample.split('/analyses/')[1].split('/')[0])

#Creates mapping tsv file with sample id and fasta id   
table = create_table_and_save_tsv(sample_list, links)

#print(table)

# Save the DataFrame to a TSV file without column names
table.to_csv(filename, sep='\t', index=False, header=False)
print(f"Table saved as {filename}")


print("DONE")
