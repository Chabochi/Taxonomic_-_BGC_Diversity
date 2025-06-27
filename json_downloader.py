"""
Script to download a FASTA (.fasta.gz) file from a JSON-derived URL.

Given a single tab-separated line containing a sample ID and its corresponding download URL, 
this script downloads the FASTA file and saves it to a specified output directory. The file is 
named using the sample ID and saved in `.fasta.gz` format.

Usage:
    python script.py "<SAMPLE_ID>\t<URL>" -o output_directory

Arguments:
    line (str): A tab-separated string with a sample ID and download URL.
    -o, --output (str): Path to the output directory where the file will be saved.

Example:
    python script.py "SAMPLE123\thttps://example.com/contigs.fasta.gz" -o ./downloads
"""


#Function for downloading fasta file from json url
def download_fasta(line, output):
    
    #Output folder to store the download files
    output_folder = output

    #Create the output folder if it doesn't exist
    os.makedirs(output_folder, exist_ok = True)
        
    print("The line is", line)
    # split data by tab 
    if len(line.split('\t')) > 1:
        link = line.split('\t')[1]
        print("The link is", link)
        index = 0
    
        try:
            #Parse tge URL to separate components
            parsed_url = urlparse(link)

            #Extract the filename from the URL
            filename = line.split('\t')[0] + ".fasta.gz"

            #Full path to save the file
            output_path = os.path.join(output_folder, filename)

            #Send HTTP GET request to download the file
            response = requests.get(link)

            #Check if request was succesful
            if response.status_code == 200:

                #Open the output file in binary mode and write the content
                with open(output_path, 'wb') as file:
                    
                    #file.write(response.content)
                    file.write(response.content)

                index = index + 1
                print(f"Downloaded {index}: {filename} to {output_path}")

            else:
                print(f"Failed to download {link}. Status code: {response.status_code}")
        
        except Exception as e:
            print(f"Error ocurred while downloading {link}: {e}")
    
    else:
        print("Link does not exist")

    return

#MAIN
import argparse
import os
from urllib.parse import urlparse
import requests

parser = argparse.ArgumentParser(description = "A script that takes a tab separated line with id and url and donwloads the processed contigs.")
parser.add_argument("line", help="The line containing the url of assemblies as txt")
parser.add_argument("-o", "--output", type = str, help="Output directory")
args = parser.parse_args()

print("Running")

in_line = args.line
in_out = args.output       

#Downloads fasta files to output directory
download_fasta(in_line, in_out)

print("DONE")
