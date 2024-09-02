#!/usr/bin/python

# author : Gyan Prakash Mishra 
# Email : j12mishra@gmail.com

import os
import sys
from Bio import Entrez
import pandas as pd
import xml.etree.ElementTree as ET
import argparse
import GEOparse
import logging
import subprocess

# Suppress GEOparse messages by setting logging level to ERROR
logging.getLogger("GEOparse").setLevel(logging.ERROR)

# Enter your email and API key
Entrez.email = ()
Entrez.api_key = ()

def fetch_geo_samples(geo_accession,GSMid):
    """Fetch sample names and GSM IDs for a given GEO accession (e.g., GSEXXXXX)."""
    handle = Entrez.esearch(db="gds", term=geo_accession)
    record = Entrez.read(handle)
    handle.close()
    
    geo_id = record['IdList'][0]
    handle = Entrez.esummary(db="gds", id=geo_id)
    summary = Entrez.read(handle)
    handle.close()
    
    #print(summary[0])
    samples = []
    if "All" not in GSMid:
        samples = [item for item in summary[0]['Samples'] if item["Accession"] in  GSMid.split(",")]
    else:
        samples = summary[0]['Samples']
    

    samples_info = []

    # Extract Sample details from the GEO summary
    for sample in samples: #summary[0]['Samples']:
        sample_name = sample['Title']
        gsm_id = sample['Accession']
    
            
        samples_info.append((sample_name, gsm_id))

    return samples_info

def map_gsm_to_sra(gsm_id):
    """Map a GSM ID to its corresponding SRA ID."""
    handle = Entrez.esearch(db="sra", term=gsm_id)
    record = Entrez.read(handle)
    handle.close()
    
    sra_ids = record['IdList']
    #print(gsm_id)
    
    if not sra_ids:
        return None  # No SRA ID found for this GSM ID
    
    # Fetch the first SRA ID (you can modify this to get all)
    sra_id = sra_ids[0]
    
    # Fetch detailed information about this SRA ID
    handle = Entrez.efetch(db="sra", id=sra_id, rettype="xml")
    xml_data = handle.read()
    handle.close()
    
    # Parse the XML data to extract the SRA run accession
    root = ET.fromstring(xml_data)
    for run in root.iter("RUN"):
        return run.attrib['accession']  # Return the first SRA accession found
    
    return None

def map_gsm_to_org(gsm_id,GEOaccession):
    # Download and parse the SOFT file from GEO
    gse = GEOparse.get_GEO(geo=GEOaccession, destdir=".")

    # Extract the sample data
    # Assuming we want to extract and view the sample metadata
    # Extract sample metadata
    samples = gse.gsms
    
    samples = {k: v for k, v in samples.items() if k == gsm_id}

    ## Create a DataFrame to store sample information
    #sample_data = []
    for sample_name, sample in samples.items():
    # Extract relevant information from the sample
        sample_info = {
        #"Sample Name": sample_name,
        #"Title": sample.metadata.get("title", ""),
        #"Description": sample.metadata.get("description", ""),
        #"Source Name": sample.metadata.get("source_name_ch1", ""),
        "Organism": sample.metadata.get("organism_ch1"),
    }

    return(";".join(sample_info['Organism']))

def download_sra(sra_id, directory):
    """
    Downloads the SRA file if it's not already present in the specified directory.

    Parameters:
    sra_id (str): The ID of the SRA file to download.
    directory (str): The directory where the SRA file should be downloaded.
    """
    # Construct the file path where the SRA file would be saved
    sra_file_path = os.path.join(directory, f"{sra_id}.sra")

    # Check if the SRA file already exists
    if os.path.exists(sra_file_path):
        print(f"The file '{sra_file_path}' already exists.")
    else:
        # If the file does not exist, download it using fasterq-dump
        print(f"Downloading {sra_id} to {directory} .....")
        try:
            # Execute the fasterq-dump command to download the SRA file
            subprocess.run(["fasterq-dump", sra_id, "--outdir", directory], check=True)
            print(f"Successfully downloaded {sra_id} to {directory}.")
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while downloading {sra_id}: {e}")
        except Exception as e:
            print(f"An unexpected error occurred: {e}")

def main():
    # Create the parser
    parser = argparse.ArgumentParser(description='A script to download SRA fles and metadata from GEO accession number')

    # Add arguments
    parser.add_argument(
        '-GEO','--GEOaccession', 
        type=str, 
        required=True, 
        help='GEO Accession required'
    )
    
    parser.add_argument(
        '-GSM', '--GSMid', 
        type=str,
        required=True,
        help='default: "All"; GSM id (e.g) "GSM4297142,GSM4297143"'
    )

    parser.add_argument(
        '-d', '--outDir',
        type=str, 
        required=False,
        help='Enter directory to downlod SRA files (optional)'
    )
    parser.add_argument(
        '-e', '--EntrezEmail',
        type=str, 
        required=True,
        help='Enter Entrez email id'
    )
    parser.add_argument(
        '-key', '--EntrezKey',
        type=str, 
        required=True,
        help='Enter Entrez key'
    )


    # Parse the arguments
    args = parser.parse_args()
    
    Entrez.email = args.EntrezEmail
    Entrez.api_key = args.EntrezKey
    
    # Fetch sample names and GSM IDs from GEO
    geo_samples = fetch_geo_samples(args.GEOaccession,args.GSMid)
    #print(type(geo_samples))
    
    # Create a list to hold sample names and their corresponding SRA IDs
    data = []

    for sample_name, gsm_id in geo_samples:

        sra_id = map_gsm_to_sra(gsm_id)
        organism = map_gsm_to_org(gsm_id,args.GEOaccession)
        #data.append((sample_name,gsm_id,sra_id,organism))
        sra_file_path = ()
        reads = []
        fastq_file_1 = ()
        fastq_file_2 = ()
        single_end_fastq =()
        # Check if the directory exists
        if args.outDir:
            if not os.path.isdir(downloadDir):
                print(f"The directory '{downloadDir}' does not exist. Please provide a valid directory.")
            else:
            # Download the SRA file if the directory exists
                sra_file_path = os.path.join(downloadDir, f"{sra_id}.sra")
                download_sra(sra_id, downloadDir)

                # Construct expected fastq file paths
                fastq_file_1 = os.path.join(downloadDir, f"{sra_id}_1.fastq")
                fastq_file_2 = os.path.join(downloadDir, f"{sra_id}_2.fastq")
                single_end_fastq = os.path.join(downloadDir, f"{sra_id}.fastq")

            # Check if files already exist (for both paired and single-end data)
            if os.path.exists(fastq_file_1) and os.path.exists(fastq_file_2):
                print(f"Paired-end files '{fastq_file_1}' and '{fastq_file_2}' already exist.")
                subprocess.run(["gzip", fastq_file_1, fastq_file_2], check=True)
                data.append((sample_name, gsm_id,sra_id,organism,fastq_file_1 + ".gz", fastq_file_2 + ".gz"))
            elif os.path.exists(single_end_fastq):
                print(f"Single-end file '{single_end_fastq}' already exists.")
                subprocess.run(["gzip", single_end_fastq], check=True)
                data.append((sample_name, gsm_id,sra_id,organism,single_end_fastq + ".gz"))
        else:
            data.append((sample_name, gsm_id,sra_id,organism))

    # Convert the data into a DataFrame for easy viewing
    if args.outDir:
        if os.path.exists(single_end_fastq):
            df = pd.DataFrame(data, columns=["Sample Name", "GSM_ID", "SRA_ID", "Organism", "Read_1"])
            df.to_csv('data.csv', index=False)
        elif os.path.exists(fastq_file_1) and os.path.exists(fastq_file_2):
            df = pd.DataFrame(data, columns=["Sample Name", "GSM_ID", "SRA_ID", "Organism", "Read_1","Read_2"])
            df.to_csv('data.csv', index=False)
    else:
        df = pd.DataFrame(data, columns=["SampleName", "GSM_ID", "SRA_ID", "Organism"])
        df.to_csv('data.csv', index=False)

if __name__ == "__main__":
    main()
