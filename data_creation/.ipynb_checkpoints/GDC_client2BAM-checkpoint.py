import os
import subprocess
from concurrent.futures import ProcessPoolExecutor
import time
import argparse, sys, logging

logging.basicConfig(level=logging.INFO)

# Set the paths to your manifest and token files
TOKEN_PATH =  sys.argv[1]
MANIFEST_PATH = sys.argv[2]
DOWNLOAD_DIR =  '/data/projects/GDC_Cancer_Wise/'+sys.argv[3]+'/Data/BAM_NEED_TO_BE_PROCESSED/'  # Adjust as needed





def download_and_process_bam(file_id):
    # Download the BAM file using gdc-client
    print(file_id)
    cmd_download = [
        'gdc-client',
        'download',
        '-t', TOKEN_PATH,
        '-d', DOWNLOAD_DIR,
        file_id
    ]
    subprocess.call(cmd_download)
    
    bam_file_path = None
    for root, dirs, files in os.walk(os.path.join(DOWNLOAD_DIR, file_id)):
        for file in files:
            if file.endswith('.bam'):
                bam_file_path = os.path.join(root, file)
                break
    if bam_file_path:
        print(bam_file_path, "downloaded..!!!")
    else:
        print(f"No BAM file found for {file_id}")
    


    

if __name__ == "__main__":
    if not os.path.exists(DOWNLOAD_DIR):
        os.makedirs(DOWNLOAD_DIR)

        
    
    
    # Read the manifest file to get a list of file IDs
    with open(MANIFEST_PATH, 'r') as f:
        lines = f.readlines()[1:]  # Skip the header
        file_ids = [line.split('\t')[0] for line in lines]
    
    # Download the BAM files sequentially and convert in parallel
    # num_workers =int(os.cpu_count() /2) +10
    # print("Number of process parrallely working ", num_workers)
    for file_id in file_ids:
        download_and_process_bam(file_id)
