import os
import subprocess
from concurrent.futures import ProcessPoolExecutor

# Set the paths to your manifest and token files
MANIFEST_PATH = '/home/pdutta/Data/GDC/RAMANA/Manifest/gdc_manifest_seven.txt'
TOKEN_PATH = '/home/pdutta/Data/GDC/RAMANA/Token/gdc-user-token.2023-08-18T17_34_29.807Z.txt'
DOWNLOAD_DIR = '/home/pdutta/Data/GDC/RAMANA/data/BAM/'  # Adjust as needed

def download_and_process_bam(file_id):
    # Download the BAM file using gdc-client
    bam_file_path = os.path.join(DOWNLOAD_DIR, f"{file_id}.bam")
    print(bam_file_path)
    cmd_download = [
        'gdc-client',
        'download',
        '-m', MANIFEST_PATH,
        '-t', TOKEN_PATH,
        '-d', DOWNLOAD_DIR,
        file_id
    ]
    subprocess.call(cmd_download)
    
#     # Convert BAM to VCF using samtools and bcftools
#     vcf_file_path = bam_file_path.replace('.bam', '.vcf')
#     print(vcf_file_path)
#     cmd_bam_to_vcf = [
#         'samtools', 'mpileup', '-uf', 'path_to_reference_genome.fasta',
#         bam_file_path, '|', 'bcftools', 'call', '-mv', '>', vcf_file_path
#     ]
#     subprocess.call(' '.join(cmd_bam_to_vcf), shell=True)
    
#     # Remove the BAM file
#     os.remove(bam_file_path)

if __name__ == "__main__":
    # Read the manifest file to get a list of file IDs
    with open(MANIFEST_PATH, 'r') as f:
        lines = f.readlines()[1:]  # Skip the header
        file_ids = [line.split('\t')[0] for line in lines]
    
    # Use ProcessPoolExecutor to download and process BAM files in parallel
    with ProcessPoolExecutor(max_workers=2) as executor:
        executor.map(download_and_process_bam, file_ids)
