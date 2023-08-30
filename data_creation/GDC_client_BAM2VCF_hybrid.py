import os
import subprocess
from concurrent.futures import ProcessPoolExecutor
import time
import argparse, sys, logging

logging.basicConfig(level=logging.INFO)

path_to_reference_genome = sys.argv[1]
manifest_file_name = sys.argv[2]

# Set the paths to your manifest and token files
MANIFEST_PATH = '/home/pdutta/Data/GDC/RAMANA/Manifest/'+ manifest_file_name
TOKEN_PATH = '/home/pdutta/Data/GDC/RAMANA/Token/gdc-user-token.2023-08-18T17_34_29.807Z.txt'
DOWNLOAD_DIR = '/home/pdutta/Data/GDC/RAMANA/data/BAM_multi/'  # Adjust as needed
VCF_DIR = '/home/pdutta/Data/GDC/RAMANA/data/VCF_multi/'



CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ['chrX', 'chrY']

def bam_to_vcf_chrom(bam_file_path, chrom):
    logging.info(f"Started processing for {chrom}")
    vcf_filename = os.path.basename(bam_file_path).replace('.bam', f'_{chrom}.vcf.gz')  # Appending chromosome to the filename
    print(vcf_filename , "is generating !!!")
    vcf_file_path = os.path.join(VCF_DIR, vcf_filename)

    # Convert BAM to VCF using bcftools with subprocess.Popen for the specific chromosome
    cmd_mpileup = ['bcftools', 'mpileup', '-r', chrom, '-f', path_to_reference_genome , bam_file_path]
    cmd_call = ['bcftools', 'call', '-mv', '-Oz', '-o', vcf_file_path]
    
    mpileup = subprocess.Popen(cmd_mpileup, stdout=subprocess.PIPE)
    call = subprocess.Popen(cmd_call, stdin=mpileup.stdout)
    mpileup.stdout.close()
    call.communicate()
    print(vcf_filename , "is generated !!!")
    logging.info(f"Finished processing for {chrom}")

def wrapper_function(chrom, bam_file_path):
    try:
        logging.info(f"Wrapper function called for {chrom}")
        bam_to_vcf_chrom(bam_file_path, chrom)
    except Exception as e:
        logging.error(f"An error occurred: {e}")
    

    
def process_chromosome(chrom, bam_file_path):
    wrapper_function(chrom, bam_file_path)
    
    
def concat_vcf_files(bam_file_path):
    base_name = os.path.basename(bam_file_path).replace('.bam', '')
    vcf_files = [os.path.join(VCF_DIR, f"{base_name}_{chrom}.vcf.gz") for chrom in CHROMOSOMES]
    
    merged_vcf = os.path.join(VCF_DIR, f"{base_name}.vcf.gz")
    
    cmd_concat = ['bcftools', 'concat', '-o', merged_vcf, '-Oz'] + vcf_files
    
    subprocess.call(cmd_concat)
    
    # Remove individual chromosome VCF files to save space
    for vcf_file in vcf_files:
        os.remove(vcf_file)

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
        # Process each chromosome in parallel
        try:
            with ProcessPoolExecutor(max_workers=24) as executor:
                logging.info("Starting processing.")
                executor.map(process_chromosome, CHROMOSOMES, [bam_file_path]*len(CHROMOSOMES))
                logging.info("Finished processing.")
        except Exception as e:
            logging.error(f"An error occurred: {e}")

        # with ProcessPoolExecutor(max_workers=24) as executor:
        #     executor.map(lambda chrom: bam_to_vcf_chrom(bam_file_path, chrom), CHROMOSOMES)
        
        # After processing all chromosomes, merge VCFs into one
        concat_vcf_files(bam_file_path)

        os.remove(bam_file_path)
    else:
        print(f"No BAM file found for {file_id}")
    
    
#     if bam_file_path:
#         # Convert BAM to VCF
#         vcf_filename = os.path.basename(bam_file_path).replace('.bam', '.vcf.gz')
#         vcf_file_path = os.path.join(VCF_DIR, vcf_filename)
#         print(bam_file_path, vcf_file_path)

#         # Convert BAM to VCF using samtools and bcftools
#         # Convert BAM to VCF using samtools and bcftools with subprocess.Popen
#         cmd_mpileup = ['bcftools', 'mpileup', '-f', '/data/projects/Resources/Gencode_genome_annotation/GRCh38.primary_assembly.genome.fa', bam_file_path]
#         cmd_call = ['bcftools', 'call', '-mv', '-Oz', '-o', vcf_file_path]
#         mpileup = subprocess.Popen(cmd_mpileup, stdout=subprocess.PIPE)
#         call = subprocess.Popen(cmd_call, stdin=mpileup.stdout)
#         mpileup.stdout.close()
#         call.communicate()
#         subprocess.call(' '.join(cmd_bam_to_vcf), shell=True)
        
#         # Optional: Remove the downloaded folder to save space
#         #shutil.rmtree(os.path.join(DOWNLOAD_DIR, file_id))
#         os.remove(bam_file_path)
#     else:
#         print(f"No BAM file found for {file_id}")
    
    
    # cmd_bam_to_vcf = (
    # 'samtools mpileup -uf /home/pdutta/Data/Human_Genome_Data/GRCh38_latest_genomic.fna ' + bam_file_path +
    # ' | bcftools call -mv > ' + vcf_file_path
    # )
    # subprocess.call(cmd_bam_to_vcf, shell=True)



    

if __name__ == "__main__":
    if not os.path.exists(DOWNLOAD_DIR):
        os.makedirs(DOWNLOAD_DIR)
    if not os.path.exists(VCF_DIR):
        os.makedirs(VCF_DIR)
        
    
    
    # Read the manifest file to get a list of file IDs
    with open(MANIFEST_PATH, 'r') as f:
        lines = f.readlines()[1:]  # Skip the header
        file_ids = [line.split('\t')[0] for line in lines]
    
    # Download the BAM files sequentially and convert in parallel
    # num_workers =int(os.cpu_count() /2) +10
    # print("Number of process parrallely working ", num_workers)
    for file_id in file_ids:
        download_and_process_bam(file_id)
