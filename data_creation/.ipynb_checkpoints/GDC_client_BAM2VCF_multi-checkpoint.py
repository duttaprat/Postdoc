import os
import subprocess
from concurrent.futures import ProcessPoolExecutor

# Set the paths to your manifest and token files
MANIFEST_PATH = '/home/pdutta/Data/GDC/RAMANA/Manifest/bam_manifest.txt'
TOKEN_PATH = '/home/pdutta/Data/GDC/RAMANA/Token/gdc-user-token.2023-08-18T17_34_29.807Z.txt'
DOWNLOAD_DIR = '/home/pdutta/Data/GDC/RAMANA/data/BAM/'  # Adjust as needed
VCF_DIR = '/home/pdutta/Data/GDC/RAMANA/data/VCF/'

def download_and_process_bam(file_id):
    # Download the BAM file using gdc-client

    cmd_download = [
        'gdc-client',
        'download',
        '-m', MANIFEST_PATH,
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
        # Convert BAM to VCF
        vcf_filename = os.path.basename(bam_file_path).replace('.bam', '.vcf')
        vcf_file_path = os.path.join(VCF_DIR, vcf_filename)
        print(bam_file_path, vcf_file_path)

        # Convert BAM to VCF using samtools and bcftools
        cmd_bam_to_vcf = [
            'bcftools', 'mpileup', '-f', '/home/pdutta/Rekha_LabWork/Collaborator_work/RNA-seq/RNAseq_CommonFiles/STAR/Genome/GRCh38.primary_assembly.genome.fa',
            bam_file_path, '|', 'bcftools', 'call', '-mv', '>', vcf_file_path
        ]
        subprocess.call(' '.join(cmd_bam_to_vcf), shell=True)
        
        # Optional: Remove the downloaded folder to save space
        #shutil.rmtree(os.path.join(DOWNLOAD_DIR, file_id))
        os.remove(bam_file_path)
    else:
        print(f"No BAM file found for {file_id}")
    
    
    # cmd_bam_to_vcf = (
    # 'samtools mpileup -uf /home/pdutta/Data/Human_Genome_Data/GRCh38_latest_genomic.fna ' + bam_file_path +
    # ' | bcftools call -mv > ' + vcf_file_path
    # )
    # subprocess.call(cmd_bam_to_vcf, shell=True)

    

if __name__ == "__main__":
    # Read the manifest file to get a list of file IDs
    with open(MANIFEST_PATH, 'r') as f:
        lines = f.readlines()[1:]  # Skip the header
        file_ids = [line.split('\t')[0] for line in lines]
    
    # Use ProcessPoolExecutor to download and process BAM files in parallel
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        executor.map(download_and_process_bam, file_ids)
