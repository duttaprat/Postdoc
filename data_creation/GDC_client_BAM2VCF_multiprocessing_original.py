import os
import subprocess
from concurrent.futures import ProcessPoolExecutor

# Set the paths to your manifest and token files
MANIFEST_PATH = '/data/projects/BAM_files/manifest_files/bam_manifest.txt'
TOKEN_PATH = '/data/projects/BAM_files/token/gdc-user-token.2023-08-18T17_34_29.807Z.txt'
DOWNLOAD_DIR = '/data/projects/BAM_files/data/BAM/'  # Adjust as needed
VCF_DIR = '/data/projects/BAM_files/data/VCF/'

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
        # Convert BAM to VCF
        vcf_filename = os.path.basename(bam_file_path).replace('.bam', '.vcf.gz')
        vcf_file_path = os.path.join(VCF_DIR, vcf_filename)
        print(bam_file_path, vcf_file_path)

        # Convert BAM to VCF using samtools and bcftools
        # Convert BAM to VCF using samtools and bcftools with subprocess.Popen
        cmd_mpileup = ['bcftools', 'mpileup', '-f', '/data/projects/Resources/Gencode_genome_annotation/GRCh38.primary_assembly.genome.fa', bam_file_path]
        cmd_call = ['bcftools', 'call', '-mv', '-Oz', '-o', vcf_file_path]
        mpileup = subprocess.Popen(cmd_mpileup, stdout=subprocess.PIPE)
        call = subprocess.Popen(cmd_call, stdin=mpileup.stdout)
        mpileup.stdout.close()
        call.communicate()
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
    num_workers =os.cpu_count() /2 +10
    print("Number of process parrallely working ", num_workers)
    # Read the manifest file to get a list of file IDs
    with open(MANIFEST_PATH, 'r') as f:
        lines = f.readlines()[1:]  # Skip the header
        file_ids = [line.split('\t')[0] for line in lines]
    
    # Use ProcessPoolExecutor to download and process BAM files in parallel
    with ProcessPoolExecutor() as executor:
        executor.map(download_and_process_bam, file_ids)
