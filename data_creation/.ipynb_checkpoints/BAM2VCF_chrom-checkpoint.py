import os
import subprocess
from concurrent.futures import ProcessPoolExecutor
import logging

logging.basicConfig(level=logging.INFO)


parent_dir = "/home/pdutta/Data/GDC/RAMANA/data/BAM_multi/"
VCF_DIR = '/home/pdutta/Data/GDC/RAMANA/data/VCF_multi/'

CHROMOSOMES = [f"chr{i}" for i in range(1, 23)] + ['chrX', 'chrY']  # Your list of chromosomes

def bam_to_vcf_chrom(bam_file_path, chrom):
    vcf_filename = os.path.basename(bam_file_path).replace('.bam', f'_{chrom}.vcf.gz')
    vcf_file_path = os.path.join(VCF_DIR, vcf_filename)

    cmd_mpileup = ['bcftools', 'mpileup', '-r', chrom, '-f', '/home/pdutta/Rekha_LabWork/Collaborator_work/RNA-seq/RNAseq_CommonFiles/STAR/Genome/GRCh38.primary_assembly.genome.fa', bam_file_path]
    cmd_call = ['bcftools', 'call', '-mv', '-Oz', '-o', vcf_file_path]
    
    mpileup = subprocess.Popen(cmd_mpileup, stdout=subprocess.PIPE)
    call = subprocess.Popen(cmd_call, stdin=mpileup.stdout)
    mpileup.stdout.close()
    call.communicate()

def concat_vcf_files(bam_file_path):
    base_name = os.path.basename(bam_file_path).replace('.bam', '')
    vcf_files = [os.path.join(VCF_DIR, f"{base_name}_{chrom}.vcf.gz") for chrom in CHROMOSOMES]
    
    merged_vcf = os.path.join(VCF_DIR, f"{base_name}.vcf.gz")
    
    cmd_concat = ['bcftools', 'concat', '-o', merged_vcf, '-Oz'] + vcf_files
    
    subprocess.call(cmd_concat)
    
    # Remove individual chromosome VCF files to save space
    for vcf_file in vcf_files:
        os.remove(vcf_file)

def wrapper_function(chrom):
    logging.info(f"Processing {chrom}")
    bam_to_vcf_chrom(bam_file_path, chrom)

if __name__ == "__main__":
    # if not os.path.exists(VCF_DIR):
    #     os.makedirs(VCF_DIR)
    
    for root, dirs, files in os.walk(parent_dir):
        for file in files:
            if file.endswith(".bam"):
                bam_file_path = os.path.join(root, file)
                print(bam_file_path)
#                 try:
#                     with ProcessPoolExecutor(max_workers=24) as executor:
#                         logging.info("Starting processing.")
#                         executor.map(wrapper_function, CHROMOSOMES)
#                         logging.info("Finished processing.")
#                 except Exception as e:
#                     logging.error(f"An error occurred: {e}")


#                 concat_vcf_files(bam_file_path)
#                 os.remove(bam_file_path)

    
    




