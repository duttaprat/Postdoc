import pandas as pd
import numpy as np
import os, io
import pysam
import pickle, subprocess
import pybedtools
print(pybedtools.__file__)
pybedtools.helpers.set_tempdir('/data/private/pdutta_new/pybedtools_temp')
from IPython.display import display, HTML
import multiprocessing




top_models_df= pd.read_csv("/home/campus.stonybrook.edu/pdutta/Github/Postdoc/DNABERT_data_processing/TFBS/important_tfbs.txt", sep="\t", header=None)
for index, row in top_models_df.iterrows():
    print(index, row[0])
#print(top_models_df)
input()
cancer_type="Brain" 

vcf_folder_path="/data/projects/GDC_Cancer_Wise/New_data/{}/Downloaded_files/VCF".format(cancer_type)
reference_genome_path="/data/projects/Resources/Gencode_genome_annotation/GRCh38.primary_assembly.genome.fa"
chip_seq_base_path= "/data/projects/DNABERT_snv/Manuscript_11_2023/TFBS_fine_tune_data/Data_Jan_2024/cleaned_names_0524/1_1_pos_neg"
output_base_path = "/data/projects/GDC_Cancer_Wise/New_data/{}/Generated_files/Intersected_Data/300bp_TFBS".format(cancer_type)



reference_fasta = pysam.FastaFile(reference_genome_path)



def get_vcf_gz_files_except_logs(root_folder):
    all_files = []

    # Walk through the directory tree
    for dirpath, dirnames, filenames in os.walk(root_folder):
        # If "logs" is in dirnames, remove it to avoid traversing it
        if 'logs' in dirnames:
            dirnames.remove('logs')

        # Add only the filenames with the extension .vcf.gz in the current directory to the all_files list
        for filename in filenames:
            if filename.endswith('.vcf.gz'):
                all_files.append(os.path.join(dirpath, filename))
    return all_files





def reindex_vcf(vcf_path):
    try:
        subprocess.run(['tabix', '-p', 'vcf', vcf_path], check=True)
        print(f"Re-indexed successfully: {vcf_path}")
    except subprocess.CalledProcessError as e:
        print(f"Failed to re-index {vcf_path}: {e}")





def vcf_to_dataframe(vcf_path):
    """
    Convert a .vcf.gz file into a pandas DataFrame.

    Parameters:
    - vcf_path (str): path to the .vcf.gz file

    Returns:
    - pd.DataFrame: VCF data as a DataFrame
    """

    try:
        # Open the VCF file
        vcf_file = pysam.VariantFile(vcf_path)
    except ValueError as e:
        print(f"Failed to open file {vcf_path}: {e}")
        with open("/data/projects/GDC_Cancer_Wise/Brain/Data/"+vcf_path.split("/")[-2]+"_error_log.txt", "a") as error_log:
            error_log.write(f"{vcf_path}\n")
        return None

    # Extracting the data and the columns
    data = []
    
    columns =  ["CHROM", "START_POS", "ID", "REF", "ALT", "QUAL", "FILTER"] + list(vcf_file.header.info.keys()) 
    for record in vcf_file:
        basic_data = [record.chrom, record.pos, record.id, record.ref, 
                      ','.join(str(alt) for alt in record.alts), record.qual, record.filter.keys()[0] if record.filter.keys() else 'PASS']
        row_data = [record.info.get(key) for key in vcf_file.header.info.keys()]
        data.append(basic_data + row_data)

    df = pd.DataFrame(data, columns=columns)
    df["START_POS"] = df["START_POS"]-1
    end = df["START_POS"] + df['REF'].str.len()
    df.insert(2, 'END_POS', end)

    # Close the VCF file
    vcf_file.close()

    return df


# In[19]:




# In[ ]:


def process_tfbs_data(sub_df):
    for index, row in sub_df.iterrows():
        file_path = f"{chip_seq_base_path}/{row[0]}/300bp_unique_balanced.bed"
        print(file_path)
        tfbs_output_path = f"{output_base_path}/{row[0]}"
        if not os.path.exists(tfbs_output_path):
            os.makedirs(tfbs_output_path)
            print(tfbs_output_path, " folder is generated..!!! \n")
        df_tfbs = pd.read_csv(file_path, sep="\t")
        #print(df_tfbs.head())
        #input()
        TFBS_bed = pybedtools.BedTool.from_dataframe(df_tfbs)
        intersected_vcf_data = {}
        df_statistics = pd.DataFrame(columns=["filename", "Patient_ID", 'work_flow', 'VCF_instance', 'VCF_feature', 'Intersected_instances', 'VCF_column_names'])

        for file_path in files:
            file_name = file_path.split('/')[-1]
            parts = file_name.split('.')
            #print(parts)
            #reindex_vcf(file_path)
            df_vcf = vcf_to_dataframe(file_path)
            if df_vcf is None:
                continue
            vcf_bed = pybedtools.BedTool.from_dataframe(df_vcf)
            intersect_vcf_acceptor = TFBS_bed.intersect(vcf_bed, wa=True, wb=True)
            column_list = df_tfbs.columns.to_list() + df_vcf.columns.to_list()
            df_intersection = intersect_vcf_acceptor.to_dataframe(names=column_list)
            df_intersection = df_intersection[(df_intersection['REF'].str.len() < 10) & (df_intersection['ALT'].str.len() < 10)]

            patient_ID = parts[0]
            work_flow = parts[2]
            instance = df_vcf.shape[0]
            feature = df_vcf.shape[1]
            intersected_instance = df_intersection.shape[0]
            feature_list = list(df_vcf.columns)
            print(patient_ID, work_flow, instance, feature, intersected_instance, feature_list)
            df_statistics.loc[len(df_statistics)] = [file_name, patient_ID, work_flow, instance, feature, intersected_instance, feature_list]
            intersected_vcf_data[patient_ID+"_"+work_flow] = df_intersection

        df_statistics.to_csv(tfbs_output_path+"/VCF_statistics.tsv", sep="\t", index=False)
        with open(tfbs_output_path+"/intersected_vcf_data.pkl", "wb") as file:
            pickle.dump(intersected_vcf_data, file)



def split_dataframe(df, num_splits):
    return np.array_split(df, num_splits)
            
files = get_vcf_gz_files_except_logs(vcf_folder_path)
print(len(files))




# Number of processes you want to run in parallel
num_processes = 34
sub_dfs = split_dataframe(top_models_df, num_processes)


# Creating a pool of processes
pool = multiprocessing.Pool(processes=num_processes)

# Distribute the data processing function to the pool
pool.map(process_tfbs_data, sub_dfs)

# Close the pool and wait for the work to finish
pool.close()
pool.join()