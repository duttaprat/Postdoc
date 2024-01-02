import os, io
from collections import defaultdict
import pandas as pd
import pysam
import pickle
import argparse, sys, logging
import pybedtools
pybedtools.helpers.set_tempdir('/data/projects/temp')

# In[2]:
non_coding_region = "Donor"
folder_name = sys.argv[1]
parent_folder_path  = "/data/projects/GDC_Cancer_Wise/Brain/Data/all_vcf_subfolders/"+ folder_name
reference_file_path = "/data/projects/Resources/Gencode_genome_annotation/GRCh38.primary_assembly.genome.fa"


reference_fasta = pysam.FastaFile(reference_file_path)
df_donor_bed = pd.read_csv("/data/private/pdutta/DNABERT_data/Splice_sites/{}_raw_files/Unique_{}_RAW_Combined_BED.tsv".format(non_coding_region, non_coding_region), sep="\t")
donor_bed = pybedtools.BedTool.from_dataframe(df_donor_bed)

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


# In[5]:


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


# In[6]:


files = get_vcf_gz_files_except_logs(parent_folder_path)
intersected_vcf_data = {}
df_statistics  = pd.DataFrame(columns= ['filename', 'gdc_tag', 'VCF_instance','VCF_feature','Intersected_instances', 'VCF_column_names' ])
for file_path in files:
    file_name =  file_path.split('/')[-1]
    parts = file_name.split('_')
    print(parts)
    df_vcf = vcf_to_dataframe(file_path)
    if df_vcf is None:
        continue
    vcf_bed = pybedtools.BedTool.from_dataframe(df_vcf)
    intersect_vcf_donor = donor_bed.intersect(vcf_bed, wa=True, wb=True)
    column_list  = df_donor_bed.columns.to_list() + df_vcf.columns.to_list()
    df_intersection = intersect_vcf_donor.to_dataframe(names=column_list)
    #print(df_intersection)
    print(df_statistics)
    #display(df_vcf)
    gdc_tag = parts[-4]
    instance =  df_vcf.shape[0]
    feature = df_vcf.shape[1]
    intersected_instance = df_intersection.shape[0]
    feature_list = list(df_vcf.columns)
    #print(file_name, gdc_tag, instance, feature,intersected_instance,  feature_list)
    df_statistics.loc[len(df_statistics)] = [file_name, gdc_tag, instance, feature,intersected_instance,  feature_list]
    intersected_vcf_data[file_name] = df_intersection





output_path = "/data/projects/GDC_Cancer_Wise/Brain/Data/Preprocessed_file/{}/".format(non_coding_region) + folder_name
if not os.path.exists(output_path):
    os.makedirs(output_path)


# In[12]:


df_statistics.to_csv(output_path+"/VCF_brain_statistics.tsv", sep="\t", index=False)


# In[13]:


with open(output_path+"/vcf_data.pkl", "wb") as file:
    pickle.dump(intersected_vcf_data, file)