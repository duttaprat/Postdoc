import pandas as pd
import gzip
from tqdm import tqdm
import os
import io
import csv


gz_file_path = '/data/projects/VCF_files/DBSNP/dbSnp155Details.tab.gz'
csv_folder_path = "/data/projects/VCF_files/DBSNP/ChrWise_Data"
output_path = "/data/private/pdutta/VCF_data/Merged_DBSNP_Chrwise/"


# DBSNP Columns

column_names = ['rs_ID', 'freqSourceCount', 'alleleCounts', 'alleleTotals', 'soTermCount', 'soTerms', 'clinVarCount', 'clinVarAccs', 'clinVarSigs', 'submitterCount', 'submitters', 'pubMedIdCount' ,'pubMedIds' ]




# In[8]:
def load_csv_data(file_path):
    """Load CSV data into a dictionary for quick lookup."""
    csv_data = {}
    headers = None
    with open(file_path, 'r') as f:
        reader = csv.reader(f)
        headers =next(reader)  # get header
        for row in reader:
            key = row[3]
            csv_data[key] = row
    return csv_data, headers


# In[9]:
def process_dbsnp_file(dbsnp_file_path, csv_data, csv_headers, outfile_path):
    """Process DBSNP file line by line and merge."""
    with gzip.open(dbsnp_file_path, 'rt') as dbsnp_file, open(outfile_path, 'w', newline='') as outfile:
        writer = csv.writer(outfile)
        dbsnp_headers = ['freqSourceCount', 'alleleCounts', 'alleleTotals', 'soTermCount', 'soTerms', 'clinVarCount', 'clinVarAccs', 'clinVarSigs', 'submitterCount', 'submitters', 'pubMedIdCount' ,'pubMedIds' ]
        writer.writerow(csv_headers + dbsnp_headers)
        for line in dbsnp_file:
            dbsnp_data = line.strip().split('\t')  # Assuming DBSNP data is tab-separated
            key = dbsnp_data[0]  # Assuming the column to merge on is the first one
            if key in csv_data:
                merged_row = csv_data[key] + dbsnp_data[1:]
                writer.writerow(merged_row)
                del csv_data[key]
                if not csv_data:
                    break


# In[10]:


files_generated = os.listdir(output_path)
print(files_generated, output_path)
output_files_names  = [ f.split('_')[0] for f in files_generated]
print(output_files_names)
print('#####')


files = os.listdir(csv_folder_path)
for file in files:
    if file.endswith('.csv') and 'M' not in file and file.split('_')[0] not in output_files_names:
        print("@@",file)
        file_path = os.path.join(csv_folder_path, file)
        csv_data, csv_headers = load_csv_data(file_path)
        # Create the output file path
        base_name = os.path.splitext(file)[0]
        outfile_path = output_path+base_name + '_merged.csv'
        # Process DBSNP file and merge
        process_dbsnp_file(gz_file_path, csv_data, csv_headers, outfile_path)


