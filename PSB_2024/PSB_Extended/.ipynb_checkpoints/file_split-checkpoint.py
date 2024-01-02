import os
import pandas as pd

# Load the TSV file
file_path = '/data/projects/DNABERT_snv/Manuscript_11_2023/DBSNP_output/DNABERT_data/core_prom/dev.tsv'  # Update this with the path to your dev.tsv file
data = pd.read_csv(file_path, sep='\t', header=None)  # Assuming the file has no header
print("Data loaded..")
print(data.head())
print(data.shape)
# Calculate the size of each split
split_size = len(data) // 20

# Create 20 folders and save parts of the file
for i in range(20):
    # Define folder name
    #folder_path_name = f'/data/projects/DNABERT_snv/Manuscript_11_2023/DBSNP_output/DNABERT_data/core_prom/folder_{i+1}'
    #os.makedirs(folder_path_name, exist_ok=True)

    # Calculate start and end indices for the split
    start = i * split_size
    end = (i + 1) * split_size if i < 19 else len(data)  # Handle the last split

    # Extract the part of the data
    part = data.iloc[start:end]
    print(part.head())
    print(part.shape)

    # Save this part to a new TSV in the corresponding folder
    #part.to_csv(f'{folder_path_name}/dev.tsv', sep='\t', index=False, header=False)

print("Splitting completed.")