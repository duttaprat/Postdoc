import os
import pandas as pd

# List all chromosomes
chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y']


# Loop over all chromosome folders
for chr_name in chromosomes:
    dir_name = f"/data/projects/PSB/DNABERT_data/CoreProm/chr{chr_name}"
    outdir_name = f"/data/projects/PSB/DNABERT_data/CoreProm/chr{chr_name}"
    file_name = os.path.join(dir_name, "dev.tsv")
    df = pd.read_csv(file_name, sep='\t')
    print(df.shape)
    df['Sequence'] = df['Sequence'].str.upper()
    print(chr_name)
    df.to_csv(file_name, sep='\t',  index=False)

#     # Check if file exists
#     if os.path.isfile(file_name):
#         # Read the file in chunks
#         chunk_size = 500000  # You may need to adjust this value
#         chunks = pd.read_csv(file_name, sep='\t', chunksize=chunk_size)

#         folder_index = 0
#         file_index = 0
#         for chunk in chunks:
#             if file_index % (chunk_size*10) == 0:  # Create a new subfolder for every 10 chunks
#                 folder_index += 1

#             # Create subfolder if it doesn't exist
#             subfolder = os.path.join(outdir_name, f"subfolder{folder_index}")
#             os.makedirs(subfolder, exist_ok=True)

#             # Write chunk to new dev.tsv in subfolder
#             chunk.to_csv(os.path.join(subfolder, "dev.tsv"), sep='\t', mode='a', index=False)

#             file_index += chunk_size
#     else:
#         print(f"File not found: {file_name}")
