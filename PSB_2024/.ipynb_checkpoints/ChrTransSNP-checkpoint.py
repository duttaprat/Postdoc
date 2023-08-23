import pandas as pd
import pysam, csv
from intervaltree import Interval, IntervalTree
import os


transcript_file_path = "/data/private/pdutta/PSB_Data/all_tc_input.csv"
output_path = "/data/private/pdutta/VCF_data/MergedTransVCF"
dbsnp_path = "/data/private/pdutta/VCF_data/Merged_DBSNP_Chrwise"


df_transcript = pd.read_csv(transcript_file_path, sep =",")
df_transcript['chromosome_name'] = 'chr' + df_transcript['chromosome_name'].astype(str)
print(df_transcript.shape)


print(df_transcript.groupby(['chromosome_name'])['chromosome_name'].count())



# Add a new column for mutations
df_transcript['mutations'] = ''

# Create a dictionary to store transcripts by chromosome
transcripts_by_chromosome = {chr_name: df for chr_name, df in df_transcript.groupby('chromosome_name')}

# # Create a dictionary to store IntervalTrees by chromosome
# trees_by_chromosome = {chr_name: IntervalTree() for chr_name in transcripts_by_chromosome}

# # In[ ]:


# Iterate over all the DBSNP files
for DBSNP_file in os.listdir(dbsnp_path):
    # Construct the full path to the DBSNP file
    DBSNP_path = os.path.join(dbsnp_path, DBSNP_file)
    print(DBSNP_path)

    # Extract the chromosome number from the DBSNP file name
    # Assumes the file name is in the format "DBSNP_chr#.csv"
    chr_name = DBSNP_file.split('_')[0]
    print("DBSNP_Chrosome: ", chr_name)

    # Check if there are any transcripts for this chromosome
    if chr_name not in transcripts_by_chromosome:
        continue

    # Get the IntervalTree for this chromosome
    tree = IntervalTree()

    # Read the DBSNP file line by line
    with open(DBSNP_path, 'r') as file:
        reader = csv.reader(file)
        next(reader)  # Skip the header row
        for row in reader:
            #print("%^",  int(row[1]), int(row[2]), row[3])
            #input()
            start, end, rs_id =  int(row[1]), int(row[2]), row[3]
            # Add a small amount to the end position if it's the same as the start position
            if start == end:
                end += 1
            # Add the interval to the tree
            tree[start:end] = rs_id

    

    # Get the transcripts for this chromosome
    df_chr = transcripts_by_chromosome[chr_name]
    print(df_chr)
    
     # Iterate over the transcripts
    for index, row in df_chr.iterrows():
        # Find the matching intervals in the DBSNP data
        
        matches = tree[row['transcript_start']:row['transcript_end']]

        # Add the mutation information to the transcript
        df_transcript.at[index, 'mutations'] = ', '.join([str(match.data) for match in matches])



    # Save the merged dataframe for this chromosome to a new csv file
    df_transcript.loc[df_transcript['chromosome_name'] == chr_name].to_csv(f'{output_path}/merged_{chr_name}.csv', index=False)









            
#     # Fill the tree with intervals from the DBSNP file
#     for index, row in df_DBSNP.iterrows():
#         # Check if the start and end positions are the same
#         if row['start'] == row['end']:
#             # Add a small amount to the end position
#             tree[row['start']:row['end'] + 1] = row['name']
#         else:
#             tree[row['start']:row['end']] = row['name']


#     chr_transcripts = transcripts_by_chromosome.get_group(chr_name)
#     # Iterate over the transcripts
#     for index, row in chr_transcripts.iterrows():
#         # Find the matching intervals in the DBSNP data
#         matches = tree[row['start']:row['end']]

#         # Add the mutation information to the transcript
#         df_transcripts.at[index, 'mutations'] = ', '.join([str(match.data) for match in matches])
#     df_transcripts.loc[df_transcripts['chr_name'] == chr_name].to_csv(f'{output_path}/merged_{chr_name}.csv', index=False)


# # In[ ]:




