import sys, os
import pandas as pd
import numpy as np
from tqdm import tqdm

chr_file = sys.argv[1]
chr_name = chr_file.split('_')[0]
input_path = "/data/projects/Enhancer/Data_processing/Chromosome_Wise_Files/Merged_Files/"
input_file_path  = input_path + chr_file
output_path =  "/data/projects/Enhancer/Data_processing/Chromosome_Wise_Files/Smaller_chunks_removed_merged_files"


print(chr_name)
df   = pd.read_csv(input_file_path, sep= ",")

#df = pd.read_csv('/data/projects/Enhancer/Data_processing/TE-DBs.bed', sep= '\t', names= ['chr', 'start', 'end'])
print(df.shape)
print(df.head())
df = df.sort_values(by=['start'], ascending=True).reset_index(drop=True)
df['length'] = df['end']-df['start']
print(df.head())
print(df.shape)
df= df[df['length']<5000].reset_index(drop=True)
print(df.head())
print("Number of rows after removing greater than 5000 sequence length", df.shape)

end_list = df['end'].to_list()
count_list = []   

for i in tqdm(range(len(end_list))):
    count = 0
    for j in reversed(end_list[:i]):
        if end_list[i] < j:  
            count += 1
            break
    count_list.append(count)


df['count_list']=count_list
df_remove = df[df['count_list']==0]
print(df_remove.shape)

df_remove.to_csv(output_path+"/"+chr_name+"_drop_small_chunks.csv", sep= ",", index=False)






