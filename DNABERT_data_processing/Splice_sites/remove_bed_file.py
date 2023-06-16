from pybedtools import BedTool
import os, glob
import pandas as pd


os.chdir("/home/pdutta/Data/Splicesites_data/Non_donor_splice_sites/Coordinates")
# Load your bed file
file_extension = '.txt'
all_filenames = [i for i in glob.glob(f"*{file_extension}")]



# Define a function to check if start is greater than end
def filter_fn(interval):
    return interval.start <= interval.end

for file in all_filenames:
    df= pd.read_csv(file, sep="\t", names= ["Chr_name","Start","End"])
    #print(df.head(5))
    df_temp= df[df['Start'] < df['End']]
    #print (df_temp.head(5))
    #print(df.shape, df_temp.shape)
    df_temp.to_csv(file.split('.')[0]+'.bed', sep='\t', header=False, index=False)
