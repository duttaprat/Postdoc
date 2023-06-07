import pandas as pd
import sys
import re
import os
import pyBigWig
#from gtfparse import read_gtf
#pd.set_option('display.max_columns', None)


chr_num=sys.argv[1]
data_df=pd.read_csv("/home/pdutta/Data/Splicesites_data/Acceptor_splice_sites/"+chr_num+"_acceptor_splice_set_sequence.csv", sep=",")
# print("{} Data is loaded".format(chr_num))



# # GET NON-ENHANCER COORDINATE FROM GENCODE FILE
gencodeFile= pyBigWig.open("/home/pdutta/Rekha_LabWork/Source_Folder/gencodeV38.bb")
# print("Gencode annotation is loaded")


chr_length= gencodeFile.chroms(chr_num)
# print("{} length is ".format(chr_num), chr_length)

data_df= data_df.head(10)

data_df = data_df.sort_values(by='coordinates')

# print(data_df.head())
start_pos=(data_df['coordinates']-200).to_list()
end_pos = ((data_df['coordinates']+200).to_list())
# print(start_pos)
# print(end_pos)

last=len(end_pos)
first_entry=chr_num+"\t"+"1"+"\t"+str(start_pos[0]-1)
print(first_entry)
for i in range(0,len(end_pos)-1):
    mid_entry=chr_num+"\t"+str(end_pos[i]+1)+"\t"+str(start_pos[i+1]-1)
    print(mid_entry)
last_entry=chr_num+"\t"+str(end_pos[last-1]+1)+"\t"+str(chr_length)
print(last_entry)