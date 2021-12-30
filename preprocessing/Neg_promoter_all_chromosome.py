import os, tarfile
import re
import multiprocessing
import pandas as pd
import pyBigWig
import pybedtools
import zipfile, gzip
from gtfparse import read_gtf
import pandas as pd
pd.set_option('display.max_columns', None)
from Bio.Seq import Seq
import sys




reference_genome_path = "/mnt/data05/shared/pdutta_data/Human_Genome_Data/UCSC/genome_by_sequence/"
data_path= "/mnt/data05/shared/pdutta_data/Human_Genome_Data/UCSC/gencode_annotation/"
promoter_path = "/mnt/data05/shared/pdutta_data/Human_Genome_Data/chromosome_wise_sequence_3000bp/Promoters/"
non_promoter_path = "/mnt/data05/shared/pdutta_data/Human_Genome_Data/chromosome_wise_sequence_3000bp/Non_promoters/"
file_name="gencodeV38.bb"
gencode_path= "/mnt/data05/shared/pdutta_data/Human_Genome_Data/gencode_annotation/"
gtf_annotation_filename= "gencode.v38.annotation.gtf"
gff_annotation_filename= "gencode.v38.annotation.gff3"
file_path= data_path+file_name
gencode_gtf_file_path= gencode_path + gtf_annotation_filename
gencode_gff_file_path= gencode_path + gff_annotation_filename
#reference_genome_file_path = reference_genome_path + reference_genome_name
data_save_path= "/mnt/data05/shared/pdutta_data/Human_Genome_Data/region_wise_sequence/"

chr_files= os.listdir(reference_genome_path)

bb= pyBigWig.open(file_path)

promoter_coordinates = [-300, 2700]



# In[5]:


def remove_data(elements, distance):
    #elements = sorted(elements) # sorting the user input
    new_elements = [elements[0]] # make a new list for output
    for element in elements[1:]: # Iterate over the remaining elements...
        if element - new_elements[-1] >= distance: 
            # this is the condition you described above
            new_elements.append(element)
    return new_elements



for file in chr_files:
    if file.endswith('.fa'):
        print ("File name:", file)
        reference_genome_file_path = reference_genome_path + file
        chr_name = file.split('.')[0]
        print (chr_name)


        length= bb.chroms(chr_name)
        print (length)
        fasta = pybedtools.example_filename(reference_genome_file_path)
        temp_seq= pybedtools.BedTool.seq(chr_name+':1-'+str(length), fasta)
        print ("length:", len(temp_seq))



        df_promoter = pd.read_csv(promoter_path+"/"+chr_name+"_promoter.csv")

        print (df_promoter.head())

        print(df_promoter.shape)



        for index, row in df_promoter.iterrows():
            if (index%5==0):
                print (index)
            if(row['strand']=='+'):
                start = row['TSS'] -300
                end = row['TSS'] + 2700
            else:
                start = row['TSS'] -2700
                end = row['TSS'] + 300

            #print (start, end)
            sub_string = temp_seq[start:end]
            temp_seq = temp_seq.replace(sub_string, 3001*'$')


        Non_promoter_sequence = re.sub('[$]', '', temp_seq)
        print ("@@##$$", len(Non_promoter_sequence))

    
        list_non_promoter = [Non_promoter_sequence[i:i+3001] for i in range(0, len(Non_promoter_sequence), 3001)]
        pd_df = pd.DataFrame(list_non_promoter)



        print (pd_df.head(45))


        print(pd_df.shape)



        pd_promoter = pd_df.drop_duplicates().reset_index(drop=True)
        print(pd_promoter.shape)
        pd_promoter.to_csv(non_promoter_path+"/"+chr_name+"_non_promoter.csv")







