import os, tarfile
import re
import multiprocessing
import pandas as pd
import pyBigWig
import pybedtools
import zipfile, gzip
from gtfparse import read_gtf
import pandas as pd
from Bio.Seq import Seq
import sys


reference_genome_path = "/home/pdutta/Data/Human_Genome_Data/UCSC/genome_by_sequence/"
data_path= "/home/pdutta/Data/Human_Genome_Data/UCSC/gencode_annotation/"
acceptor_path = "/home/pdutta/Data/Splicesites_data/Acceptor_splice_sites/"
non_acceptor_path = "/home/pdutta/Data/Splicesites_data/Non_acceptor_splice_sites/"
file_name="gencodeV38.bb"
gencode_path= "/home/pdutta/Data/Human_Genome_Data/UCSC/gencode_annotation/"
gtf_annotation_filename= "gencode.v38.annotation.gtf"
gff_annotation_filename= "gencode.v38.annotation.gff3"
file_path= data_path+file_name
gencode_gtf_file_path= gencode_path + gtf_annotation_filename
gencode_gff_file_path= gencode_path + gff_annotation_filename

print("FIRST CHECK")

print(file_path)

print (sys.argv[0], sys.argv[1]) 
chr_files=[]
chr_files.append(str(sys.argv[1]))
#chr_files = ["chrY.fa"]
#chr_files.append(str(sys.argv[1]))
print(chr_files)


bb= pyBigWig.open(file_path)
print("BigBed file uploaded...")


df = read_gtf(gencode_gtf_file_path)

print("Gencode annotation file uploaded....")


print(df.columns)
print(df.shape)


# In[7]:


print (df.head())




promoter_coordinates = [-249, 50]
core_promoter_coordinates= [-34, 35]
donor_coordinates= [-199, 200]
acceptor_coordinates= [-199, 200]





def remove_data(elements, distance):
    #elements = sorted(elements) # sorting the user input
    new_elements = [elements[0]] # make a new list for output
    for element in elements[1:]: # Iterate over the remaining elements...
        if element - new_elements[-1] >= distance: 
            # this is the condition you described above
            new_elements.append(element)
    return new_elements





def string_subs(row, temp_seq):
    start = row['coordinates'] -199
    end = row['coordinates'] + 200
    #print (start, end)
    sub_string = temp_seq[start:end]
    #print(sub_string)
    temp_seq = temp_seq.replace(sub_string, 70*'$')
    return temp_seq









# In[17]:


#global df_promoter, df_core_promoter, df_donor_splice, df_acceptor_splice

for file in chr_files:
    print (file)
    if file.endswith('.fa'):
        reference_genome_file_path = reference_genome_path + file
        chr_name = file.split('.')[0]
        print (chr_name)
        df_chr= df[df['seqname']==chr_name]
        df_new = df_chr[df_chr['feature']=='exon'].reset_index(drop=True)
        print (df_new.shape)
        #print (reference_genome_file_path)
        length= bb.chroms(chr_name)
        print (length)
        fasta = pybedtools.example_filename(reference_genome_file_path)
        temp_seq= pybedtools.BedTool.seq(chr_name+':1-'+str(length), fasta)
        print ("length:", len(temp_seq))
        
        pattern_acceptor= ['AG', 'ag']
        
        df_acceptor = pd.read_csv(acceptor_path+chr_name+"_acceptor_splice_set_sequence.csv")
        print(df_acceptor.shape)
        
        
       
        
        for index, row in df_acceptor.iterrows():
            if(index%500 ==0):
                print (index)
            #temp_seq = string_subs(row, temp_seq)
            start = row['coordinates'] -199
            end = row['coordinates'] + 200
            #print (start, end)
            sub_string = temp_seq[start:end]
            #print(sub_string)
            temp_seq = temp_seq.replace(sub_string, 400*'$')
        
        
        Non_acceptor_sequence = re.sub('[$]', '', temp_seq)
        print ("After removing all acceptor areas..", len(Non_acceptor_sequence))
        non_acceptor_indices = []
        for j in pattern_acceptor:
            non_acceptor_indices.append([match.start() for match in re.finditer(j, Non_acceptor_sequence)])
            
        print (len(non_acceptor_indices[0]))
        flat_list = [item for sublist in non_acceptor_indices for item in sublist]
        flat_list.sort()
        print ("@#@",len(flat_list))
        removed_flat_list = remove_data(flat_list, 400)
        print ("$$$", len(removed_flat_list))
        
        list_non_acceptor = []
        for i in removed_flat_list:
            temps =  Non_acceptor_sequence[i-40:i+39]
            list_non_acceptor.append(temps)
        pd_df = pd.DataFrame(list_non_acceptor)
        print (pd_df)
        pd_df.to_csv(non_acceptor_path+"/"+chr_name+"_non_acceptor.csv")
        


