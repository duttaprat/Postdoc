import os, tarfile
import re
import multiprocessing
import pandas as pd
import pyBigWig
import pybedtools
import zipfile, gzip
from gtfparse import read_gtf
import pandas as pd
#pd.set_option('display.max_columns', None)
from Bio.Seq import Seq
import sys


# In[2]:


reference_genome_path = "/mnt/data05/shared/pdutta_data/Human_Genome_Data/UCSC/genome_by_sequence/"
data_path= "/mnt/data05/shared/pdutta_data/Human_Genome_Data/UCSC/gencode_annotation/"
donor_path = "/mnt/data05/shared/pdutta_data/Human_Genome_Data/Donor_splice_sites/"
non_donor_path = "/mnt/data05/shared/pdutta_data/Human_Genome_Data/Non_donor_splice_sites/"
file_name="gencodeV38.bb"
gencode_path= "/mnt/data05/shared/pdutta_data/Human_Genome_Data/gencode_annotation/"
gtf_annotation_filename= "gencode.v38.annotation.gtf"
gff_annotation_filename= "gencode.v38.annotation.gff3"
file_path= data_path+file_name
gencode_gtf_file_path= gencode_path + gtf_annotation_filename
gencode_gff_file_path= gencode_path + gff_annotation_filename
#reference_genome_file_path = reference_genome_path + reference_genome_name
data_save_path= "/mnt/data05/shared/pdutta_data/Human_Genome_Data/region_wise_sequence/"


bb= pyBigWig.open(file_path)
print("BigBed file uploaded...")


df = read_gtf(gencode_gtf_file_path)

print("Gencode annotation file uploaded....")


print(df.columns)
print(df.shape)



# In[6]:


df.head()


# In[8]:


promoter_coordinates = [-249, 50]
core_promoter_coordinates= [-34, 35]
donor_coordinates= [-199, 200]
acceptor_coordinates= [-199, 200]


# In[13]:


def remove_data(elements, distance):
    #elements = sorted(elements) # sorting the user input
    new_elements = [elements[0]] # make a new list for output
    for element in elements[1:]: # Iterate over the remaining elements...
        if element - new_elements[-1] >= distance: 
            # this is the condition you described above
            new_elements.append(element)
    return new_elements


# In[14]:


def string_subs(row, temp_seq):
    start = row['coordinates'] -35
    end = row['coordinates'] + 34
    #print (start, end)
    sub_string = temp_seq[start:end]
    #print(sub_string)
    temp_seq = temp_seq.replace(sub_string, 70*'$')
    return temp_seq


# In[9]:





# In[18]:


#chr_file= "chr2.fa"
chr_file = sys.argv[1]

print (chr_file)
reference_genome_file_path = reference_genome_path + chr_file
chr_name = chr_file.split('.')[0]
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

pattern_donor= ['GT', 'gt']

df_donor = pd.read_csv(donor_path+chr_name+"_donor_splice_set.csv")


print(df_donor.shape)



for index, row in df_donor.iterrows():
    if (index%5==0):
        print (index)
    start = row['coordinates'] -35
    end = row['coordinates'] + 34
    #print (start, end)
    sub_string = temp_seq[start:end]
    #print(sub_string)
    temp_seq = temp_seq.replace(sub_string, 70*'$')


Non_donor_sequence = re.sub('[$]', '', temp_seq)
print ("@@##$$", len(Non_donor_sequence))
non_donor_indices = []
for j in pattern_donor:
    non_donor_indices.append([match.start() for match in re.finditer(j, Non_donor_sequence)])

print (len(non_donor_indices[0]))
flat_list = [item for sublist in non_donor_indices for item in sublist]
print ("test")
flat_list.sort()
print ("test1")
print ("@#@",len(flat_list))
removed_flat_list = remove_data(flat_list, 70)
print ("$$$", len(removed_flat_list))
print (removed_flat_list)

list_non_donor = []
for i in removed_flat_list:
    temps =  Non_donor_sequence[i-35:i+34]
    list_non_donor.append(temps)
pd_df = pd.DataFrame(list_non_donor)
print (pd_df)
pd_df.to_csv(non_donor_path+chr_name+"_non_donor.csv")


