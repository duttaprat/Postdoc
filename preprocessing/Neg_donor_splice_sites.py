#!/usr/bin/env python
# coding: utf-8

# In[1]:


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


# In[2]:


reference_genome_path = "/mnt/data05/shared/pdutta_data/Human_Genome_Data/UCSC/genome_by_sequence/"
data_path= "/mnt/data05/shared/pdutta_data/Human_Genome_Data/UCSC/gencode_annotation/"
donor_path = "/mnt/data05/shared/pdutta_data/Human_Genome_Data/Donor_splice_sites/"
non_donor_path = "/mnt/data05/shared/pdutta_data/Human_Genome_Data/Non_donnor_splice_sites/"
file_name="gencodeV38.bb"
gencode_path= "/mnt/data05/shared/pdutta_data/Human_Genome_Data/gencode_annotation/"
gtf_annotation_filename= "gencode.v38.annotation.gtf"
gff_annotation_filename= "gencode.v38.annotation.gff3"
file_path= data_path+file_name
gencode_gtf_file_path= gencode_path + gtf_annotation_filename
gencode_gff_file_path= gencode_path + gff_annotation_filename
#reference_genome_file_path = reference_genome_path + reference_genome_name
data_save_path= "/mnt/data05/shared/pdutta_data/Human_Genome_Data/region_wise_sequence/"


# In[3]:


bb= pyBigWig.open(file_path)


# In[4]:


df = read_gtf(gencode_gtf_file_path)


# In[5]:


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


chr_file= "chr2.fa"
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



for index, row in df_donor.head(15).iterrows():
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
for j in pattern_acceptor:
    non_donor_indices.append([match.start() for match in re.finditer(j, Non_donor_sequence)])

print (len(non_donor_indices[0]))
flat_list = [item for sublist in non_donor_indices for item in sublist]
print ("test")
flat_list.sort()
print ("test1")
print ("@#@",len(flat_list))
removed_flat_list = remove_data(flat_list, 70)
print ("$$$", len(removed_flat_list))
print (removed_flat_list[0:15])

list_non_donor = []
for i in removed_flat_list[0:10]:
    temps =  Non_donor_sequence[i-35:i+34]
    list_non_donor.append(temps)
pd_df = pd.DataFrame(list_non_donor)
print (pd_df)
#pd_df.to_csv(non_acceptor_path+"/"+chr_name+"_non_acceptor.csv")


# In[46]:


test_seq = pybedtools.BedTool.seq('chr1:11628-11927', fasta)
test_seq


# In[47]:


seq = Seq(test_seq)

print (seq.reverse_complement())


# In[14]:


reference_genome_file_path = '/mnt/data05/shared/pdutta_data/Human_Genome_Data/UCSC/genome_by_sequence/chr1.fa'
fasta = pybedtools.example_filename(reference_genome_file_path)


# In[125]:


a = a.sequence(fi=fasta)


# In[33]:


temp_seq = pybedtools.BedTool.seq('chr1:1-248956422', fasta)


# In[34]:


len(temp_seq)


# In[68]:


import sys
sys.getsizeof(temp_seq,'U6') 


# In[70]:


np_test_seq  = np.array(list(temp_seq))


# In[72]:


sys.getsizeof(np_test_seq,'U6') 


# In[115]:


test_str= "HGGVHGJBKJGNBJHBNBHJBNHGJHBHB"


# In[69]:


import numpy as np
test_str_new = np.array(test_str)


# In[61]:


test_str_new


# In[36]:


position = [[10,12], [16,20]]


# In[ ]:





# In[119]:


for i in position:
    print (i)
    sub_string = test_str[i[0]:i[1]+1]
    print(sub_string)
    substitute_str = "$"* int(i[1]+1-i[0])
    print (substitute_str)
    test_str = test_str.replace(sub_string, substitute_str)
    print (test_str)
line = re.sub('[$]', '', test_str)
print (line)


# In[51]:


# test_str_new[position[0][0]:position[0][1]] = 'NN'
test_str_new[0][0]


# In[12]:


test_str = test_str[0:position[0]]+ test_str[position[1]:]


# In[94]:


test_str = test_str.replace(test_str[0], '$')
test_str


# In[ ]:


for 


# In[54]:


check_sr= ""
for i in range(10):
    check_sr= check_sr+ '0'


# In[55]:


check_sr


# In[56]:


check_sr1= ""
check_sr1 = ["0" for i in range(10)  ]


# In[57]:


check_sr1


# In[58]:


checksr2= ''.join(["0" for i in range(10)])


# In[59]:


checksr2


# In[63]:


A = np.array([0,1,2,3,4,5,6])


# In[64]:


A


# In[66]:


A[0:2] = "N"


# In[77]:


for  i in range([10:12]):
    print(i)


# In[78]:


str1 = "0"*5


# In[79]:


str1


# In[84]:


p = range(1,4)


# In[85]:


p


# In[86]:


p[0]


# In[98]:


a_string = "Welme to Kite! CO"

index = a_string.find("co")

print(index)


# In[100]:


import re

string = 'i think that that that that student wrote there is not that right'
pattern = 'that'

print([match.start() for match in re.finditer(pattern, string)])


# In[104]:


pattern1= 'GT'


# In[106]:


Non_GT_donor = [match.start() for match in re.finditer(pattern1, temp_seq)]
len(Non_GT_donor)


# In[112]:


Non_GT_donor[0:10]


# In[8]:


list1= [10472, 10499, 10518, 10534, 10549, 10580, 10591, 10595, 10620, 10813, 10819, 10827, 10887, 10895, 10903]


# In[51]:


list1.sort()
print(list1)


# In[48]:


index = 1
while index < len(list1):
    print(f"comparing {list1[index-1]} with {list1[index]} nums list {list1}")
    if list1[index] - list1[index - 1] < 70:
        del list1[index]
    else:
        index += 1

print(list1)


# In[49]:


def remove_data(sorted_list, span):
    index = 1
    while index < len(sorted_list):
        #print(f"comparing {sorted_list[index-1]} with {sorted_list[index]} nums list {sorted_list}")
        if sorted_list[index] - sorted_list[index - 1] < span:
            del list1[index]
        else:
            index += 1
    
    print(list1)


# In[9]:


remove_data(list1, 70)


# In[10]:


def foo(elements, distance):
    #elements = sorted(elements) # sorting the user input
    new_elements = [elements[0]] # make a new list for output
    for element in elements[1:]: # Iterate over the remaining elements...
        if element - new_elements[-1] >= distance: 
            # this is the condition you described above
            new_elements.append(element)
            
    return new_elements


# In[11]:


chek_list = foo(list1, 70)


# In[12]:


chek_list


# In[ ]:




