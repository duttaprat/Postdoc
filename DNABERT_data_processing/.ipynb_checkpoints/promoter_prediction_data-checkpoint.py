import os, tarfile
import re
import pandas as pd
import numpy as np
import zipfile, gzip
import pandas as pd
from sklearn.model_selection import train_test_split




promoter_data_path = "/mnt/data05/shared/pdutta_data/Human_Genome_Data/chromosome_wise_sequence_2000bp/Promoters/"
non_promoter_data_path = "/mnt/data05/shared/pdutta_data/Human_Genome_Data/chromosome_wise_sequence_2000bp/Non_promoters/"
data_path= "/mnt/data05/shared/pdutta_data/DNABERT_data/Promoter_prediction/2000BP_promoter_prediction_data/1:1/"





def seq2kmer(seq, k):
    kmer = [seq[x:x+k] for x in range(len(seq)+1-k)]
    #print (kmer)
    kmers = " ".join(kmer)
    #print (kmers)
    return kmers




# In[37]:


def df2kmer(dataframe, name):
    temp_df= pd.DataFrame(columns= ['sequence','label' ])
    for index, row in dataframe.iterrows():
        temp_df.loc[index,['sequence']]= seq2kmer(row['sequence'], 6)
        temp_df.loc[index,['label']]= row['label']
        #print(temp)
        # if(index==10):
        #     break
    temp_df.to_csv(data_path+name, index=False,  sep='\t')
    #return temp_df


# In[6]:


def create_csv(path):
    file_names = os.listdir(path)
    for i in file_names:
        print (i)
    li = []
    for filename in file_names:
        df = pd.read_csv(path+filename, index_col=None, header=0)
        li.append(df)
    frame = pd.concat(li, axis=0, ignore_index=True)
    print ("New", frame.shape)
    return frame






# In[7]:


promoter_csv = create_csv(promoter_data_path)
non_promoter_csv = create_csv(non_promoter_data_path)







promoter_sequence = promoter_csv[['sequence']]
promoter_sequence['label']= 1
print (promoter_sequence.head())
non_promoter_sequence = non_promoter_csv[['0']]
non_promoter_sequence = non_promoter_sequence.rename(columns={'0': 'sequence'})
non_promoter_sequence['label']= 0
print (non_promoter_sequence.head())


# In[14]:


promoter_sequence = promoter_sequence[~promoter_sequence.sequence.str.contains("N")]
promoter_sequence['sequence'] = promoter_sequence['sequence'].str.upper()
non_promoter_sequence = non_promoter_sequence[~non_promoter_sequence.sequence.str.contains("N")]
non_promoter_sequence['sequence'] = non_promoter_sequence['sequence'].str.upper()


# In[15]:
print (promoter_sequence.head())
print (non_promoter_sequence.head())



print (promoter_sequence.shape, non_promoter_sequence.shape)


non_promoter_rest, non_promoter_data = train_test_split(non_promoter_sequence, test_size=0.046)
non_promoter_rest.reset_index(drop=True, inplace=True)
non_promoter_data.reset_index(drop=True, inplace=True)
print (non_promoter_rest.shape , non_promoter_data.shape)


print("Promoter and Non-promoter sizes", promoter_sequence.shape, non_promoter_data.shape)
# In[32]:



print (len(promoter_sequence['sequence'].loc[1]), len(non_promoter_data['sequence'].loc[1]))
non_promoter_rest.head()


# In[33]:


all_data = pd.concat([promoter_sequence, non_promoter_data], ignore_index=True)
train, test = train_test_split(all_data, test_size=0.10)
#dev, test= train_test_split(test, test_size=0.5)
train.reset_index(drop=True, inplace=True)
test.reset_index(drop=True, inplace=True)
#dev.reset_index(drop=True, inplace=True)


# In[34]:


print (train.shape, test.shape)

# In[35]:



df2kmer(test, "dev.tsv")
df2kmer(train, "train.tsv")
#df2kmer(dev,"dev.tsv")


# In[ ]:


#df2kmer(non_promoter_rest, "non_promoter_rest.tsv")


# In[ ]:




