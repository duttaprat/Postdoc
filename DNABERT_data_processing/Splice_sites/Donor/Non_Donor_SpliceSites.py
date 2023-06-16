#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os, tarfile
import re
import pandas as pd
import numpy as np
import zipfile, gzip
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.utils import shuffle


# In[2]:


folder_path = "/home/pdutta/Data/Splicesites_data/Non_donor_splice_sites/Data"


# ## Creation of Non-acceptor Set

# In[4]:


def concatenate_csv_files(folder_path):
    # Find all csv files in the folder
    csv_files = [f for f in os.listdir(folder_path) if f.endswith('.txt')]

    # Read each csv file and store them in a list
    dataframes = []
    for csv in csv_files:
        dataframes.append(pd.read_csv(os.path.join(folder_path, csv),  names=['Sequence']))
        print(csv," is loaded")

    # Concatenate all dataframes
    combined_df = pd.concat(dataframes, ignore_index=True)
    return combined_df


# In[ ]:


concat_negative_df = concatenate_csv_files(folder_path)
print(concat_negative_df)
print(concat_negative_df.shape)

concat_negative_df = concat_negative_df.drop_duplicates().reset_index(drop=True)
print(concat_negative_df)
print(concat_negative_df.shape)

concat_negative_df = concat_negative_df[~concat_negative_df['Sequence'].str.contains('N')]
concat_negative_df['Label'] =0
print(concat_negative_df.shape)




concat_negative_df.to_csv("/home/pdutta/Data/Splicesites_data/Non_donor_concat.tsv", sep="\t", index=False)




