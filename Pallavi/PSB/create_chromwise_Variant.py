import pandas as pd
import gzip
from tqdm import tqdm
import os
import io


# In[2]:
file_path = '/home/pdutta/Data/DBSNP/dbSnp155Details.tab.gz'
output_path = "/home/pdutta/Data/DBSNP/Chromosome_wise"




total_size = os.path.getsize(file_path)
print(total_size)

column_names = ['rs_ID', 'freqSourceCount', 'alleleCounts', 'alleleTotals', 'soTermCount', 'soTerms', 'clinVarCount', 'clinVarAccs', 'clinVarSigs', 'submitterCount', 'submitters', 'pubMedIdCount' ,'pubMedIds' ]

# initialize progress bar
progress_bar = tqdm(total=total_size, unit='B', unit_scale=True)


# In[7]:


def read_in_chunks(file_object, chunk_size=1024 * 1024):  # read in chunks of 1 MB
    while True:
        data = file_object.read(chunk_size)
        if not data:
            break
        progress_bar.update(len(data))
        yield data


# In[9]:


buffer = ''
with gzip.open(file_path, 'rt') as f:
    for chunk in read_in_chunks(f):
        buffer += chunk
    df = pd.read_csv(io.StringIO(buffer), delimiter='\t')

progress_bar.close()


print(df)
print(df.columns)



# In[ ]:






