import os, glob
import numpy as np
import pandas as pd
import math
import pickle
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy import stats
from statsmodels.stats.weightstats import ztest
from IPython.display import display, HTML
from multiprocessing import Pool




# In[2]:


base_dir = "/data/private/pdutta/DNABERT_data/Genetic_variant/Brain/Acceptor"


# In[3]:


df_acceptor = pd.read_csv("/data/private/pdutta/DNABERT_data/Splice_sites/Acceptor_raw_files/Unique_Acceptor_RAW_Combined_BED.tsv", sep="\t")
#df_acceptor


# In[4]:


with open("/data/projects/GDC_Cancer_Wise/Brain/Data/raw_acceptor_vcf_data.pkl", "rb") as file:
    loaded_dictionary = pickle.load(file)


# In[5]:


def log_oods_ratio(p1,p2):
    term1 = math.log2(p1/(1-p1))
    term2 = math.log2(p2/(1-p2))
    return (term1-term2)
def calculate_diff_probab(p1,p2):
    return((p2-p1)*max(p1,p2))


# In[7]:


patient_dfs= []
patient_id = []
for dirpath, dirnames, filenames in os.walk(base_dir):
    #print(dirpath, dirnames, filenames)
    for dir_name in dirnames:
        dir_path = dirpath +"/" + dir_name
        print(dir_path)
        patient_id.append(dir_name)
        df_temp = loaded_dictionary[dir_name]
        concat_probab = np.load(dir_path+"/Prediction_result/pred_results.npy")
        print(len(concat_probab))
        Reference_probab = concat_probab[::2]
        Alternative_probab = concat_probab[1::2]
        df_temp['Ref_probab']=Reference_probab
        df_temp['Alt_probab']=Alternative_probab
        df_temp['ScoreChange'] = df_temp.apply(lambda row: calculate_diff_probab(row['Ref_probab'], row['Alt_probab']), axis=1)
        df_temp['LogOddRatio'] = df_temp.apply(lambda row: log_oods_ratio(row['Ref_probab'], row['Alt_probab']), axis=1)
        #display(df_temp)
        patient_dfs.append(df_temp)
        #input()
    break


# In[8]:


# Concatenate the patient dataframes and drop duplicates to get all unique regions across the patients
unique_patient_regions = pd.concat(patient_dfs).drop_duplicates(subset=['chr', 'Transcript_ID', 'Acceptor_start', 'Acceptor_end'])
unique_patient_regions


# In[9]:


unique_patient_regions[unique_patient_regions['Acceptor_start']==28452437]


# In[ ]:


# Step 3: Create a new dataframe with three rows for each patient, with columns as the acceptor coordinates
data_ref = {f'{region.chr}_{region.Acceptor_start}_{region.Acceptor_end}_{region.varinat_start}_{region.variant_end}': [] for _, region in unique_patient_regions.iterrows()}
data_alt = {f'{region.chr}_{region.Acceptor_start}_{region.Acceptor_end}_{region.varinat_start}_{region.variant_end}': [] for _, region in unique_patient_regions.iterrows()}
data_log_odd = {f'{region.chr}_{region.Acceptor_start}_{region.Acceptor_end}_{region.varinat_start}_{region.variant_end}': [] for _, region in unique_patient_regions.iterrows()}

# Loop through each patient and for each region, get the Alt_probab score for that region for that patient
for df in patient_dfs:
    print(df.shape)
    #print(df.columns)
    for _, region in unique_patient_regions.iterrows():
        for data, score_column in zip([data_ref, data_alt, data_log_odd], ['Ref_probab', 'Alt_probab', 'LogOddRatio']):
            score_series = df.loc[(df['chr'] == region.chr) & (df['Acceptor_start'] == region.Acceptor_start) & (df['varinat_start'] == region.varinat_start) , score_column]
            if score_series.empty:
                score = np.nan
            else:
                score = score_series.max()
                #print(score_series , score)
                #input()
            data[f'{region.chr}_{region.Acceptor_start}_{region.Acceptor_end}_{region.varinat_start}_{region.variant_end}'].append(score)


# In[ ]:


result_df_ref = pd.DataFrame(data_ref)
result_df_ref.insert(0, "Patient_ID", patient_id)
#result_df_ref


# In[ ]:


result_df_alt = pd.DataFrame(data_alt)
result_df_alt.insert(0, "Patient_ID", patient_id)
#result_df_alt


# In[ ]:


result_df_log_odd = pd.DataFrame(data_log_odd)
result_df_log_odd.insert(0, "Patient_ID", patient_id)
#result_df_log_odd


# In[ ]:


result_df_ref.to_csv("/data/projects/GDC_Cancer_Wise/Brain/Data/DNABERT_result/Acceptor/df_ref_score.tsv", sep="\t", index=False)


# In[ ]:


result_df_alt.to_csv("/data/projects/GDC_Cancer_Wise/Brain/Data/DNABERT_result/Acceptor/df_alt_score.tsv", sep="\t", index=False)


# In[ ]:


result_df_log_odd.to_csv("/data/projects/GDC_Cancer_Wise/Brain/Data/DNABERT_result/Acceptor/df_log_odd_score.tsv", sep="\t", index=False)


# In[ ]:





# In[ ]:




