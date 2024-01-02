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



non_coding_region= "donor"
cancer_type= "Brain"


base_dir = "/data/private/pdutta/DNABERT_data/Genetic_variant/Brain/{}".format(non_coding_region)
df_donor = pd.read_csv("/data/private/pdutta/DNABERT_data/Splice_sites/Donor_raw_files/Unique_Donor_RAW_Combined_BED.tsv", sep="\t")


with open("/data/projects/GDC_Cancer_Wise/{}/Data/Sequence_generation/raw_{}_vcf_data.pkl".format(cancer_type, non_coding_region), "rb") as file:
    loaded_dictionary = pickle.load(file)



def log_oods_ratio(p1,p2):
    term1 = math.log2(p1/(1-p1))
    term2 = math.log2(p2/(1-p2))
    return (term1-term2)
def calculate_diff_probab(p1,p2):
    return((p2-p1)*max(p1,p2))


# In[41]:


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
        df_temp.rename(columns={'Acceptor_start': 'donor_start', 'Acceptor_end': 'donor_end'}, inplace=True)
        patient_dfs.append(df_temp)
        #input()
    break


# In[43]:


# Concatenate the patient dataframes and drop duplicates to get all unique regions across the patients
unique_patient_regions = pd.concat(patient_dfs).drop_duplicates(subset=['chr', 'Transcript_ID', 'donor_start', 'donor_end'])
print(unique_patient_regions)




# Initialize the dictionaries
data_ref, data_alt, data_log_odd, data_score_change = {}, {}, {},  {}

# Loop over each region in the unique_patient_regions
for _, region in unique_patient_regions.iterrows():
    key = f'{region.chr}_{region["{}_start".format(non_coding_region)]}_{region["{}_end".format(non_coding_region)]}_{region.varinat_start}_{region.variant_end}'
    data_ref[key] = []
    data_alt[key] = []
    data_log_odd[key] = []
    data_score_change[key] = []

# Loop through each patient
for df in patient_dfs:
    print(df.shape)
    for _, region in unique_patient_regions.iterrows():
        for data, score_column in zip([data_ref, data_alt, data_log_odd, data_score_change], ['Ref_probab', 'Alt_probab', 'LogOddRatio', 'ScoreChange']):
            score_series = df.loc[(df['chr'] == region.chr) & (df['{}_start'.format(non_coding_region)] == region["{}_start".format(non_coding_region)]) & (df['varinat_start'] == region.varinat_start), score_column]
            if score_series.empty:
                score = np.nan
            else:
                score = score_series.max()
            data[f'{region.chr}_{region["{}_start".format(non_coding_region)]}_{region["{}_end".format(non_coding_region)]}_{region.varinat_start}_{region.variant_end}'].append(score)
    print(df)




result_df_ref = pd.DataFrame(data_ref)
result_df_ref.insert(0, "Patient_ID", patient_id)


result_df_alt = pd.DataFrame(data_alt)
result_df_alt.insert(0, "Patient_ID", patient_id)


result_df_log_odd = pd.DataFrame(data_log_odd)
result_df_log_odd.insert(0, "Patient_ID", patient_id)



result_df_score_change = pd.DataFrame(data_score_change)
result_df_score_change.insert(0, "Patient_ID", patient_id)



result_df_ref.to_csv("/data/projects/GDC_Cancer_Wise/Brain/Data/DNABERT_result/{}/df_ref_score.tsv".format(non_coding_region), sep="\t", index=False)
result_df_alt.to_csv("/data/projects/GDC_Cancer_Wise/Brain/Data/DNABERT_result/{}/df_alt_score.tsv".format(non_coding_region), sep="\t", index=False)
result_df_log_odd.to_csv("/data/projects/GDC_Cancer_Wise/Brain/Data/DNABERT_result/{}/df_log_odd_score.tsv", sep="\t".format(non_coding_region), index=False)
result_df_score_change.to_csv("/data/projects/GDC_Cancer_Wise/Brain/Data/DNABERT_result/{}/df_log_odd_score.tsv", sep="\t".format(non_coding_region), index=False)