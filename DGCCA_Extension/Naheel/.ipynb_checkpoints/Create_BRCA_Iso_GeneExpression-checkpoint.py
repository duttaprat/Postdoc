import pandas as pd


df_clinical = pd.read_csv("/home/pdutta/Data/DGCCA/clinical_files/BRCA/nationwidechildrens.org_clinical_patient_brca.txt", sep='\t')
df_clinical = df_clinical.iloc[2: , :]
df_clinical['bcr_patient_barcode'] = df_clinical['bcr_patient_barcode'].astype(str)+'-01' 



print(len(set(df_clinical['bcr_patient_barcode'].to_list())))


# In[28]:


df_PAN_isoform = pd.read_csv("/home/pdutta/Data/DGCCA/DGCCA_extension/tcga_RSEM_isoform_fpkm", sep="\t") 


# In[29]:


print(df_PAN_isoform)
print(df_clinical)


# In[34]:


set_PAN = set(df_PAN_isoform)
set_BRCA_clinical = set(df_clinical['bcr_patient_barcode'].to_list())

print(set_PAN)
print(set_BRCA_clinical)
# In[38]:


common_list = list(set_PAN.intersection(set_BRCA_clinical))
print("FGFGF", len(common_list))


# In[40]:


common_list.insert(0, 'sample')
print(common_list)

# In[41]:


print(df_PAN_isoform[common_list])


df_isoform_BRCA = df_PAN_isoform[common_list]

df_isoform_BRCA.T.to_csv("/home/pdutta/Data/DGCCA/DGCCA_extension/tcga_BRCA_RSEM_isoform_fpkm.csv")
# In[ ]:




