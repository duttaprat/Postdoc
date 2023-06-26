import os
import pandas as pd
import numpy as np



file_path = "/data/projects/Tissue_specific_data/subsample/Data/all_input_stratified.feather"
output_path = "/data/projects/Tissue_specific_data/Tissue_specific_transcript"


df= pd.read_feather(file_path)

tissue_names_list  = list(set(df['sample_class']))
df_mean = pd.DataFrame()
for i in tissue_names_list:
    df_temp = df[df['sample_class']==i] 
    #print(len(df_temp))
    df_mean["{}".format(i)] = df_temp.mean(axis=0)


print(df_mean)


tissue_list=  tissue_names_list
threshold = 10
df_store = pd.DataFrame()
for tissue_name in tissue_list:
    #input()
    print(tissue_name)
    max_val = int(df_mean[tissue_name].max())
    print(max_val)
    folder_path = output_path + "/"+ tissue_name
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
    for threshold in range(10,max_val):
        if (threshold%1000==0):
            print(threshold)
        df_temp = df_mean[df_mean[tissue_name]>threshold]
        #print(df_temp.head())
        df_other = df_temp.drop([tissue_name], axis=1)
        df_specific =  df_temp[(df_other<threshold).all(axis=1)]
        #print(df_specific.shape)
        length = len(df_specific)
        #print('df_thr='+str(threshold)+'_num='+str(len(df_specific)))
        if (len(df_specific)>0):
            if not (df_specific.equals(df_store)):
                df_store  = df_specific
                df_specific.to_csv(folder_path+'/df_thr='+str(threshold)+'_num='+str(len(df_specific))+'_transcripts_list.csv', sep= ",")

