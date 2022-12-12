import sys
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind, wilcoxon, mannwhitneyu
import time
import os


tissue_name = sys.argv[1]
tissue_name  = tissue_name.strip()
print(tissue_name)
no_of_iteration = 500
sample_size = 0.7
file_path = "/data/projects/Tissue_specific_data/subsample/Data/all_input_stratified.feather"
OUTPUT_PATH = "/data/projects/Tissue_specific_data/subsample/results/statistical_tests/"+tissue_name
OUTPUT_FILE_NAME_ITERATION = tissue_name+"_iter="+ str(no_of_iteration)+"_sample_size="+str(sample_size)+"_iteration_p.csv"
OUTPUT_FILE_NAME_EMPIRICAL = tissue_name+"_iter="+ str(no_of_iteration)+"_sample_size="+str(sample_size)+"_empirical_p.csv"

if not os.path.exists(OUTPUT_PATH):
    os.makedirs(OUTPUT_PATH)

df= pd.read_feather(file_path)

print(df.shape)
#print(list(df.columns))

print (df.columns[2:])
print (len(df.columns[2:]))
print(df.groupby(['sample_class']).size())
# groups_name = list(df.groupby([tissue_name]).groups.keys())
# print(groups_name)

df_1s_t_test = pd.DataFrame()
df_2s_t_test= pd.DataFrame()
df_1sg_wilcox= pd.DataFrame()
df_1ss_wilcox= pd.DataFrame()
df_2s_wilcox= pd.DataFrame()

df_1s_t_test['Transcript_ID'] = df.columns[2:]
df_2s_t_test['Transcript_ID'] = df.columns[2:]
df_1sg_wilcox['Transcript_ID'] = df.columns[2:]
df_1ss_wilcox['Transcript_ID'] = df.columns[2:]
df_2s_wilcox['Transcript_ID'] = df.columns[2:]



start_time = time.time()
for i in range(no_of_iteration):
    print('#$@',i)
    #stratified_sample(df=df, iteration =i, tissue_of_interest = tissue_name, strata=['sample_class'], size_percentage =0.7)
    stratified_df = df.groupby('sample_class', group_keys=False).apply(lambda x: x.sample(frac=0.7)).reset_index(drop=True)
    print(tissue_name)
    cat1_df = stratified_df[stratified_df['sample_class']==tissue_name]
    cat2_df = stratified_df[stratified_df['sample_class']!=tissue_name]
    print(cat1_df.shape , cat2_df.shape)
    one_tailed_p_value_list=[]
    two_tailed_p_value_list=[]
    two_tailed_Wil_p_value_list=[]
    one_tailed_Wil_smaller_p_value_list=[]
    one_tailed_Wil_greater_p_value_list=[]
    for index, row in df_1s_t_test.iterrows():
        if ((index%10000)==0):
            print("##", index)
        # print("##",row['Transcript_ID'])
        # print("$$#",cat1_df[row['Transcript_ID']].shape)
        # print("$$",cat2_df[row['Transcript_ID']].shape)
        one_tailed_p_value_list.append(ttest_ind(cat1_df[row['Transcript_ID']], cat2_df[row['Transcript_ID']], equal_var=False)[1]) ## One-tailed t-test
        two_tailed_p_value_list.append(ttest_ind(cat1_df[row['Transcript_ID']], cat2_df[row['Transcript_ID']], equal_var=True)[1]) ## Two-tailed t-test
        two_tailed_Wil_p_value_list.append(mannwhitneyu(cat1_df[row['Transcript_ID']], cat2_df[row['Transcript_ID']], alternative = 'two-sided')[1]) ## Two-tailed Wilcox test
        one_tailed_Wil_greater_p_value_list.append(mannwhitneyu(cat1_df[row['Transcript_ID']], cat2_df[row['Transcript_ID']], alternative = 'greater')[1]) ## One-tailed geater Wilcox test
        one_tailed_Wil_smaller_p_value_list.append(mannwhitneyu(cat1_df[row['Transcript_ID']], cat2_df[row['Transcript_ID']], alternative = 'less')[1]) ## One-tailed geater Wilcox test
    


    df_1s_t_test['ITER_{}'.format(i)]= one_tailed_p_value_list
    df_2s_t_test['ITER_{}'.format(i)]= two_tailed_p_value_list
    df_1sg_wilcox['ITER_{}'.format(i)]= one_tailed_Wil_greater_p_value_list
    df_1ss_wilcox['ITER_{}'.format(i)]= one_tailed_Wil_smaller_p_value_list
    df_2s_wilcox['ITER_{}'.format(i)]= two_tailed_Wil_p_value_list
    #print (df_p_values)
    
    print("--- %s seconds ---" % (time.time() - start_time)) 
    if((i+1)%100==0):
        df_1s_t_test.to_csv(OUTPUT_PATH +"/"+ tissue_name+"_1TailedTtest_iter="+ str(i)+"_sample_size="+str(sample_size)+"_iteration_p.csv", sep=",")
        df_2s_t_test.to_csv(OUTPUT_PATH +"/"+ tissue_name+"_2TailedTtest_iter="+ str(i)+"_sample_size="+str(sample_size)+"_iteration_p.csv", sep=",")
        df_1sg_wilcox.to_csv(OUTPUT_PATH +"/"+ tissue_name+"_1TailedWilCoxGreater_iter="+ str(i)+"_sample_size="+str(sample_size)+"_iteration_p.csv", sep=",")
        df_1ss_wilcox.to_csv(OUTPUT_PATH +"/"+ tissue_name+"_1TailedWilCoxSmaller_iter="+ str(i)+"_sample_size="+str(sample_size)+"_iteration_p.csv", sep=",")
        df_2s_wilcox.to_csv(OUTPUT_PATH +"/"+ tissue_name+"_2TailedWilCox_iter="+ str(i)+"_sample_size="+str(sample_size)+"_iteration_p.csv", sep=",")
        
    
    
    


start_time = time.time()

def count(df):
    #n=df.shape[1]
    #print(df)
    
    #print(no_of_iteration)
    r = np.sum(df >= 0.05)
    df['empirical_p_value'] = (r+1) /(no_of_iteration+1)
    #df['r'] = np.sum(df == 'r')
    #print(df['empirical_p_value'])
    #print(temp_list)
    temp_list = df['empirical_p_value'].tolist()
    #print(temp_list)
    # with open(OUTPUT_PATH +"/"+'your_file.txt', 'w') as f:
    #     np.savetxt(f, temp_list)
    return temp_list

temp_df_empirical = pd.DataFrame()
temp_df_empirical['Transcript_ID'] = df.columns[2:]

temp_df_empirical['1SidedTTest_empirical_p_value'] = df_1s_t_test.apply(count, axis = 1).values
temp_df_empirical['2SidedTTest_empirical_p_value'] = df_2s_t_test.apply(count, axis = 1).values
temp_df_empirical['1SidedWilCoxGreater_empirical_p_value'] = df_1sg_wilcox.apply(count, axis = 1).values
temp_df_empirical['1SidedWilCoxSmaller_empirical_p_value'] = df_1ss_wilcox.apply(count, axis = 1).values
temp_df_empirical['2SidedWilCox_empirical_p_value'] = df_2s_wilcox.apply(count, axis = 1).values

temp_df_empirical.to_csv(OUTPUT_PATH +"/"+OUTPUT_FILE_NAME_EMPIRICAL, sep=",")


print("--- %s seconds ---" % (time.time() - start_time))
# df_mean = df.groupby(tissue_name).mean(numeric_only=True)
# index_list = list(df_mean.index.values)
# df_mean.loc['Fold_change_'+index_list[0]+'/'+index_list[1]] = df_mean.iloc[0]/df_mean.iloc[1]
# df_mean.loc['Fold_change_'+index_list[1]+'/'+index_list[0]] = df_mean.iloc[1]/df_mean.iloc[0]
# df_mean
# temp_df_empirical['Fold_change_no_tissue/yes_tissue'] = df_mean.loc['Fold_change_no_tissue/yes_tissue'].values[:]
# temp_df_empirical['Fold_change_yes_tissue/no_tissue'] = df_mean.loc['Fold_change_yes_tissue/no_tissue'].values[:]

# temp_df_empirical.iloc[: , -3:].to_csv(OUTPUT_PATH +"/"+OUTPUT_FILE_NAME_EMPIRICAL, sep=",")