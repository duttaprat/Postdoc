
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
import time
import os


tissue_name = "Kidney"
no_of_iteration = 500
sample_size = 0.8
file_path = "/data/projects/shared_data/tissue_specific_data/subsample/all_input.feather"
OUTPUT_PATH = "/data/projects/shared_data/tissue_specific_data/subsample/results/"+tissue_name
OUTPUT_FILE_NAME_ITERATION = tissue_name+"_"+ str(no_of_iteration)+"_"+str(sample_size)+"_iteration_p.csv"
OUTPUT_FILE_NAME_EMPIRICAL = tissue_name+"_"+ str(no_of_iteration)+"_"+str(sample_size)+"_empirical_p.csv"

if not os.path.exists(OUTPUT_PATH):
    os.makedirs(OUTPUT_PATH)

df= pd.read_feather(file_path)

print(df.shape)
#print(list(df.columns))

print (df.columns[31:])
print (len(df.columns[31:]))
print(df.groupby([tissue_name]).size())
groups_name = list(df.groupby([tissue_name]).groups.keys())
print(groups_name)

df_p_values = pd.DataFrame()
df_p_values['Transcript_ID'] = df.columns[31:]


def stratified_sample(df, iter, strata, size=None,  seed=None, keep_index= True):
    population = len(df)
    #print(population, size)
    size= np.round(population * size)
    tmp = df.loc[:,strata]
    tmp['size'] = 1
    tmp_grpd = tmp.groupby(strata).count().reset_index()

    tmp_grpd['samp_size'] = np.round(size/population * tmp_grpd['size']).astype(int)

    # controlling variable to create the dataframe or append to it
    first = True 
    for i in range(len(tmp_grpd)):
        # query generator for each iteration
        qry=''
        for s in range(len(strata)):
            stratum = strata[s]
            value = tmp_grpd.iloc[i][stratum]
            n = tmp_grpd.iloc[i]['samp_size']

            if type(value) == str:
                value = "'" + str(value) + "'"
            
            if s != len(strata)-1:
                qry = qry + stratum + ' == ' + str(value) +' & '
            else:
                qry = qry + stratum + ' == ' + str(value)
        
        # final dataframe
        if first:
            stratified_df = df.query(qry).sample(n=n, random_state=seed).reset_index(drop=(not keep_index))
            first = False
        else:
            tmp_df = df.query(qry).sample(n=n, random_state=seed).reset_index(drop=(not keep_index))
            stratified_df = stratified_df.append(tmp_df, ignore_index=True)
    #print(stratified_df.head(), tissue_name, groups_name)
    cat1_df = stratified_df[stratified_df[tissue_name]==groups_name[0]]
    cat2_df = stratified_df[stratified_df[tissue_name]==groups_name[1]]
    print(iter)
    for index, row in df_p_values.iterrows():
        #print(row['Transcript_ID'])
        t, p = ttest_ind(cat1_df[row['Transcript_ID']], cat2_df[row['Transcript_ID']], equal_var=False)
        df_p_values.loc[index, ['ITER_{}'.format(iter)]]= p
    return None

def __smpl_size(population, size):
    '''
    A function to compute the sample size. If not informed, a sampling 
    size will be calculated using Cochran adjusted sampling formula:
        cochran_n = (Z**2 * p * q) /e**2
        where:
            - Z is the z-value. In this case we use 1.96 representing 95%
            - p is the estimated proportion of the population which has an
                attribute. In this case we use 0.5
            - q is 1-p
            - e is the margin of error
        This formula is adjusted as follows:
        adjusted_cochran = cochran_n / 1+((cochran_n -1)/N)
        where:
            - cochran_n = result of the previous formula
            - N is the population size
    '''
    if size is None:
        cochran_n = np.round(((1.96)**2 * 0.5 * 0.5)/ 0.02**2)
        n = np.round(cochran_n/(1+((cochran_n -1) /population)))
    elif size >= 0 and size < 1:
        n = np.round(population * size)
    elif size < 0:
        raise ValueError('Parameter "size" must be an integer or a proportion between 0 and 0.99.')
    elif size >= 1:
        n = size
    return n




start_time = time.time()
for i in range(no_of_iteration):
    stratified_sample(df=df, iter =i, strata=[tissue_name], size=sample_size)
print("--- %s seconds ---" % (time.time() - start_time))


df_p_values.set_index('Transcript_ID', inplace=True)
print (df_p_values)

df_p_values.to_csv(OUTPUT_PATH +"/"+OUTPUT_FILE_NAME_ITERATION, sep=",")
start_time = time.time()
def count(df):
    #n=df.shape[1]
    #print(df)
    threshold_values=[0.05]
    for threshold in threshold_values:
        r = np.sum(df >= threshold)
        df['>={}'.format(threshold)] = (r+1) /(no_of_iteration+1)
    #df['r'] = np.sum(df == 'r')
    return df

temp_df_empirical = df_p_values.apply(count, axis = 1)
print("--- %s seconds ---" % (time.time() - start_time))
df_mean = df.groupby(tissue_name).mean(numeric_only=True)
index_list = list(df_mean.index.values)
df_mean.loc['Fold_change_'+index_list[0]+'/'+index_list[1]] = df_mean.iloc[0]/df_mean.iloc[1]
df_mean.loc['Fold_change_'+index_list[1]+'/'+index_list[0]] = df_mean.iloc[1]/df_mean.iloc[0]
df_mean
temp_df_empirical['Fold_change_no_tissue/yes_tissue'] = df_mean.loc['Fold_change_no_tissue/yes_tissue'].values[:]
temp_df_empirical['Fold_change_yes_tissue/no_tissue'] = df_mean.loc['Fold_change_yes_tissue/no_tissue'].values[:]

temp_df_empirical.iloc[: , -3:].to_csv(OUTPUT_PATH +"/"+OUTPUT_FILE_NAME_EMPIRICAL, sep=",")