import sys
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind, wilcoxon, mannwhitneyu
import time
import os


tissue_name = sys.argv[1]
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
# print(df.groupby([tissue_name]).size())
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




def stratified_sample(df, iteration, tissue_of_interest,  strata, size_percentage=None, seed=None, keep_index= True):
    population = len(df)
    print(population, size_percentage)
    #size = __smpl_size(population, size)
    size = int(size_percentage * population)
    print(size)
    tmp = df.loc[:,strata]
    tmp['size'] = 1
    #print (tmp)
    tmp_grpd = tmp.groupby(strata).count().reset_index()
    #print(tmp_grpd)
    tmp_grpd['samp_size'] = np.round(size/population * tmp_grpd['size']).astype(int)
    #print(tmp_grpd)

    # controlling variable to create the dataframe or append to it
    first = True 
    for i in range(len(tmp_grpd)):
        # query generator for each iteration
        #print ("i", i)
        qry=''
        for s in range(len(strata)):
            stratum = strata[s]
            #print("##",s, stratum)
            value = tmp_grpd.iloc[i][stratum]
            #print(value)
            n = tmp_grpd.iloc[i]['samp_size']
            #print("N", n)

            if type(value) == str:
                value = "'" + str(value) + "'"
                #print(value)
            
            if s != len(strata)-1:
                qry = qry + stratum + ' == ' + str(value) +' & '
            else:
                qry = qry + stratum + ' == ' + str(value)
            #print('Query', qry)
        
        # final dataframe
        if first:
            stratified_df = df.query(qry).sample(n=n, random_state=seed).reset_index(drop=(not keep_index))
            first = False
        else:
            tmp_df = df.query(qry).sample(n=n, random_state=seed).reset_index(drop=(not keep_index))
            stratified_df = pd.concat([stratified_df,tmp_df], ignore_index=True)
    print(tissue_of_interest)
    cat1_df = stratified_df[stratified_df['sample_class']==tissue_of_interest]
    cat2_df = stratified_df[stratified_df['sample_class']!=tissue_of_interest]
    for index, row in df_1s_t_test.iterrows():
        #print(row['Transcript_ID'])
        t, pt1 = ttest_ind(cat1_df[row['Transcript_ID']], cat2_df[row['Transcript_ID']], equal_var=False) ## One-tailed t-test
        t, pt2 = ttest_ind(cat1_df[row['Transcript_ID']], cat2_df[row['Transcript_ID']], equal_var=True) ## Two-tailed t-test
        t, pw2 = mannwhitneyu(cat1_df[row['Transcript_ID']], cat2_df[row['Transcript_ID']], alternative = 'two-sided') ## Two-tailed Wilcox test
        t, pw1g = mannwhitneyu(cat1_df[row['Transcript_ID']], cat2_df[row['Transcript_ID']], alternative = 'greater') ## One-tailed geater Wilcox test
        t, pw1s = mannwhitneyu(cat1_df[row['Transcript_ID']], cat2_df[row['Transcript_ID']], alternative = 'less') ## One-tailed geater Wilcox test
        
        
        
        df_1s_t_test.loc[index, ['ITER_{}'.format(iter)]]= pt1
        df_2s_t_test.loc[index, ['ITER_{}'.format(iter)]]= pt2
        df_1sg_wilcox.loc[index, ['ITER_{}'.format(iter)]]= pw1g
        df_1ss_wilcox.loc[index, ['ITER_{}'.format(iter)]]= pw1s
        df_2s_wilcox.loc[index, ['ITER_{}'.format(iter)]]= pw2
    return None





start_time = time.time()
for i in range(no_of_iteration):
    print(i)
    stratified_sample(df=df, iteration =i, tissue_of_interest = tissue_name, strata=['sample_class'], size_percentage =0.7)
print("--- %s seconds ---" % (time.time() - start_time))


df_1s_t_test.set_index('Transcript_ID', inplace=True)
df_2s_t_test.set_index('Transcript_ID', inplace=True)
df_1sg_wilcox.set_index('Transcript_ID', inplace=True)
df_1ss_wilcox.set_index('Transcript_ID', inplace=True)
df_2s_wilcox.set_index('Transcript_ID', inplace=True)
#print (df_p_values)

df_1s_t_test.to_csv(OUTPUT_PATH +"/"+ tissue_name+"_1TailedTtest_iter="+ str(no_of_iteration)+"_sample_size="+str(sample_size)+"_iteration_p.csv", sep=",")
df_2s_t_test.to_csv(OUTPUT_PATH +"/"+ tissue_name+"_2TailedTtest_iter="+ str(no_of_iteration)+"_sample_size="+str(sample_size)+"_iteration_p.csv", sep=",")
df_1sg_wilcox.to_csv(OUTPUT_PATH +"/"+ tissue_name+"_1TailedWilCoxGreater_iter="+ str(no_of_iteration)+"_sample_size="+str(sample_size)+"_iteration_p.csv", sep=",")
df_1ss_wilcox.to_csv(OUTPUT_PATH +"/"+ tissue_name+"_1TailedWilCoxSmaller_iter="+ str(no_of_iteration)+"_sample_size="+str(sample_size)+"_iteration_p.csv", sep=",")
df_2s_wilcox.to_csv(OUTPUT_PATH +"/"+ tissue_name+"_2TailedWilCox_iter="+ str(no_of_iteration)+"_sample_size="+str(sample_size)+"_iteration_p.csv", sep=",")

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

t1 = df_1s_t_test.apply(count, axis = 1)
temp_df_empirical['1SidedTTest_empirical_p_value'] = t1.values
print()
temp_df_empirical['2SidedTTest_empirical_p_value'] = df_2s_t_test.apply(count, axis = 1).values
print(temp_df_empirical)
temp_df_empirical['1SidedWilCoxGreater_empirical_p_value'] = df_1sg_wilcox.apply(count, axis = 1).values
print(temp_df_empirical)
temp_df_empirical['1SidedWilCoxSmaller_empirical_p_value'] = df_1ss_wilcox.apply(count, axis = 1).values
print(temp_df_empirical)
temp_df_empirical['2SidedWilCox_empirical_p_value'] = df_2s_wilcox.apply(count, axis = 1).values
print(temp_df_empirical)

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