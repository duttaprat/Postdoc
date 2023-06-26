
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from multiprocessing import Pool
import time
import multiprocessing
#multiprocessing.set_start_method('fork')

tissue_name = "Bladder"
no_of_iteration = 5
sample_size = 0.8
file_path = "/data/projects/shared_data/collab_data/subsample/all_input.feather"


df= pd.read_feather(file_path)

print(df.shape)
#print(list(df.columns))
print(df.groupby([tissue_name]).size())
groups_name = list(df.groupby([tissue_name]).groups.keys())
print(groups_name)

def stratified_sample(df, strata, size=None, seed=None, keep_index= True):
    population = len(df)
    #print(population, size)
    size = __smpl_size(population, size)
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
    
    return stratified_df

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


# NUM_WORKERS = 5

# def prepare_and_save(i):
#     print(i)
#     return stratified_sample(df=df, strata=[tissue_name], size=sample_size)

# p = Pool(NUM_WORKERS)
start_time = time.time()
# output_ = p.map(prepare_and_save, np.arange(5))


# pool = Pool(processes=5)              # start 4 worker processes
# result = pool.apply_async(stratified_sample(df=df, strata=[tissue_name], size=sample_size), [5])    # evaluate "f(10)" asynchronously
# #print result.get(timeout=1)           # prints "100" unless your computer is *very* slow
# print(pool.map(stratified_sample(df=df, strata=[tissue_name], size=sample_size), range(10))) 
# print("--- %s seconds ---" % (time.time() - start_time))


jobs = []
for i in range(5):
    print(i)
    p = multiprocessing.Process(target=stratified_sample(df=df, strata=[tissue_name], size=sample_size), args=(i,))
    jobs.append(p)
    p.start()



# import multiprocessing
# from multiprocessing import Pool
# multiprocessing.set_start_method('fork')

# def f(x):
#     return x*x

# p= Pool(50)

# start_time = time.time()

# print(p.map(f, np.arange(50)))
# print("--- %s seconds ---" % (time.time() - start_time))
# start_time = time.time()
# a = []
# for i in range(50):
#     a.append(f(i))
# print(a)
# print("--- %s seconds ---" % (time.time() - start_time))

# # print(len(output))

start_time1 = time.time()
temp_list = []
for i in range(no_of_iteration):
    temp_stratified = stratified_sample(df=df, strata=[tissue_name], size=sample_size)
    temp_list.append(temp_stratified)
print("--- %s seconds ---" % (time.time() - start_time1))



