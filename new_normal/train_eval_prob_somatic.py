#!/usr/bin/env python
# coding: utf-8

import random, math
import os
from itertools import chain
from functools import partial
from multiprocessing import Pool

import pandas as pd
import numpy as np
import scipy.stats as ss
import plotly.express as px
from IPython.display import display_html 

from prob_somatic import ProbSomatic, performance_plots, make_fig, normal_curve

pd.set_option('display.max_rows', 200)

def param_table(prob_somatic):
    
    info_snp_pop = prob_somatic.info_snp_pop 
    info_snp_freq = prob_somatic.info_snp_freq 
    info_snp_depth = prob_somatic.info_snp_depth 
    min_info_snp = prob_somatic.min_info_snp 
    expand_copy_window_by = prob_somatic.expand_copy_window_by 
    max_copy_window_size = prob_somatic.max_copy_window_size
    vaf_bin_size = prob_somatic.vaf_bin_size
    
    params = [info_snp_pop, info_snp_freq, info_snp_depth, min_info_snp, expand_copy_window_by, max_copy_window_size, vaf_bin_size]
    
    return pd.DataFrame([['info_snp_pop', 'info_snp_freq', 'info_snp_depth', 'min_info_snp', 'expand_copy_window_by',
                          'max_copy_window_size', 'vaf_bin_size'], 
                          params]).T.rename(columns={0 : 'parameter', 1 : 'value'}).set_index('parameter').round(3)


def engine(trial, samples = None):
    """
    This does the calculation, and assesses performance. 
    """
    
    # Numpy will seed the same every iteration without this
    np.random.seed(int.from_bytes(os.urandom(4), byteorder = 'little'))
    
    # Sample during training
    info_snp_pop = dist['info_snp_pop'].rvs() 
    info_snp_freq = dist['info_snp_freq'].rvs() 
    info_snp_depth = dist['info_snp_depth'].rvs() 
    min_info_snp = dist['min_info_snp'].rvs() 
    expand_copy_window_by = dist['expand_copy_window_by'].rvs() 
    max_copy_window_size = dist['max_copy_window_size'].rvs()
    vaf_bin_size = dist['vaf_bin_size'].rvs()
    
    # Instantiate the class with the parameters
    prob_somatic = ProbSomatic(train,
                                   info_snp_pop = info_snp_pop, 
                                   info_snp_freq = info_snp_freq, 
                                   info_snp_depth = info_snp_depth, 
                                   min_info_snp = min_info_snp, 
                                   expand_copy_window_by = expand_copy_window_by, 
                                   max_copy_window_size = max_copy_window_size,
                                   vaf_bin_size = vaf_bin_size)
    # Estimate the probabilities!
    prob_somatic.add_samples(samples)
    
    # Calulate a score that is defined using the AUC and Average Precision Score
    full = prob_somatic.prob_somatic
    df = full[full['filter'].str.contains('PASS')]
    _, auc, aps = performance_plots(df)
    score = int(round((auc**4 + aps**4)*100)) # Heavily reward good score (or penalize bad ones)
    
    # !!! Cannot remember exactly what this is doing as far as results goes: will have to come back and comment !!!
    parameters = param_table(prob_somatic)
    parameters = dict(zip(parameters.index, parameters['value']))
    results = {}
    for parameter, value in parameters.items():
        results[parameter] = value
       
    print(' -- Trial {} complete: auc {}; aps {}; score {}'.format(trial, auc, aps, score))

    return score, results


# In[5]:


def process_sample(prob_somatic, sample):
    '''
    Called by multiprocess_samples to do one sample/core. 
    '''
    
    prob_somatic.add_samples([sample])
    
    return prob_somatic.prob_somatic


# In[6]:


def multiprocess_samples(samples, params=None):
    '''
    Ships jobs to process_samples, one sample/core
    '''
    
    
    prob_somatic = ProbSomatic(samples,
                               info_snp_pop = params['info_snp_pop'], 
                               info_snp_freq = params['info_snp_freq'], 
                               info_snp_depth = params['info_snp_depth'], 
                               min_info_snp = params['min_info_snp'], 
                               expand_copy_window_by = params['expand_copy_window_by'], 
                               max_copy_window_size = params['max_copy_window_size'],
                               vaf_bin_size = params['vaf_bin_size'])
        
    dfs = []
    with Pool() as p:
        dfs.extend(p.map(partial(process_sample, prob_somatic), samples.index.tolist()))

    return pd.concat(dfs)


# In[7]:


# Initialize the parameter ranges...
param_ranges = {
    'info_snp_pop':(0.001, 0.05),
    'info_snp_freq':(0.90, 0.980),
    'info_snp_depth':(5, 50),
    'min_info_snp':(5, 25),
    'expand_copy_window_by':(0.001, 0.05),
    'max_copy_window_size':(0.5, 2.5),
    'vaf_bin_size':(0.01, 0.1)
}
# Make a dictionary of uniform distributions--defined by the above parameter ranges--for each parameter. 
# This is because the 'engine' samples to get parameters
dist = dict([(parameter, ss.uniform(rang[0], rang[1] - rang[0])) for parameter, rang in param_ranges.items()])


# In[8]:


data = pd.read_csv('data_splits.csv')


# In[9]:


train = data[data['split'] == 'train']
test = data[data['split'] == 'test']
validation = data[data['split'] == 'validation']


# In[40]:


# Run many trials with different parameter combinations
trials = 500
params_train = []
samples=random.sample(train.index.tolist(), train.shape[0]) # <- should do this on all training samples so I don't get (un)lucky!! (i.e., train.shape[0])
with Pool() as p: # Here were parellelize!
    params_train.extend(p.map(partial(engine, samples = samples), list(range(1,trials + 1))))


# In[43]:


results = dict([(index,score[0]) for index,score in enumerate(params_train)])
min_index = min(results, key=results.get)
min_score = results[min_index]
min_params = params_train[min_index][1]
max_index = max(results, key=results.get)
max_score = results[max_index]
max_params = params_train[max_index][1]
print('max score is {} at index {} and min score is {} at index {}'.format(max_score, max_index, min_score, min_index))


# In[44]:


max_params


# In[50]:


my_params = {'info_snp_pop' : 0.001,
             'info_snp_freq' : 0.95,
             'info_snp_depth' : 10.0,
             'min_info_snp' : 20.0,
             'expand_copy_window_by' : 0.01,
             'max_copy_window_size' : 1.0,
             'vaf_bin_size' : 0.01}


# In[25]:


# Before finding missing mutations, these were the params we used. 
#my_params = {'info_snp_pop': 0.001,
#             'info_snp_freq': 0.95,
#             'info_snp_depth': 10,
#             'min_info_snp': 25,
#             'expand_copy_window_by': 0.01,
#             'max_copy_window_size': 1.5,
#             'vaf_bin_size': 0.01}


# In[51]:


annotate_train = multiprocess_samples(train, params=my_params)


# In[52]:


df = annotate_train[annotate_train['filter'].str.contains('PASS')]
fig, auc, aps = performance_plots(df)
fig.show()


# In[53]:


annotate_test = multiprocess_samples(test, params=my_params)


# In[54]:


df = annotate_test[annotate_test['filter'].str.contains('PASS')]
fig, auc, aps = performance_plots(df)
fig.show()


# In[55]:


annotate_all = multiprocess_samples(data, params=my_params)


# In[56]:


for sample in annotate_all['sample'].unique():
    annotate_all[annotate_all['sample']==sample].to_csv('/bioinformatics/Users/masicdx/tcga/maf/prob_somatic/' + sample + '.maf', sep='\t', index_label=False)


# In[ ]:




