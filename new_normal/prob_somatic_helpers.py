#!/usr/bin/env python
# coding: utf-8

import math
import os
from itertools import chain
from functools import partial
from multiprocessing import Pool

import pandas as pd
import numpy as np

from .prob_somatic import ProbSomatic

pd.set_option('display.max_rows', 200)

def process_sample(prob_somatic, sample):
    '''
    Called by multiprocess_samples to do one sample/core. 
    '''
    prob_somatic.add_samples([sample])
    
    return prob_somatic.prob_somatic

def multiprocess_samples(samples, params = None, max_cores = None):
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
    with Pool(max_cores) as p:
        dfs.extend(p.map(partial(process_sample, prob_somatic), samples.index.tolist()))

    return pd.concat(dfs)

def get_prob_somatic_mafs(project_dir, max_cores = None): 
    data = pd.read_csv(os.path.join(project_dir,'data_splits.csv'))
    
    # these ProbSomatic hyperparameters are reasonably close to the optimal parameters found using random search.
    # Using these hyperparams, the resulting output works quite well as a feature in the downstream machine learning model. 
    prob_somatic_hyperparams = {'info_snp_pop' : 0.001, 'info_snp_freq' : 0.95,
                                'info_snp_depth' : 10.0, 'min_info_snp' : 20.0,
                                'expand_copy_window_by' : 0.01, 'max_copy_window_size' : 1.0, 'vaf_bin_size' : 0.01}
    
    prob_somatic_annotated = multiprocess_samples(data, params=prob_somatic_hyperparams, max_cores = max_cores)
    prob_somatic_output_dir = os.path.join(project_dir, 'mafs/prob_somatic/') 
    if not os.path.exists(prob_somatic_output_dir):
        os.mkdir(prob_somatic_output_dir)
     
    # save individual mafs to prob_somatic dir
    for sample in prob_somatic_annotated['sample'].unique():
        prob_somatic_annotated[prob_somatic_annotated['sample'] == sample].to_csv(os.path.join(prob_somatic_output_dir, sample + '.maf'), sep='\t', index_label=False)
