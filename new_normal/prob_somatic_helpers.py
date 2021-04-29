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

PROB_SOMATIC_HYPERPARAMS = {'info_snp_pop' : 0.001, 'info_snp_freq' : 0.95,
                                'info_snp_depth' : 10.0, 'min_info_snp' : 20.0,
                                'expand_copy_window_by' : 0.01, 'max_copy_window_size' : 1.0, 'vaf_bin_size' : 0.01}

def calc_prob_somatic(row, hyperparams=None):
    
    snv_path =  row.snv
    cnv_path = row.cnv
    prob_somatic = ProbSomatic(merge_somatic((snv_path, cnv_path)), **hyperparams)
    prob_somatic._format_data()
    prob_somatic._calc_prob()
    
    return prob_somatic.prob_somatic


def multiprocess_samples(data_splits, params = None, max_cores = None):
    '''
    Ships jobs to process_samples, one sample/core
    '''
    
        
    dfs = []
    print(f'Running ProbSomatic using {max_cores} cores to process samples in parallel.')
    with Pool(max_cores) as p:
        dfs.extend(p.map(partial(calc_pro_somatic, hyperparams = hyperparams), data_splits.iterrows()))

    return pd.concat(dfs)

def get_prob_somatic_mafs(project_dir, max_cores = None, output_dir = None, hyperparams = None): 
    data_splits = pd.read_csv(os.path.join(project_dir,'data_splits.csv'))
    
    # these default ProbSomatic hyperparameters are reasonably close to the optimal parameters found using random search.
    # Using these hyperparams, the resulting output works quite well as a feature in the downstream machine learning model. 
    if hyperparams is None:
        prob_somatic_hyperparams = PROB_SOMATIC_HYPERPARAMS
    else:
        prob_somatic_hyperparams = hyperparams
    print(prob_somatic_hyperparams)
    prob_somatic_annotated = multiprocess_samples(data_splits, params=prob_somatic_hyperparams, max_cores = max_cores)
    if output_dir is None: 
        prob_somatic_output_dir = os.path.join(project_dir, 'mafs/prob_somatic/') 
    else:
        prob_somatic_output_dir = output_dir
    if not os.path.exists(prob_somatic_output_dir):
        os.mkdir(prob_somatic_output_dir)
     
    # save individual mafs to prob_somatic dir
    for sample in prob_somatic_annotated['sample'].unique():
        #print(f'prob somatic data frame for sample {sample}')
        sample_df =  prob_somatic_annotated[prob_somatic_annotated['sample'] == sample]
        sample_df_out = os.path.join(prob_somatic_output_dir, str(sample) + '.maf')
        print(f'saving {sample_df_out}')
        sample_df.to_csv(sample_df_out, sep='\t', index_label=False)
