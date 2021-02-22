"""
Evaluates tabnet across splits.
"""
import context

import os
from new_normal.model_predict_and_eval import use_trained_model 

tcga_project_dir = os.path.dirname(os.path.abspath(__file__))

use_trained_model(tcga_project_dir, pickled_model = 'tabnet_trained_tcga_no_qual.pkl', evaluate_performance = True, features_to_drop = ('qual',),splits = ['train', 'validation','test'], perf_by_sample = True) 
