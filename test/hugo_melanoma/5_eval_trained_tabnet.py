"""
Evaluates the TCGA-trained TabNet model on the Hugo 2016 Melanoma dataset.
"""
import context

import os
from new_normal.model_predict_and_eval import use_trained_model 

hugo_project_dir = os.path.dirname(os.path.abspath(__file__))
use_trained_model(hugo_project_dir, pickled_model = '../tcga_pancancer/tabnet_trained_tcga.pkl', evaluate_performance = True, input_data = 'mafs/engineered_plus_labels.csv', splits = None, perf_by_sample = True) 
