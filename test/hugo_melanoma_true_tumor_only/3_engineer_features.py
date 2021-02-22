"""
This file imports the function that creates the engineered.csv MAF of all patients, one-hot-encoded. 
again, data_splits.csv is used as input for the cohort.
"""

import context

import os
from new_normal.feature_engineering import engineer_features 

hugo_project_dir = os.path.dirname(os.path.abspath(__file__))

# generates engineered.csv
engineer_features(hugo_project_dir, use_purecn_purity = False)
