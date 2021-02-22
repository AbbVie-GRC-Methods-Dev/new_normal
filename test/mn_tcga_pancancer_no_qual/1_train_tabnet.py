"""
Trains tabnet.  Can specify a name to save the model.  Can also drop features, in this case, dropping the 'qual' feature to see how it impacts the model. 
"""
import context

import os
from new_normal.tabnet_training import train_tabnet 

tcga_project_dir = os.path.dirname(os.path.abspath(__file__))

train_tabnet(tcga_project_dir, save_model_name_as = 'tcga_no_qual', features_to_drop = ('qual',)) 
