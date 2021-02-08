import os
import numpy as np
import pandas as pd
from pytorch_tabnet.metrics import Metric
from sklearn.metrics import average_precision_score

def load_filter_engineered_df(project_dir, input_data = None):           
    if input_data is None:
        engineered_csv_path = os.path.join(project_dir, 'mafs', 'engineered.csv')
    else:
        engineered_csv_path = os.path.join(project_dir, input_data)

    df = pd.read_csv(engineered_csv_path)                     

    df = df[df['pop_max']<0.01]                               
    n_samples = len(np.unique(df['sample']))

    print(f'Unique samples: {n_samples}')                                            
    df['window_size'] = df['window_size'].fillna(2.0) # <- could compare this to the max window_size from prob_somatic optimization
    df['support'] = df['support'].fillna(0.0)                 
    return df

def get_features(use_purecn_purity):
    """
    Features list is common to both train and evaluate TabNet modules.
    """
    features = ['qual', 't_depth', 't_maj_allele', 'window_size', 'support', 'prob_somatic', 'mu', 'sigma', 'copy_ratio', 'copy_depth', 'probes',
           'weight', 'max_cosmic_count', 'inframe_indel', 'missense', 'nonsense', '100mer_mappability', 'ACC', 'ACG', 'ACT', 'ATA', 
           'ATC', 'ATG', 'ATT', 'CCA', 'CCC', 'CCG', 'CCT', 'CTA', 'CTC', 'CTG', 'CTT', 'GCA', 'GCC', 'GCG', 'GCT', 'GTA', 'GTC', 
           'GTG', 'GTT', 'TCA', 'TCC', 'TCG', 'TCT', 'TTA', 'TTC', 'TTG', 'TTT', 'non-SBS_x', 'C>G', 'C>T', 'T>A', 'T>C',
           'T>G', 'non-SBS_y', 'pop_max']
    if use_purecn_purity:
        features.append('purity')
    return features

def split_features_and_labels(df, split_name, features):
    split_df = df[df['split'] == split_name]
    n_unique_patients = len(np.unique(split_df['sample']))
    print(f'{n_unique_patients} unique patients in {split_name} data set')
    X = split_df[features]
    X = X.fillna(0).values
    y = split_df['truth'].values
    return X, y
                    
class APS(Metric):          
    """                     
    Average Precision Score.  A more reliable 'area under the precision recall curve'.
    # This custom evaluation metric APS is used for training TabNet.
    # Sklearn's f1-score does not work because it requires a 0 or 1 for the prediction,
    # and cannot handle a probability output by the neural network.
    """                     
    def __init__(self):     
        self._name = "aps"  
        self._maximize = True
                            
    def __call__(self, y_true, y_score):
        aps = average_precision_score(y_true, y_score[:, 1])
        return aps          

