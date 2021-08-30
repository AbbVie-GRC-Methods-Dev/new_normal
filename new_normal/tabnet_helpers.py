import os
import numpy as np
import pandas as pd
import glob
import pickle
from pytorch_tabnet.metrics import Metric
from sklearn.metrics import average_precision_score
import matplotlib.pyplot as plt

ALL_INTERNAL_FEATURES = ['qual', 't_depth', 't_alt_freq', 't_maj_allele', 'window_size', 'support', 'prob_somatic', 'mu', 'sigma', 'copy_ratio', 'copy_depth', 'probes', 'count',
           'weight', 'max_cosmic_count', 'inframe_indel', 'missense', 'nonsense', '100mer_mappability', 'ACC', 'ACG', 'ACT', 'ATA', 
           'ATC', 'ATG', 'ATT', 'CCA', 'CCC', 'CCG', 'CCT', 'CTA', 'CTC', 'CTG', 'CTT', 'GCA', 'GCC', 'GCG', 'GCT', 'GTA', 'GTC', 
           'GTG', 'GTT', 'TCA', 'TCC', 'TCG', 'TCT', 'TTA', 'TTC', 'TTG', 'TTT', 'non-SBS_x', 'C>G', 'C>T', 'T>A', 'T>C',
           'T>G', 'non-SBS_y', 'pop_max']

def load_filter_engineered_df(project_dir, input_data = None):           
    if input_data is None:
        engineered_csv_path = os.path.join(project_dir, 'mafs', 'engineered.csv')
    else:
        engineered_csv_path = os.path.join(project_dir, input_data)

    df = pd.read_csv(engineered_csv_path, index_col = 0)                     

    print(f'pre-pop-max filtering, n vars: {df.shape[0]}')
    df = df[df['pop_max']<0.01]                               
    # 
    scpra_keys = ['sample','chrom','pos','ref','alt']
    print(f'pre-deduplication,  n vars: {df.shape[0]}')
    if all(i in df.columns for i in scpra_keys):  
        df.drop_duplicates(subset = scpra_keys, inplace = True )
        print('Deduplicated using sample, chrom, pos, ref, and alt keys')
    else:
        print('sample, chrom, pos, ref, alt missing from df columns.  Skipping deduplication.')
    print(f'post-deduplication, n vars: {df.shape[0]}')
    #print(df.columns)
    df['sample'] = df['sample'].astype('str')
    n_samples = len(np.unique(df['sample']))
    
    print(f'Unique samples: {n_samples}')                                            
    if 'window_size' in df.columns:
        df['window_size'] = df['window_size'].fillna(2.0) # <- could compare this to the max window_size from prob_somatic optimization
    if 'support' in df.columns:
        df['support'] = df['support'].fillna(0.0)                 
    return df


def get_internal_features(use_purecn_purity: bool, internal_features_to_drop = None) -> tuple:
    """
    Features list is common to both train and evaluate TabNet modules.
    use_purecn_purity : bool
    internal_features_to_drop : tuple
    """
    if internal_features_to_drop is not None:
        print('Dropping features.')
        for e in internal_features_to_drop:
            if e not in ALL_INTERNAL_FEATURES:
                raise ValueError(f'Feature {e} not recognized.')
        features = [e for e in ALL_INTERNAL_FEATURES if e not in internal_features_to_drop]
    else:
        features = ALL_INTERNAL_FEATURES
    if use_purecn_purity:
        features.append('purity')
    print(f'Using the following features:  {features}')
    return features

def split_features_and_labels(df, split_name, features):
    split_df = df[df['split'] == split_name]
    n_unique_patients = len(np.unique(split_df['sample']))
    print(f'{n_unique_patients} unique patients in {split_name} data set')
    X = split_df[features]
    X = X.fillna(0).values
    y = split_df['truth'].values
    return X, y

def load_and_merge_new_splits(engineered_df, fold_dir):
    """
    Don't use the splits in the engineered MAF.  
    Use the splits in the data_splits.csv in the directory. 
    This is used for cross-validation across folds, where it  
    avoids creating copies of MAFs that only differ by the splits column.
    """
    if 'split' in engineered_df.columns:
        del engineered_df['split']
    new_data_splits = pd.read_csv(os.path.join(fold_dir,'data_splits.csv'))
    new_data_splits.rename(columns = {'patient' : 'sample'}, inplace = True)
    engineered_df = engineered_df.merge(new_data_splits[['sample','split']], on = 'sample' , how = 'left', validate = 'many_to_one')
    print('successfully updated splits')
    return engineered_df
                    
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

def remove_tmb_outliers(df, min_total_mutations, max_total_mutations):
    unique_samples = df['sample'].unique()
    passing_dfs = []
    n_removed = 0
    for unique_sample in unique_samples:
       sample_df = df[df['sample'] == unique_sample]
       total_mutations = sample_df[sample_df['truth'] == 1].shape[0]
       if (max_total_mutations < total_mutations) | (total_mutations < min_total_mutations):
           print(f'Rejecting sample {unique_sample} with {total_mutations} total mutations.')
           n_removed += 1
       else:
           passing_dfs.append(sample_df) 
    print(f'Removed {n_removed} total patients')
    return pd.concat(passing_dfs)

def add_truth_labels_from_mn(path_engineered_maf_tumor_only, globbed_list_of_mn_maf_paths, overwrite = False, output_file = 'mafs/engineereed_plus_labels.csv'):
    if os.path.exists(output_file) and not overwrite:
        raise FileExistsError('It looks like labels have already been added to the maf.  To overwrite, please set overwrite = True when calling add_truth_labels_from_mn.' )
    else:
        tumor_only_input = pd.read_csv(path_engineered_maf_tumor_only)
        mn_mafs_list = []
        for maf_file in globbed_list_of_mn_maf_paths:
            print(maf_file)
            mn_maf = pd.read_csv(maf_file, sep = '\t', low_memory = False)
            mn_mafs_list.append(mn_maf)
 
        full_mn_maf = pd.concat(mn_mafs_list)
        left_joined = pd.merge(tumor_only_input, full_mn_maf[['sample', 'chrom','pos','ref', 'alt', 'filter']], on = ['sample','chrom','pos','ref','alt'], how = 'left', suffixes = ['_to','_mn'] )
        left_joined['truth'] = left_joined['filter_mn'].apply(lambda x: 1 if x == 'PASS' else 0)
        left_joined.to_csv(output_file)

def print_global_feature_importance(clf, features):
    """
    An opened classifier object and a list of features used to train the classifier.
    """
    f_i = clf.feature_importances_  
    sorted_indices = np.argsort(f_i)
    print(np.array(features)[indices])
    

def plot_global_feature_importance(clf, features, save_plot_as = None):
    f_i = clf.feature_importances_  
    sorted_indices = np.argsort(f_i)
    n_features = len(features)

    # Set up plot. 
    plt.figure(figsize = (11,12))
    plt.title("TabNet global feature importances")
    plt.barh(range(n_features), f_i[sorted_indices],
           color="r", align="center")
     
    # Name the ticks.
    plt.yticks(range(n_features), np.array(features)[sorted_indices])
    plt.ylim([-1, n_features])
    plt.show()
    if save_plot_as is None:
        plt.savefig('tabnet_global_feature_importances.png')      
    elif save_plot_as[:-4:] == '.png':
        plt.savefig(save_plot_as)
    else:
        plt.savefig(save_plot_as + '.png')
    plt.close()


def save_global_feature_importance(fold_dir, features):
    figures_dir = os.path.join('global_feature_importance', fold_dir)
    if not os.path.exists(figures_dir):
        os.makedirs(figures_dir)
                
    fi_all_models = []
    for pkl in glob.glob(fold_dir + '*/*.pkl'):
        fold_name = pkl.split('/')[-2]
        clf = pickle.load(open(pkl, 'rb'))
        fi_all_models.append(clf.feature_importances_)
        plot_name = os.path.join(figures_dir, fold_name + '_global_feature_importance')
        plot_global_feature_importance(clf = clf, features = features, save_plot_as = plot_name)
                
    fi_dfs = [] 
    for i, model_fi in enumerate(fi_all_models):
        fi_dfs.append(pd.DataFrame({'feature' : features, 'importance' : model_fi, 'model_number' :[i] * len(features) }))
                
    full_feature_importance = pd.concat(fi_dfs)
    full_feature_importance.to_csv(os.path.join(figures_dir, 'full_feature_importance.csv'))

import importlib                                                                    
 
def modify_and_import(module_name, package, modification_func):
    spec = importlib.util.find_spec(module_name, package)
    source = spec.loader.get_source(module_name)
    new_source = modification_func(source)
    module = importlib.util.module_from_spec(spec)
    codeobj = compile(new_source, module.__spec__.origin, 'exec')
    exec(codeobj, module.__dict__)
    sys.modules[module_name] = module
    return module
