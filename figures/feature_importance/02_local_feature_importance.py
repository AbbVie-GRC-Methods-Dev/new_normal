import os
import sys
import pandas as pd
import numpy as np
import pickle
import pytorch_tabnet
from scipy.stats import ranksums

sys.path.append('/bioinformatics/Users/mclaurt/projects/deep_learning_tumor-only_variant_calling/')
sys.path.append('/bioinformatics/Users/mclaurt/projects/deep_learning_tumor-only_variant_calling/test/no_ps_fit_with_snp-vaf-bin/tcga_pancancer')
from train_features import vaf_bin_mode_features
import matplotlib.pyplot as plt



pickled_model = '/bioinformatics/Users/mclaurt/projects/deep_learning_tumor-only_variant_calling/test/no_ps_fit_with_snp-vaf-bin/tcga_pancancer/tabnet_trained_tcga_v3.pkl'
classifier = pickle.load(open(pickled_model, 'rb'))

results_path = '/bioinformatics/Users/mclaurt/projects/deep_learning_tumor-only_variant_calling/test/no_ps_fit_with_snp-vaf-bin/tcga_pancancer/mafs/tabnet_trained_tcga_v3_predictions_all_features.csv'

full_results = pd.read_csv(results_path)
# filter results here, i.e., to include only FPs.
# test results
test_results = full_results[full_results.split == 'test']
test_fp = test_results[(test_results.truth == 0) & (test_results.tabnet_pred == 1)]
test_tn = test_results[(test_results.truth == 0) & (test_results.tabnet_pred == 0)]
print(f'n rows test FP:  {test_fp.shape[0]}')
print(f'n rows test TN:  {test_tn.shape[0]}')

fp_cc_one_plus = test_fp[test_fp.max_cosmic_count >= 1]
tn_cc_one_plus = test_tn[test_tn.max_cosmic_count >= 1]

print(f'FP rows with nonzero CC :  {fp_cc_one_plus.shape[0]}')
print(f'TN rows with nonzero CC :  {tn_cc_one_plus.shape[0]}')

print(f'median FP CC:  {test_fp["max_cosmic_count"].median()}')
print(f'median TN CC:  {test_tn["max_cosmic_count"].median()}')

rs_fp_tn_cc = ranksums(test_fp['max_cosmic_count'].values, test_tn['max_cosmic_count'].values)
print(f'ranksums FP vs TN CC:  {rs_fp_tn_cc}')

print(f'median FP count:  {test_fp["count"].median()}')
print(f'median TN count:  {test_tn["count"].median()}')

rs_fp_tn_count = ranksums(test_fp['count'].values, test_tn['count'].values)
print(f'ranksums FP vs TN count:  {rs_fp_tn_count}')

test_fn = test_results[(test_results.truth == 1) & (test_results.tabnet_pred == 0)]
test_tp = test_results[(test_results.truth == 1) & (test_results.tabnet_pred == 1)]



print(f'median FP t alt freq: {test_fp.t_alt_freq.median()}')
print(f'median TN t alt freq: {test_tn.t_alt_freq.median()}')



#X = test_fp[vaf_bin_mode_features]
#X_values = X.fillna(0).values
#explain_matrix, masks = classifier.explain(X_values)
#print('explain_matrix:')
#print(explain_matrix)
#print('masks:')
#print(masks)
#print(X.max_cosmic_count[:50])

def four_plots(X, masks, fig_name, plot_n_features = None):
    np.random.seed(42)
    n_rows = masks[0].shape[0]
    random_indices = np.random.choice(n_rows, size=50, replace=False)
    fig, axs = plt.subplots(4, 1, figsize=(10, 10))
    for i in range(4):        
        axs[i].imshow(masks[i][random_indices].T[:plot_n_features])
        axs[i].set_title(f'mask {i}')
        axs[i].set_yticks(np.arange(len(X.columns[:plot_n_features])))
        axs[i].set_yticklabels(X.columns[:plot_n_features])
    fig_dir_out = 'local_explainability_plots'
    if not os.path.exists(fig_dir_out):
        os.mkdir(fig_dir_out) 
    plt.subplots_adjust(left = 0.3, bottom = 0.01, top = 0.50, hspace=0.0)
    plt.title(fig_name)       
    plt.savefig(os.path.join(fig_dir_out, fig_name))

def get_average_local_importance(df_subset, features):
    X = df_subset[features]#test_fp[vaf_bin_mode_features]
    X_values = X.fillna(0).values
    explain_matrix, masks = classifier.explain(X_values)
    print(masks[1].shape)
    print(df_subset.shape)
    masks_averaged = {k : np.mean(v, axis = 0) for k,v in masks.items()}
    local_importance = pd.DataFrame({'feature' : X.columns, 'feature_average' : X.mean(axis = 0).values, 'average_weight_0' : masks_averaged[0],\
        'average_weight_1' : masks_averaged[1], 'average_weight_2' : masks_averaged[2], 'average_weight_3' : masks_averaged[3] })
    local_importance = local_importance.sort_values(by = 'average_weight_0', ascending = False)
    print(local_importance.iloc[:60])
    return local_importance

#def get_patient_average_local_importance(df_subset, features):
#    X = df_subset[features]#test_fp[vaf_bin_mode_features]
#    X_values = X.fillna(0).values
#    explain_matrix, masks = classifier.explain(X_values)
#    print(masks[1].shape)
#    print(df_subset.shape)
#    masks_patient_averaged = {k : np.mean(v, axis = 0) for k,v in masks.items()}
#    masks_averaged = {k : np.mean(v, axis = 0) for k,v in masks.items()}
#    local_importance = pd.DataFrame({'feature' : X.columns, 'feature_average' : X.mean(axis = 0).values, 'average_weight_0' : masks_averaged[0],\
#        'average_weight_1' : masks_averaged[1], 'average_weight_2' : masks_averaged[2], 'average_weight_3' : masks_averaged[3] })
#    local_importance = local_importance.sort_values(by = 'average_weight_0', ascending = False)
#    print(local_importance.iloc[:60])
#    return local_importance

print('fp')
get_average_local_importance(test_fp, vaf_bin_mode_features)
print('tn')
get_average_local_importance(test_tn, vaf_bin_mode_features)
print('fn')
get_average_local_importance(test_fn, vaf_bin_mode_features)
print('tp')
get_average_local_importance(test_tp, vaf_bin_mode_features)

#X = test_fp[vaf_bin_mode_features]
#X_values = X.fillna(0).values
#explain_matrix, masks = classifier.explain(X_values)
X = test_fp[vaf_bin_mode_features]
X_values = X.fillna(0).values
explain_matrix, masks = classifier.explain(X_values)
four_plots(X, masks, 'tcga_test_fp.png', plot_n_features = 5)

X = test_tn[vaf_bin_mode_features]
X_values = X.fillna(0).values
explain_matrix, masks = classifier.explain(X_values)
four_plots(X, masks, 'tcga_test_tn.png', plot_n_features = 5)

X = test_tp[vaf_bin_mode_features]
X_values = X.fillna(0).values
explain_matrix, masks = classifier.explain(X_values)
four_plots(X, masks, 'tcga_test_tp.png', plot_n_features = 5)

X = test_fn[vaf_bin_mode_features]
X_values = X.fillna(0).values
explain_matrix, masks = classifier.explain(X_values)
four_plots(X, masks, 'tcga_test_fn.png', plot_n_features = 5)
