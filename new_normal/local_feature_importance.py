import os
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import numpy as np
from .tabnet_helpers import get_internal_features

def plot_explain(df_subset, features_to_drop, pickled_model, fig_name = 'explain_plot.png'):
    classifier = pickle.load(open(pickled_model, 'rb'))
    features = get_internal_features( use_purecn_purity = ('purity' in df_subset.columns), internal_features_to_drop = features_to_drop)
    X = df_subset[features]
    X_values = X.fillna(0).values
    explain_matrix, masks = classifier.explain(X_values)

    fig, axs = plt.subplots(4, 1, figsize=(20, 20))
    for i in range(4):
        axs[i].imshow(masks[i][:50].T)
        axs[i].set_title(f'mask {i}')
        axs[i].set_yticks(np.arange(len(X.columns)))
        axs[i].set_yticklabels(X.columns)
    fig_dir_out = 'local_explainability_plots'
    if not os.path.exists(fig_dir_out):
        os.mkdir(fig_dir_out)
    plt.title(fig_name) 
    plt.savefig(os.path.join(fig_dir_out, fig_name))
