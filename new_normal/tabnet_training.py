import os
import numpy as np
from pytorch_tabnet.tab_model import TabNetClassifier

from .tabnet_helpers import load_filter_engineered_df, get_features, split_features_and_labels, APS

import torch
import pickle

def train_tabnet(project_dir, save_model_name_as, shuffle_labels = False):
    """
    Shuffle the  labels?  This is a negative control mode used to make sure there's no information leak between train and test.
    """

    df = load_filter_engineered_df(project_dir)
                      
    features = get_features(use_purecn_purity = ('purity' in df.columns))  
    
    # shuffle variant rows.  split status will remain paired with each variant.
    df = df.sample(frac = 1)
    
    if shuffle_labels:
        # shuffle the truth column.  Should still be able to fit train, but test should be terrible.
        df.truth = np.random.permutation(df.truth.values)
        # edit:  while tabnet can fit the tcga training data nearly perfectly (~99.9% APS),
        # it cannot get above 37% APS on the same data with shuffled labels.  
        # this is suggests tabnet isn't over-parametrized.
    
    X_train, y_train = split_features_and_labels(df, 'train', features)
    X_validation, y_validation = split_features_and_labels(df, 'validation', features)
    
    print('finished building train / validation splits.')
    
    clf = TabNetClassifier(
        n_d=64, n_a=64, n_steps=5,
        gamma=1.5, n_independent=2, n_shared=2,
        lambda_sparse=1e-4, momentum=0.3, clip_value=2.,
        optimizer_fn=torch.optim.Adam,
        optimizer_params=dict(lr=2e-2),
        scheduler_params = {"gamma": 0.95,
                         "step_size": 20},
        scheduler_fn=torch.optim.lr_scheduler.StepLR, epsilon=1e-15
    )
    
    max_epochs = 500
    print('fitting the model')
    clf.fit(
        X_train=X_train, y_train=y_train,
        eval_set=[(X_train, y_train), (X_validation, y_validation)],
        eval_name=['train', 'valid'],
        eval_metric = [APS],
        max_epochs=max_epochs, patience=100,
        batch_size=4000, virtual_batch_size=256
    )
    
    pickle_out_name = os.path.join(project_dir, 'tabnet_trained_' + save_model_name_as + '.pkl') 
    if shuffle_labels:
        pickle_out_name = pickle_out_name.replace('.pkl', 'shuffled_labels.pkl') 
    pickle.dump(clf, open(pickle_out_name,'wb'))
