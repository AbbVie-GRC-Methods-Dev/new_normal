import os
import numpy as np
import pandas as pd
from pytorch_tabnet.tab_model import TabNetClassifier

from .tabnet_helpers import load_filter_engineered_df, get_internal_features, split_features_and_labels, APS, load_and_merge_new_splits

import torch
import pickle

def train_tabnet(project_dir, save_model_name_as, lambda_sparse = 1e-4,n_steps = 4, input_features = None, internal_features_to_drop = None, shuffle_labels = False, input_data = None, add_new_splits = False):
    """
    Shuffle the  labels?  This is a negative control mode used to make sure there's no information leak between train and test.
    """
    if (input_features is None) and (internal_features_to_drop is None):
        raise ValueError('Features must be supplied.')

    if input_data is None:
        df = load_filter_engineered_df(project_dir)
    else:
        df = load_filter_engineered_df(project_dir, input_data)
    if add_new_splits:
        df = load_and_merge_new_splits(engineered_df = df, fold_dir = project_dir)
                      
    # use explicitly defined features if given.
    if input_features is not None:
        features = input_features
    else:  # otherwise, use the list of features to drop from the set of internal features.
        features = get_internal_features(use_purecn_purity = ('purity' in df.columns), internal_features_to_drop = internal_features_to_drop)  
    print('The following features are being used for training:' )
    print(features)
    
    # shuffle variant rows.  split status will remain paired with each variant.
    df = df.sample(frac = 1)
    
    if shuffle_labels:
        # shuffle the truth column.  Should still be able to fit train, but test should be terrible.
        df.truth = np.random.permutation(df.truth.values)
    
    X_train, y_train = split_features_and_labels(df, 'train', features)
    X_validation, y_validation = split_features_and_labels(df, 'validation', features)
    print(f'x train shape: {X_train.shape}')
    print(f'y train shape: {y_train.shape}')
    print(f'x val shape: {X_validation.shape}')
    print(f'y val shape: {y_validation.shape}')
    print(y_validation)
    
    print('finished building train / validation splits.')
    
    clf = TabNetClassifier(
        #n_d=64, n_a=64, n_steps=5,
        #n_d=24, n_a=24, n_steps=4,
        n_d=24, n_a=24, n_steps=n_steps,
        #n_d=32, n_a=32, n_steps=10,
        device_name = 'auto',
        #device_name = 'cpu',
        gamma=1.5, n_independent=2, n_shared=2,
        #lambda_sparse=1e-4, momentum=0.3, clip_value=2.,
        # after first decent model (late feb 2021), trying increasing lambda sparse for regularization
        #lambda_sparse=1e-2, momentum=0.3, clip_value=2.,
        # new functionalized version
        lambda_sparse=lambda_sparse, momentum=0.3, clip_value=2.,
        optimizer_fn=torch.optim.Adam,
        optimizer_params=dict(lr=2e-2),
        scheduler_params = {"gamma": 0.95,
                         "step_size": 20},
        scheduler_fn=torch.optim.lr_scheduler.StepLR, epsilon=1e-15
    )
    
    max_epochs = 100
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
    # also want to save the network weights.  
    # the pickled classifier is great if you're doing inference using the 
    # same device the network was trained on. 
    #  However, the classifier cannot be transferred across devices like 
    # GPU-trained to CPU-evaluated.  
    state_dict_out_name = pickle_out_name.replace('.pkl','.pth')
    torch.save(clf.network.state_dict(), state_dict_out_name)
