import os
import numpy as np
from pytorch_tabnet.tab_model import TabNetClassifier
from pytorch_tabnet.metrics import Metric
 
import torch
import pickle
from sklearn.metrics import confusion_matrix
from .tabnet_helpers import load_filter_engineered_df, get_features, split_features_and_labels, APS

def performance_metrics(y, y_pred,  split_name):
    """
    y and y_pred are binary vectors derived from a Maf pandas DataFrame
    """
    print(f'split name:  {split_name}')
    cm = confusion_matrix(y, y_pred)
    TN, FN, TP, FP = cm[0,0], cm[1,0], cm[1,1], cm[0,1]
    print(f'TN: {TN}, FN: {FN}, TP : {TP}, FP : {FP}')
 
    sens = TP / (TP + FN)
    spec = TN / (TN + FP)
    ppv = TP / (TP + FP)
    npv = TN / (TN + FN)
 
    print(f'sens: {sens:.3f}, spec: {spec:.3f}, ppv: {ppv:.3f},  npv: {npv:.3f}')
    print(f'Jaccard Index: {(TP / (FP + TP + FN)):.3f}')


def evaluate_classifier(df, classifier, features,  splits = None, perf_by_sample = False):
    """
    Makes feature matrix and truth labels.  Prints performance via performance_metrics().  Does this for overall performance as well as performance on splits if desired.   
    """
    X_full = df[features]
    X_full = X_full.fillna(0).values
    y_full_truth = df['truth'].values

    y_full_pred = classifier.predict(X_full)
    performance_metrics(y_full_truth, y_full_pred,  'overall')

    if splits is not None:
        for split in splits:
            X_split, y_split_truth = split_features_and_labels(df, split, features)
            y_split_pred = classifier.predict(X_split)
            performance_metrics(y_split_truth, y_split_pred, split)
    if perf_by_sample:
        for unique_sample in df['sample'].unique():
            sub_df = df[df['sample'] == unique_sample]
            X_sample = sub_df[features]
            X_sample = X_sample.fillna(0).values
            y_sample_truth = sub_df['truth'].values
            y_sample_pred = classifier.predict(X_sample)
            performance_metrics(y_sample_truth, y_sample_pred,  'overall')

def use_trained_model(project_dir, pickled_model, evaluate_performance, input_data = None, splits = None, perf_by_sample = False):
    """
    Pickled model is a .pkl file produced by training tabnet on data with labels.
    Input data is for if you don't want to use the default project_dir/mafs/engineered.csv file.
    splits is an optional list of splits to evaluate separately, e.g., ['train','validation','test']
    """
    if input_data is None:
        df = load_filter_engineered_df(project_dir)
    else:
        df = load_filter_engineered_df(project_dir, input_data)
    classifier = pickle.load(open(pickled_model, 'rb'))

    features = get_features( use_purecn_purity = ('purity' in df.columns) )

    if evaluate_performance:
        evaluate_classifier(df, classifier, features, splits, perf_by_sample)
   
    X_full = df[features]
    X_full = X_full.fillna(0).values

    # make binary predictions and probabilities for saving.
    y_pred = classifier.predict(X_full)
    y_proba = classifier.predict_proba(X_full)
    
    df['tabnet_pred'] = y_pred
    df['tabnet_proba_0'] = y_proba[:,0]
    df['tabnet_proba_1'] = y_proba[:,1]
    
    predictions = os.path.join(project_dir,'mafs', pickled_model.split('/')[-1].replace('.pkl', '_predictions.csv'))
    df.to_csv(predictions)
