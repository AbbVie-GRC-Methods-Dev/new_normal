import os
import pandas as pd
import numpy as np
from pytorch_tabnet.tab_model import TabNetClassifier
from pytorch_tabnet.metrics import Metric
 
import torch
import pickle
from sklearn.metrics import confusion_matrix
from .tabnet_helpers import load_filter_engineered_df, get_internal_features, split_features_and_labels, APS, load_and_merge_new_splits

def performance_metrics(y, y_pred,  split_name):
    """
    y and y_pred are binary vectors derived from a Maf pandas DataFrame
    """
    print(f'split name:  {split_name}, total_variants = {len(y_pred)}')
    cm = confusion_matrix(y, y_pred)
    print(cm)
    TN, FN, TP, FP = cm[0,0], cm[1,0], cm[1,1], cm[0,1]
    total_somatic = TP + FN
    total_germline = TN + FP
    total_variants = total_somatic + total_germline
    print(f'TN: {TN}, FN: {FN}, TP : {TP}, FP : {FP}') 
    print(f'N somatic: {total_somatic},  N germline: {total_germline}, total = {total_variants}')
 
    sens = TP / (TP + FN)
    spec = TN / (TN + FP)
    ppv = TP / (TP + FP)
    npv = TN / (TN + FN)
    f1 = 2*sens*ppv/(sens + ppv) 
    jaccard_index = (TP / (FP + TP + FN))
    print(f'sens: {sens:.3f}, spec: {spec:.3f}, ppv: {ppv:.3f},  npv: {npv:.3f}')
    print(f'F1-score: {f1:.3f}')
    print(f'Jaccard Index: {jaccard_index:.3f}')
    performance_dict = {'sensitivity' : sens, 'specificity' : spec, 'precision' : ppv, \
        'npv': npv, 'f1' : f1, 'jaccard_index' : jaccard_index, 'TP' : TP, 'TN' : TN, \
        'FP' : FP, 'FN' : FN, 'total_somatic' : total_somatic, \
        'total_germline' : total_germline, 'total' : total_variants, 'subset' : split_name}
    return performance_dict


def evaluate_classifier(df, classifier, features,  splits = None, perf_by_sample = False):
    """
    Makes feature matrix and truth labels.  Prints performance via performance_metrics().  Does this for overall performance as well as performance on splits if desired.   
    """
    X_full = df[features]
    X_full = X_full.fillna(0).values
    y_full_truth = df['truth'].values

    y_full_pred = classifier.predict(X_full)
    # initialize list to hold performance dfs.
    classifier_performance_dicts = []

    overall_perf_dict = performance_metrics(y_full_truth, y_full_pred,  'overall')
    overall_perf_dict['is_single_sample'] = False 
    classifier_performance_dicts.append(overall_perf_dict)
    
    if splits is not None:
        for split in splits:
            X_split, y_split_truth = split_features_and_labels(df, split, features)
            y_split_pred = classifier.predict(X_split)
            split_df = performance_metrics(y_split_truth, y_split_pred, split)
            split_df['is_single_sample'] = False 
            classifier_performance_dicts.append(split_df)
    if perf_by_sample:
        for unique_sample in df['sample'].unique():
            sub_df = df[df['sample'] == unique_sample]
            X_sample = sub_df[features]
            X_sample = X_sample.fillna(0).values
            y_sample_truth = sub_df['truth'].values
            y_sample_pred = classifier.predict(X_sample)
            try:
                sample_df = performance_metrics(y_sample_truth, y_sample_pred,  str(unique_sample))
                sample_df['is_single_sample'] = True 
                classifier_performance_dicts.append(sample_df)
            except IndexError:
                print('Confusion matrix is not 2x2.')
                unique_truth_ints = np.unique(y_sample_truth) == 1
                print(unique_truth_ints)
                if len(unique_truth_ints) == 1:
                    if unique_truth_ints[0]:
                        somatic_or_germline = 'somatic'
                    else:
                        somatic_or_germline = 'germline'
                    print(f'All true variants for this sample are {somatic_or_germline}. Skipping.')
                else:
                    print('Double check the binary truth labels for the variants.')
               
    # concat list of DFs to single DF and output. 
    return pd.DataFrame(classifier_performance_dicts)

def make_filename(project_dir, pickled_model, new_suffix):
    """
    new_suffix is the such as '_predictions_plus_features.csv'
    """ 
    file_basename = pickled_model.split('/')[-1].replace('.pkl', new_suffix)
    new_filename = os.path.join(project_dir,'mafs',file_basename) 
    return new_filename

def use_trained_model(project_dir, pickled_model, evaluate_performance, input_features = None, \
                      internal_features_to_drop = None,  input_data = None, splits = None, \
                      add_new_splits = False, perf_by_sample = False):
    """
    Pickled model is a .pkl file produced by training tabnet on data with labels.
    Input data is for if you don't want to use the default project_dir/mafs/engineered.csv file.
    splits is an optional list of splits to evaluate separately, e.g., ['train','validation','test']
    """
    # first make sure features were supplied.
    if (input_features is None) and (internal_features_to_drop is None):
        raise ValueError('Features must be supplied, either through input_features argument \
                        (recommended) or by dropping from list of internal features.')
    # see if path to input data is supplied.  if not, load from relative directory. 
    if input_data is None:
        print('Using default path to engineered input data.')
        df = load_filter_engineered_df(project_dir)
    else:
        df = load_filter_engineered_df(project_dir, input_data)
    # use new splits created in cross-validation fold?
    # if so, use the data_splits in project_dir to define the splits
    print(f'input df has {df.shape[0]} variants.')
    if add_new_splits:
        df = load_and_merge_new_splits(engineered_df = df, fold_dir = project_dir)
                
    classifier = pickle.load(open(pickled_model, 'rb'))

    if input_features is not None: 
        # check if any features are missing from df.
        for input_feature in input_features:
            if input_feature not in df.columns:
                raise ValueError(f'input feature {input_feature} is missing from dataframe.')
        features = input_features 
    else:
        # drop from list of all internal features
        features = get_internal_features( use_purecn_purity = ('purity' in df.columns), internal_features_to_drop = internal_features_to_drop)
    features = list(features)
    print('finished selecting features.')
    print(f'df.columns: {df.columns}')
    print(f'features: {features}')
    print(f'Number of input variants to classify (before eval):  {df.shape[0]}')
    if evaluate_performance:
        performance_df = evaluate_classifier(df, classifier, features, splits, perf_by_sample)
        performance_df.to_csv(make_filename(project_dir, pickled_model, \
            '_full_performance_across_groups.csv'))
        single_samples = performance_df[performance_df['is_single_sample']]
        if single_samples.shape[0] > 0: # if there are some single samples
            print('Performance averaged across single samples.')
            # It wouldn't make sense to average over these two categorical columns
            mean_by_sample = single_samples.drop(columns = ['subset','is_single_sample']).mean()
            print(mean_by_sample)
            mean_by_sample.to_csv(make_filename(project_dir, pickled_model, \
                '_average_performance_stats_by_sample.csv'))
 
    X_full = df[features]
    #print(f'df.shape:  {df.shape}')
    X_full_values = X_full.fillna(0).values
    #print(f'X_full_values.shape: {X_full_values.shape}')

    # make binary predictions and probabilities for saving.
    y_pred = classifier.predict(X_full_values)
    y_proba = classifier.predict_proba(X_full_values)
    #print(f'y_pred.shape: {y_pred.shape}')
    # df already has features plus truth values; add y_pred prediction and class probabilities.
    df['tabnet_pred'] = y_pred
    df['tabnet_proba_0'] = y_proba[:,0]
    df['tabnet_proba_1'] = y_proba[:,1]

    # uncomment if you need to generate the full matrix including unused features.
    output_predictions_path = make_filename(project_dir, pickled_model, '_predictions_all_features.csv')
    print(output_predictions_path)
    df.to_csv(output_predictions_path)
    # drop superfluous unused features from dataframe.
    #if input_features is not None:
         # no good way to do this.  skipping!
    #    df_dropped_unused_features = df[input_features]
    #else:
    #    intersect_features_to_drop = set(df.columns).intersection(set(internal_features_to_drop))
    #    df_dropped_unused_features = df.drop(columns = list(intersect_features_to_drop))
    #df_dropped_unused_features.reset_index(drop = True).to_csv(make_filename(project_dir, pickled_model, '_predictions.csv'))
