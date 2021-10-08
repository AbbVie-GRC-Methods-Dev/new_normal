import os
import time
import pandas as pd
from .gather_input_data import wes_pipeline_to_ml_input
from .default_features import prob_somatic_hyperparams 
from .prob_somatic_helpers import get_prob_somatic_mafs
from .feature_engineering import engineer_features 
from .model_predict_and_eval import use_trained_model 

def benchmark(trained_model, input_features, truth_labels = None):
    if not os.path.exists(trained_model):
        raise FileNotFoundError
    perf_outdir = 'single_vcfs'   
    if not os.path.exists(perf_outdir):
        os.mkdir(perf_outdir)
    
    all_data = pd.read_csv('data_splits.csv')
    
    unique_samples = all_data['patient'].unique()
    for s in unique_samples:
       str_s = str(s)
       start = time.time()
       print(f'Sample:  {str_s}')
       sample_outdir = os.path.join(perf_outdir, str_s)
       if not os.path.exists(sample_outdir):
           os.mkdir(sample_outdir)
       os.chdir(sample_outdir)
       all_data.to_csv('all_data.csv')
       sub_df = all_data[all_data['patient'] == s] 
       sub_df.to_csv('data_splits.csv')

       print(f'Formatting input.')
       wes_pipeline_to_ml_input(sample_outdir, kit = 'unknown', \
                          indication = 'unknown', 
                          split = 'full_cohort', 
                          prefilter = True, 
                          build_splits = False)

       print(f'Integrating CNV and mutation data.')
       get_prob_somatic_mafs(os.getcwd(), max_cores = 1,\
                             output_dir = 'mafs/prob_somatic/',\
                             hyperparams = prob_somatic_hyperparams)

       print(f'One-hot encoding categorical features.')
       engineer_features(os.getcwd())
       
       if truth_labels is not None:
           # if truth labels is provided, we can evaluate accuracy
           # add truth labels
           evaluate_accuracy = True
           engineered_input = 'mafs/engineered_plus_labels.csv'
       else:
           evaluate_accuracy = False
           engineered_input = 'mafs/engineered.csv'
       print('Applying new_normal')
       try:
           use_trained_model(os.getcwd(), pickled_model = trained_model,
                         evaluate_performance = evaluate_accuracy,
                         input_features = tuple(input_features),
                         input_data = engineered_input,
                         splits = None, perf_by_sample = False)
        
           elapsed = time.time() - start
           print(f'{elapsed} seconds elapsed for patient sample {str_s}.')
       except ValueError:
           print('features missing for this patient. skipping.')
       os.chdir('../..')

#calculate features based on info_snps:
# prob_somatic or snp vaf histograms 
