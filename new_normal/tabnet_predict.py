import os
import sys

# sets the path for the new_normal package.
import context
 
from new_normal.fioseq_reformat import fioseq_to_ml_input
from new_normal.prob_somatic_helpers import get_prob_somatic_mafs
from new_normal.feature_engineering import engineer_features
from new_normal.model_predict_and_eval import use_trained_model                    

# features and models used in the paper.
paper_tcga_training_dir = '/bioinformatics/Users/mclaurt/projects/deep_learning_tumor-only_variant_calling/test/no_ps_fit_with_snp-vaf-bin/tcga_pancancer/'
     
# tcga-training-dependent hyperparams and features
sys.path.append(paper_tcga_training_dir)
from train_features import prob_somatic_hyperparams, vaf_bin_mode_features
    
paper_trained_model = paper_tcga_training_dir + 'tabnet_trained_tcga_v3.pkl'
 

def tabnet_predict(fio_dir, kit_name, indication_name, \
                       tcga_training_dir = paper_tcga_training_dir,\
                       prob_somatic_hyperparams = prob_somatic_hyperparams,\
                       tabnet_features = vaf_bin_mode_features,\
                       trained_model = paper_trained_model): 
    """
    Uses TCGA-trained tabnet to make somatic/germline predictions on samples.
    Does not evaluate performance, but saves results in mafs/tabnet_predictions.
    kit_name does not determine PoN or anything, 
    it's just a name that goes in datasplits.   Same with indication.
    All CNV must be in the same fioSeq input folder as mutation dir.
    Combintions of the optional arguments have not been tested and may not work.
    """

    print('Building datasplits and dropping alt chroms and structural variants from tumor-only MAFs.')
    fioseq_to_ml_input(project_dir = fio_dir, kit = kit_name, indication = indication_name, prefilter = True)
     
    working_dir = os.getcwd()
    get_prob_somatic_mafs(working_dir, max_cores = 64, hyperparams = prob_somatic_hyperparams) 
     
    print(f'Engineering features.')
    print(f'Feaatures used: {tabnet_features}.')
    # generates engineered.csv
    engineer_features(working_dir, use_purecn_purity = False)
    
    print(f'Using trained model {trained_model}.')
    # using trained model
    use_trained_model(working_dir, pickled_model = trained_model, evaluate_performance = False, input_features = tuple(tabnet_features), input_data = 'mafs/engineered.csv', splits = None) 
    print('Done!')
