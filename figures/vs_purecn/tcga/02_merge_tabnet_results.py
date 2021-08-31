import pandas as pd

tabnet_results = '/bioinformatics/Users/mclaurt/projects/deep_learning_tumor-only_variant_calling/test/no_ps_fit_with_snp-vaf-bin/tcga_pancancer/mafs/tabnet_trained_tcga_v3_predictions_all_features.csv'

purecn_results = 'merged_fio_purecn_filtered.csv'

tn = pd.read_csv(tabnet_results)
pcn = pd.read_csv(purecn_results)
print(pcn._merge.value_counts())

scpra = ['sample','chrom','pos','ref','alt']
tn_minimal = tn[scpra + ['tabnet_pred', 'tabnet_proba_1', 'truth']]

merged = pd.merge(pcn, tn_minimal, indicator = '_tn_merge', on = scpra, how = 'left')
merged.to_csv('pure_tabnet_merge.csv')
