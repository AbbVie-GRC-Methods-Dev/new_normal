import sys
import os
import glob
import pandas as pd

# strategy:  merge formatted mafs with individual purecn variant outputs.

purecn_results_dir = '/bioinformatics/Users/mclaurt/projects/deep_learning_tumor-only_variant_calling/purecn/with_cnvkit_seg/purecn_results/'

fmafs_tcga = glob.glob('../../../test/no_ps_fit_with_snp-vaf-bin/tcga_pancancer/mafs/formatted/*.maf')
fmafs_hugo = glob.glob('/bioinformatics/Users/mclaurt/projects/deep_learning_tumor-only_variant_calling/test/05_13_earlier_dm-splits/with_snp-vaf-bin/hugo_melanoma/mafs/formatted/*.maf')

scpra = ['sample','chrom','pos','ref','alt']
purecn_scpra = ['Sampleid','chr','start','REF','ALT']

df_list = []
dataset = sys.argv[1]
if dataset == 'tcga':
    fmafs = fmafs_tcga
elif dataset == 'hugo':
    fmafs = fmafs_hugo

if not os.path.exists('purecn_only_vars'):
   os.mkdir('purecn_only_vars')

for i, fmaf in enumerate(fmafs):
   fmaf_df = pd.read_csv(fmaf, sep = '\t')
   patient_id = os.path.basename(fmaf).split('.')[0]
   print(i+1, patient_id)
   pcn_vars = purecn_results_dir + str(patient_id) + '_variants.csv'
   print(pcn_vars)
   try:
      pcn_df = pd.read_csv(pcn_vars)
      merged = pd.merge(fmaf_df, pcn_df, left_on  = scpra, right_on = purecn_scpra, how = 'outer',  indicator = True)
      # get purecn-only variants.  only after filtering will we care about fioSeq only variants.
      purecn_only_variants = merged[merged._merge == 'right_only']
      purecn_only_variants.to_csv(f'purecn_only_vars/{patient_id}.csv')
      merged = pd.merge(fmaf_df, pcn_df, left_on  = scpra, right_on = purecn_scpra, how = 'outer', indicator = True)
      print(merged._merge.value_counts())
      # this part will discard the purecn--only variants, hence why they were saved above.
      merged = merged[(merged['filter'].str.contains('PASS')) & (merged['ontology'].isin(['missense', 'frameshift_indel', 'inframe_indel', 'nonsense'])) & (merged['fpfilter'] == 'PASS') & (merged['pop_max'] < 0.01)]
      df_list.append(merged)
   except FileNotFoundError:
      print(f'patient {patient_id} not found in PureCN results.')
      continue 

out = pd.concat(df_list)
out.to_csv('merged_fio_purecn_filtered.csv')
