"""
This file loads up the MAF corresponding to the  23 patient subset from the Hugo melanoma study. 
The MAF was produced by our internal variant calling pipeline.  The MAF is formatted before adding copy number info via probSomatic.
It basically removes old classifications from our production pipeline. No matched normals were used.   

R. Tyler McLaughlin 
"""
import pandas as pd
import glob
import os

project_dir = '/fioSeq/projects/output/dna/hugo_2016/ucla_samples/'

# build data splits, containing the maf path, cnr path,
row_list = []
# e.g. project_dir/mutation/pt1/pt1.maf
for patient_maf in glob.glob(project_dir + 'mutation/*/*.maf'):
    patient_id = os.path.basename(patient_maf).split('.')[0]
    cnv_file = os.path.join(project_dir, 'cnv', patient_id, patient_id + "_tumor.cns" ) 
    row = {'patient' : patient_id, 'kit' :  'nimblegen_seqcap_v3_kit', 'indication' : 'melanoma', 'split' : 'full_cohort', 'cnv' : cnv_file, 'input_maf' : patient_maf}
    row_list.append(row)

data_splits_ = pd.DataFrame(row_list)
data_splits_.to_csv('data_splits_.csv')

if not os.path.exists('mafs/formatted'):
    os.makedirs('mafs/formatted')

# Drop the old prob_somatic and xgboost_class labels
for idx in data_splits_.index:
    single_patient_output_maf_filename = 'mafs/formatted/' + data_splits_.loc[idx, 'patient']+'.maf' 
    if not os.path.exists(single_patient_output_maf_filename):
        maf = pd.read_csv(data_splits_.loc[idx,'input_maf'], sep='\t')
        maf.drop(columns=['info_snp', 'fioseq_prob_somatic','xgboost_class'], inplace=True)
        # add ';PASS;tumor_only' so that we don't lose these germline mutations 
        maf.loc[maf['filter'].isin(['alt_allele_in_normal;germline_risk', 'germline_risk;alt_allele_in_normal', 'germline_risk', 'alt_allele_in_normal']), 'filter'] = maf['filter'] + ';PASS;tumor_only'
        maf.to_csv(single_patient_output_maf_filename, index_label=False, sep='\t')
    else:
        print(single_patient_output_maf_filename + ' already exists. Skipping.')

# Finally, save data_splits to the directory 'formatted'. Downstream, these are loaded into prob_somatic.
data_splits_['snv'] = os.path.dirname(os.path.abspath(__file__)) + '/mafs/formatted/' + data_splits_['patient'] + '.maf'
data_splits_.to_csv('data_splits.csv', index_label=False)
# delete the intermediate file.
os.remove('data_splits_.csv')
