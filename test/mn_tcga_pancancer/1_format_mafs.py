"""
This file builds the 'data_splits' table for TCGA data, 
that contains designated train / validation / test as well as the paths to input mutation data, copy number data. 
Authors:  Dave Masica and R. Tyler McLaughlin 
"""
import pandas as pd
import warnings
import os

warnings.filterwarnings('ignore')

# Use the directory structure to get all the info I need to define the train-test-validation splits
path = '/fioSeq/projects/input/dna/original/benchmark/tcga/'
kits = [path + kit for kit in os.listdir(path)]
indications = [kit + '/' + indication for kit in kits for indication in os.listdir(kit)]
patients = set([indication + '/' + patient.split('_tumor')[0] for indication in indications for patient in os.listdir(indication) if '_tumor' in patient])


# Make the data splits df
data = []
for patient in patients:
    kit, indication, patient = patient.split('/')[-3:]
    data.append([patient,kit,indication])
    
data_splits = pd.DataFrame(data, columns=['patient', 'kit', 'indication'])
data_splits['split'] = data_splits['kit'].apply(lambda x: 'train' if x == 'agilent_custom_v2_kit' else \
              ('validation' if x == 'hgsc_vcrome_kit' else 'test'))
data_splits['cnv'] = '/fioSeq/projects/output/dna/tcga/' + data_splits['patient'] + '/cnv/' + \
         data_splits['patient'] + '/' + data_splits['patient'] + '_tumor.cns'
data_splits['input_maf'] = '/fioSeq/projects/output/dna/tcga/' + data_splits['patient'] + '/mutation/' + \
                  data_splits['patient'] + '/' + data_splits['patient'] + '.maf'



data_splits.to_csv('data_splits_.csv', index_label=False)


#maf = pd.read_csv(data_splits.loc[0,'input_maf'], sep='\t')

if not os.path.exists('mafs/formatted'):
    os.makedirs('mafs/formatted')

# Drop the old prob_somatic and xgboost_class labels
for idx in data_splits.index:
    single_patient_output_maf_filename = 'mafs/formatted/' + data_splits.loc[idx, 'patient']+'.maf' 
    if not os.path.exists(single_patient_output_maf_filename):
        maf = pd.read_csv(data_splits.loc[idx,'input_maf'], sep='\t')
        maf.drop(columns=['info_snp', 'fioseq_prob_somatic','xgboost_class'], inplace=True)
        maf.loc[maf['filter'].isin(['alt_allele_in_normal;germline_risk', 'germline_risk;alt_allele_in_normal', 'germline_risk', 'alt_allele_in_normal']), 'filter'] = maf['filter'] + ';PASS;tumor_only'
        maf.to_csv(single_patient_output_maf_filename, index_label=False, sep='\t')
    else:
        print(single_patient_output_maf_filename + ' already exists. Skipping.')



# Finally, save data_splits to the directory 'formatted'. Downstream, these are loaded into prob_somatic.
data_splits['snv'] = os.path.dirname(os.path.abspath(__file__)) + '/mafs/formatted/' + data_splits['patient'] + '.maf'
data_splits.to_csv('data_splits.csv', index_label=False)
# delete the intermediate file.
os.remove('data_splits_.csv')

