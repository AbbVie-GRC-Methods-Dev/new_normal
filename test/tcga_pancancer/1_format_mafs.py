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
    
df = pd.DataFrame(data, columns=['patient', 'kit', 'indication'])
df['split'] = df['kit'].apply(lambda x: 'train' if x == 'agilent_custom_v2_kit' else ('validation' if x == 'hgsc_vcrome_kit' else 'test'))
df['cnv'] = '/bioinformatics/Users/masicdx/tcga/cnv/' + df['patient'] + '.cns'
df['input_maf'] = '/bioinformatics/Users/masicdx/tcga/maf/' + df['patient'] + '.maf'
df['input_maf'] = '/fioSeq/projects/output/dna/tcga/' + df['patient'] + '/mutation/' + df['patient'] + '/' + df['patient'] + '.maf'


#df.to_csv('data_splits.csv', index_label=False)


#maf = pd.read_csv(df.loc[0,'input_maf'], sep='\t')

if not os.path.exists('mafs/formatted'):
    os.mkdir('mafs/formatted')

# Drop the old prob_somatic and xgboost_class labels
for idx in df.index:
    maf = pd.read_csv(df.loc[idx,'input_maf'], sep='\t')
    maf.drop(columns=['info_snp', 'fioseq_prob_somatic','xgboost_class'], inplace=True)
    maf.loc[maf['filter'].isin(['alt_allele_in_normal;germline_risk', 'germline_risk;alt_allele_in_normal', 'germline_risk', 'alt_allele_in_normal']), 'filter'] = maf['filter'] + ';PASS;tumor_only'
    maf.to_csv(df.loc[idx, 'patient']+'.maf', index_label=False, sep='\t')


# Check one of them. Then I'll move the original mafs to the directory 'sentieon' 
pd.read_csv('TCGA-ZM-AA0N.maf', sep='\t')


# Finally, save data_splits to the directory 'formatted'. Downstream, these are loaded into prob_somatic.
df['snv'] = '/bioinformatics/Users/masicdx/tcga/maf/formatted/' + df['patient'] + '.maf'
df.to_csv('data_splits.csv', index_label=False)


# Check it!
#pd.read_csv('data_splits.csv')

