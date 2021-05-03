#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import os
import warnings

warnings.filterwarnings('ignore')

def merge_data(snv, cnv):
    chroms = set(snv['chrom']).intersection(set(cnv['chromosome']))
    dfs = []
    for chrom in chroms:
        cnv_chrom = cnv[cnv['chromosome'] == chrom]
        v = cnv_chrom.loc[:, 'start':'end'].apply(tuple, 1).tolist()
        idx = pd.IntervalIndex.from_tuples(v, closed='neither')
        
        snv_chrom = snv[snv['chrom'] == chrom]
        snv_chrom['seg_chr'] = cnv_chrom.iloc[idx.get_indexer(snv_chrom.pos),:]['chromosome'].values
        snv_chrom['seg_start'] = cnv_chrom.iloc[idx.get_indexer(snv_chrom.pos),:]['start'].values
        snv_chrom['seg_end'] = cnv_chrom.iloc[idx.get_indexer(snv_chrom.pos),:]['end'].values
        snv_chrom['copy_ratio'] = 2**(cnv_chrom.iloc[idx.get_indexer(snv_chrom.pos),:]['log2'].values)
        snv_chrom['copy_depth'] = cnv_chrom.iloc[idx.get_indexer(snv_chrom.pos),:]['depth'].values
        snv_chrom['probes'] = cnv_chrom.iloc[idx.get_indexer(snv_chrom.pos),:]['probes'].values
        snv_chrom['weight'] = cnv_chrom.iloc[idx.get_indexer(snv_chrom.pos),:]['weight'].values

        # Only include variants/mutations that actually fell into a copy-number seg/bin
        snv_chrom = snv_chrom[snv_chrom['pos'].between(snv_chrom['seg_start'], snv_chrom['seg_end'])]
        dfs.append(snv_chrom)

    df = pd.concat(dfs)

    return df

def engineer_features(project_dir, use_purecn_purity):
    data = pd.read_csv(os.path.join(project_dir,'data_splits.csv'))
    data['snv'] = data['snv'].str.replace('formatted/', 'prob_somatic/')
    data.set_index('patient', inplace=True)
    
    dfs = []
    idx = 0
    for patient in data.index:
        # do exception handling, because in the rare case, 
        # after dropping SNPs from the prob somatic subroutine, 
        # the patient won't have variants.
        try:
            maf = pd.read_csv(data.loc[patient]['snv'], sep='\t')
            cnv = pd.read_csv(data.loc[patient]['cnv'], sep='\t')
            df = merge_data(maf,cnv)
            dfs.append(df[(df['filter'].str.contains('PASS')) & (df['ontology'].isin(['missense', 'frameshift_indel', 'inframe_indel', 'nonsense'])) &                   (df['fpfilter']=='PASS')])
                          #((df['fpfilter'].str.contains('PASS')) | (df['fpfilter'].isnull()))])                 
            idx += 1
            print(patient, idx)
        except FileNotFoundError:
            print(f'Sample {patient} not found. Probably has no variants after prob_somatic informatic SNP filtering.  Skipping.')
    
    df = pd.concat(dfs)
    df.reset_index(inplace=True)
    
    # all the features we're going to keep (and maybe one-hot-encode)
    engineered = df[['sample', 'chrom', 'pos', 'id', 'ref', 'alt', 'qual', 'filter', 'gene', 't_depth', 't_maj_allele', 'window_size', 't_alt_freq',
                     'support', 'prob_somatic', 'mu', 'sigma', 'seg_chr', 'seg_start', 'seg_end', 'copy_ratio', 'copy_depth', 'probes', 'weight', 
                     'max_cosmic_count', 'pop_max']]
    
    engineered = engineered.merge(pd.get_dummies(df['ontology'], drop_first=True), how='inner', left_index=True, right_index=True)
    engineered['100mer_mappability'] = df['100mer_mappability']
    # add count feature -- the number of mutations for each sample.
    engineered['count'] = engineered.groupby('sample')['sample'].transform('count')

    engineered = engineered.merge(pd.get_dummies(df['trinucleotide_context'], drop_first=True), how='inner', left_index=True, right_index=True)
    engineered = engineered.merge(pd.get_dummies(df['mutation_change'], drop_first=True), how='inner', left_index=True, right_index=True)
    if use_purecn_purity:
        purity = pd.read_csv(os.path.join(project_dir, 'purecn_purity.csv'))
        engineered = engineered.merge(purity, how='left', left_on='sample', right_on='sample')
    engineered['truth'] = engineered['filter'].apply(lambda x: 1 if x == 'PASS' else 0)
    engineered = engineered.merge(data[['split']], how='left', left_on='sample', right_index=True)
    # drop any duplicates remaining.
    engineered = engineered.drop_duplicates(subset = ['sample','chrom','pos','ref','alt'])
    
    print(engineered['truth'].value_counts())
    print(engineered['filter'].value_counts())
    engineered_output_file = os.path.join(project_dir, 'mafs/', 'engineered.csv')
    print(f'Saving engineered.csv to {engineered_output_file}')
    engineered.to_csv(engineered_output_file, index_label=False)
