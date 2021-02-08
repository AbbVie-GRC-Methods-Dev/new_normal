import pandas as pd
import glob




def run_to_mn_test():
    # to make sure filters look good.
    test_mn_maf = pd.read_csv('/bioinformatics/Users/mclaurt/mountables/pt23.maf', sep ='\t')

    tumor_only_input = pd.read_csv('mafs/engineered.csv')
    test_to_maf = tumor_only_input[tumor_only_input['sample'] == 'pt23']
    
    test_mn_maf['filter'].value_counts()
    
    minimal_columns = ['chrom','pos','ref','alt','filter','pop_max']
    debug_columns = minimal_columns + ['fpfilter','fioseq_prob_somatic','xgboost_class']
    
    left_joined = pd.merge(test_to_maf[minimal_columns], test_mn_maf[debug_columns], on = ['chrom','pos','ref','alt'], how = 'left', suffixes = ['_to','_mn'] )
    left_joined.filter_to.value_counts()
    # PASS;tumor_only 1112
    print('Filters for the matched normal frame.  PASS are true somatic mutations, the rest are normal mutations')
    print(left_joined.filter_mn.value_counts())
    #PASS;tumor_only         857
    #PASS                    240
    #alt_allele_in_normal     14
    #germline_risk             1
    #Name: filter_mn, dtype: int64


if __name__ == '__main__':
    tumor_only_input = pd.read_csv('mafs/engineered.csv')
    matched_normal_dir = '/fioSeq/projects/output/dna/hugo_2016/ucla_samples/mutation/'
    matched_normal_mafs = glob.glob(matched_normal_dir + '*/*.maf') 
    mn_mafs_list = []
    for maf_file in matched_normal_mafs:
        print(maf_file)
        mn_maf = pd.read_csv(maf_file, sep = '\t')
        mn_mafs_list.append(mn_maf)
    
    full_mn_maf = pd.concat(mn_mafs_list)    
    left_joined = pd.merge(tumor_only_input, full_mn_maf[['sample', 'chrom','pos','ref', 'alt', 'filter']], on = ['sample','chrom','pos','ref','alt'], how = 'left', suffixes = ['_to','_mn'] )
    left_joined['truth'] = left_joined['filter_mn'].apply(lambda x: 1 if x == 'PASS' else 0)
    left_joined.to_csv('mafs/engineered_plus_labels.csv')

