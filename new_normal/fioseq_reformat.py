import glob
import os
import pandas as pd

def fioseq_to_ml_input(project_dir, kit = 'unknown', indication = 'unknown', split = 'full_cohort', prefilter = False):
    """
    Makes simple data_splits table and saves mafs without old classification data.  Prefiltering removes non-standard chromosomes (keeps chroms 1 thru 22 + X) and structural variants. 
    """
    # build data splits, containing the maf path, cnr path,
    row_list = []
    # e.g. project_dir/mutation/pt1/pt1.maf
    for patient_maf in glob.glob(project_dir + 'mutation/*/*.maf'):
        patient_id = os.path.basename(patient_maf).split('.')[0]
        cnv_file = os.path.join(project_dir, 'cnv', patient_id, patient_id + "_tumor.cns" ) 
        row = {'patient' : patient_id, 'kit' :  kit, 'indication' : indication, 'split' : split, 'cnv' : cnv_file, 'input_maf' : patient_maf}
        row_list.append(row)
    
    data_splits_ = pd.DataFrame(row_list)
    data_splits_.to_csv('data_splits_.csv')
    
    # Example row of desired output:
    #TCGA-DX-A1KZ,nimblegen_seqcap_ez_v3_kit,SARC,test,/fioSeq/projects/output/dna/tcga/TCGA-DX-A1KZ/cnv/TCGA-DX-A1KZ/TCGA-DX-A1KZ_tumor.cns,/fioSeq/projects/output/dna/tcga/TCGA-DX-A1KZ/mutation/TCGA-DX-A1KZ/TCGA-DX-A1KZ.maf,/bioinformatics/Users/mclaurt/projects/deep_learning_tumor-only_variant_calling/test/tcga_pancancer/mafs/formatted/TCGA-DX-A1KZ.maf
    
    if not os.path.exists('mafs/formatted'):
        os.makedirs('mafs/formatted')
    
    # for the optional prefilter step
    standard_chroms = ['chr' + str(x) for x in list(range(1,23)) + ['X'] ]

    ## Drop the old prob_somatic and xgboost_class labels
    for idx in data_splits_.index:
        single_patient_output_maf_filename = 'mafs/formatted/' + data_splits_.loc[idx, 'patient']+'.maf' 
        if not os.path.exists(single_patient_output_maf_filename):
            maf = pd.read_csv(data_splits_.loc[idx,'input_maf'], sep='\t', low_memory = False)
            maf.drop(columns=['info_snp', 'fioseq_prob_somatic','xgboost_class'], inplace=True)
            if prefilter:
                maf = maf[maf['chrom'].isin(standard_chroms)]
                maf = maf[maf['ontology'] != 'structural_variant']
            maf.to_csv(single_patient_output_maf_filename, index_label=False, sep='\t')
        else:
            print(single_patient_output_maf_filename + ' already exists. Skipping.')
    
    # Finally, save data_splits to the directory 'formatted'. Downstream, these are loaded into prob_somatic.
    data_splits_['snv'] = os.getcwd() + '/mafs/formatted/' + data_splits_['patient'] + '.maf'
    data_splits_.to_csv('data_splits.csv', index_label=False)
    # delete the intermediate file.
    os.remove('data_splits_.csv')
