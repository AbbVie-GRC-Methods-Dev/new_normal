sens_features = ['t_alt_freq', 't_maj_allele',  'count', 'max_cosmic_count', 'inframe_indel', 'missense', 'nonsense', 'ACC', 'ACG', 'ACT', 'ATA', 'ATC', 'ATG', 'ATT', 'CCA', 'CCC', 'CCG', 'CCT', 'CTA', 'CTC', 'CTG', 'CTT', 'GCA', 'GCC', 'GCG', 'GCT', 'GTA', 'GTC', 'GTG', 'GTT', 'TCA', 'TCC', 'TCG', 'TCT', 'TTA', 'TTC', 'TTG', 'TTT', 'non-SBS_x', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G', 'non-SBS_y', 'pop_max', 'purity']

sens_features_no_purity = sens_features[:-1]

vaf_bin_mode_features = sens_features_no_purity + ['snp_vaf_bin_' + str(x).zfill(2) for x in range(0,20)]

prob_somatic_hyperparams = {'info_snp_pop' : 0.01, 
                            'info_snp_freq' : 0.95,
                            'info_snp_depth' : 10.0, 
                            'min_info_snp' : 20.0,
                            'expand_copy_window_by' : 0.01,
                            'max_copy_window_size' : 1.0,
                            'vaf_bin_size' : 0.01,
                            'snp_vaf_bin_mode' : True}
