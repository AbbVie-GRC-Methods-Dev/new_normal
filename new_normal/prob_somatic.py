import random
import time
import string
import warnings

import pandas as pd
import numpy as np
import scipy.stats as ss

warnings.filterwarnings('ignore')


class ProbSomatic:
    '''
    
    '''
    
    def __init__(self, somatic, info_snp_pop=0.025, info_snp_freq=0.95, info_snp_depth=20, 
                 min_info_snp=10, expand_copy_window_by=0.01, max_copy_window_size=2.0,
                 vaf_bin_size=0.05):
        
        self.raw = somatic
        self.info_snp_pop = info_snp_pop # info snps have a pop_max greater than this
        self.info_snp_freq = info_snp_freq # maximmum allele frequency for heterozygous SNP.
        self.info_snp_depth = info_snp_depth # this is just a qc variable for 'info_snp'
        self.min_info_snp = min_info_snp # minimum number of informative snps to calculate prob_somatic
        self.expand_copy_window_by = expand_copy_window_by # if not enough info snps in segment, expand copy number ratio by this ammount.  look in chromosome segments in this new window. 
        self.max_copy_window_size = max_copy_window_size # to increase info_snp count, what is the max change in copy-ratio threshold
        self.vaf_bin_size = vaf_bin_size # variant of interest VAF +/- freq, used to calculate probability from MLE-fit Gaussian. 
        
    def _format_data(self):
        ''' From the raw input, make an object of informative snps, 
        one for variants we are interested in classifying, and one
        for variants we are interested in classifying but are not 
        also informative snps (i.e., we still have to classify those 
        variants.)
        '''
        
        # Get all variants meeting the user-defined requirements to be an info_snp, and create a field of True info_snps
        self.info_snps = self.raw[(self.raw['t_maj_allele'] <= self.info_snp_freq) & \
                                  (self.raw['pop_max'] >= self.info_snp_pop) & \
                                 (self.raw['t_depth'] >= self.info_snp_depth)]
        self.info_snps['info_snp'] = True
        
        # Get the variants we want to classify, making a dataframe object for those we definitely need to apply prob_somatic to (unknown)
        # !!!AM I USING FPFILTER CORRECTLY HERE!!!
        self.variants = self.raw[(self.raw['filter'].str.contains('PASS')) & \
                                 (self.raw['ontology'].isin(['missense', 'frameshift_indel', 'inframe_indel', 'nonsense'])) & \
                                 (self.raw['fpfilter']=='PASS')]
                                 #((self.raw['fpfilter'].str.contains('PASS')) | (self.raw['fpfilter'].isnull()))]
        self.variants = self.variants.assign(info_snp=np.where(self.variants.index.isin(self.info_snps.index), True, False))
        self.variants.unknown = pd.concat([self.variants, self.info_snps, self.info_snps]).drop_duplicates(keep=False)
        
    def _get_segment(self, chrom, seg_start):
        ''' Private function takes starting coordinate for a copy segment,
        then builds dataframes of informative snps and variants to classify
        specific to that segment. 
        '''
        
        # Get the informative snps for this copy segment
        info_snps = self.info_snps[(self.info_snps['chrom'] == chrom) & (self.info_snps['seg_start'] == seg_start)]
        #print(info_snps.shape[0])
        info_snps['prob_somatic'], info_snps['mu'], info_snps['sigma'] = 0.0, np.nan, np.nan # don't calc stats for info_snps!
        
        # Get the variants we want to classify. Note the stats initialization! If nothing else, is this messing up visualization
        unknown = self.variants.unknown[(self.variants.unknown['chrom'] == chrom) & (self.variants.unknown['seg_start'] == seg_start)]
        unknown['info_snp'] = False
        unknown['prob_somatic'], unknown['mu'], unknown['sigma'] = 0.5, np.nan, np.nan # Initialize w/0.5 as to not bias xgboost with 0 or 1? 
        
        # Do all this to get the 'known' variants (those we want to classify, but are germline info_snps 
        variants = self.variants[(self.variants['chrom'] == chrom) & (self.variants['seg_start'] == seg_start)]
        # !! Instead do 1 - t_maj_allele !! ????
        variants['prob_somatic'], variants['mu'], variants['sigma'] = 0.0, np.nan, np.nan
        known = pd.concat([variants, unknown]).drop_duplicates(subset=['chrom', 'pos', 'ref', 'alt'], keep=False)
    
        window_size = 0.0
        if unknown.shape[0]:
        
            copy_ratio = unknown.copy_ratio.values[0]
            # If a segment has too few info snps, lets find some on the same chrom, where the copy ratio is close to targe. 
            while info_snps.shape[0] <= self.min_info_snp and window_size <= self.max_copy_window_size:
                
                window_size += self.expand_copy_window_by
                info_snps = self.info_snps[(self.info_snps['chrom'] == chrom) & \
                                           (self.info_snps['copy_ratio'].between(copy_ratio - window_size, copy_ratio + window_size))]
                    
        
        random.seed(1000 * time.time()) # Barcode all variants to keep track of info snps included by copy ratio window expansion
        barcode = ''.join(random.choices(string.ascii_letters + string.digits, k=16))
        
        # These could be useful features to track for training the downstream learner. 
        info_snps['window_size'], info_snps['support'], info_snps['barcode'] = np.nan, np.nan, barcode
        unknown['window_size'], unknown['support'], unknown['barcode'] = round(window_size, 2), info_snps.shape[0], barcode
        known['window_size'], known['support'], known['barcode'] = np.nan, np.nan, barcode
    
        return info_snps, known, unknown
    
    def _calc_prob(self):
        ''' For all variants we want to classify--that aren't deemed germline
        by virtue of also being informative snps--compute the probability of
        being somatic and ultimately make the prob somatic dataframe object. 
        '''
        
        dfs = [] # These need to be computed one segment at a time
        segments = pd.concat([self.info_snps, self.variants])
        segments = segments[['chrom', 'seg_start']].drop_duplicates().sort_values(['chrom', 'seg_start'])
        for segment in segments.index:
            
            chrom, seg_start = segments.loc[segment].values
            info_snps, known, unknown = self._get_segment(chrom, seg_start)
            dfs.append(pd.concat([info_snps, known]).drop_duplicates()) # Nothing else to do with the informative snps. 
            #dfs.append(known) # Nothing left to do becuase we've already said your germline. 
            
            # If there are any unknown variants left in this segment, calc prob somatic
            if unknown.shape[0]:
                # if enough info snps, calculate prob somatic the proper way
                if info_snps.shape[0] >= self.min_info_snp: 
                
                    # Get shape parameters from a maximum liklihood estimate (MLE)
                    mu, sigma = ss.norm.fit(info_snps.t_maj_allele)

                    # Compute probability if unknown VAFs belonging to a Gaussian defined by the MLE-derived shape paramters
                    prob_germ = ss.norm(mu, sigma).cdf(unknown.t_maj_allele + self.vaf_bin_size) - \
                        ss.norm(mu, sigma).cdf(unknown.t_maj_allele - self.vaf_bin_size)
                    
                    # Convert the above-calulated probability of being germline to probability of being somatic. 
                    unknown['prob_somatic'] = 1 - abs(prob_germ) 
                    # Store these shape parameters for training machine-learning models
                    unknown['mu'] = mu
                    unknown['sigma'] = sigma
                # if unknown variant exist in this segment, but not enough info snps were found, we should still output the variant.
                # this is an edge case bug fix.
                else:  
                    unknown['prob_somatic'], unknown['mu'], unknown['sigma'] = 0.5, 0.5, 1 
                
                # append unknown variants if there were any.
                dfs.append(unknown)
        
        self.prob_somatic = pd.concat(dfs)
        self.prob_somatic.loc[self.prob_somatic['info_snp'], 'prob_somatic'] = 0.0 # !!! I DONT KNOW HOW THESE ARE SLIPPING THRU !!!
        self.prob_somatic.loc[self.prob_somatic['t_maj_allele'] > self.info_snp_freq, 'prob_somatic'] = 0.0


def merge_somatic(paths):
    ''' For a given sample, merge all copy-number and variant info by chromosome and position.
    '''

    snv = pd.read_csv(paths[0], sep='\t')
    cnv = pd.read_csv(paths[1], sep='\t')

    chroms = set(snv['chrom']).intersection(set(cnv['chromosome']))
    dfs = []
    for chrom in chroms:
        cnv_chrom = cnv[cnv['chromosome'] == chrom]
        intervals = cnv_chrom.loc[:, 'start':'end'].apply(tuple, 1).tolist()
        idx = pd.IntervalIndex.from_tuples(intervals, closed='neither')

        snv_chrom = snv[snv['chrom'] == chrom]
        snv_chrom['seg_chr'] = cnv_chrom.iloc[idx.get_indexer(snv_chrom.pos),:]['chromosome'].values
        snv_chrom['seg_start'] = cnv_chrom.iloc[idx.get_indexer(snv_chrom.pos),:]['start'].values
        snv_chrom['seg_end'] = cnv_chrom.iloc[idx.get_indexer(snv_chrom.pos),:]['end'].values
        snv_chrom['copy_ratio'] = 2**(cnv_chrom.iloc[idx.get_indexer(snv_chrom.pos),:]['log2'].values)
        snv_chrom['copy_depth'] = cnv_chrom.iloc[idx.get_indexer(snv_chrom.pos),:]['depth'].values
        snv_chrom['probes'] = cnv_chrom.iloc[idx.get_indexer(snv_chrom.pos),:]['probes'].values
        snv_chrom['weight'] = cnv_chrom.iloc[idx.get_indexer(snv_chrom.pos),:]['weight'].values

        dfs.append(snv_chrom)
    out = pd.concat(dfs)
    out.to_csv('merge_somatic.csv')
    return out
