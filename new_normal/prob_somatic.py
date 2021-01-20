import math
from math import pi, exp, sqrt
import string
import time
import random
import warnings
warnings.filterwarnings('ignore')

import pandas as pd
import numpy as np
from scipy.stats import norm
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from sklearn.metrics import roc_auc_score, roc_curve, average_precision_score, precision_recall_curve

def get_segments(df):
    ''' Given a dataframe that at least contains a "chrom" and "seg_start" field,
        return a dataframe that is just the none-redundant and sorted coordinates
        for each copy segment start point. 
    '''
    
    return df[['chrom', 'seg_start']].drop_duplicates().sort_values(['chrom', 'seg_start'])

class ProbSomatic:
    '''
    
    '''
    
    def __init__(self, cohort, 
                 info_snp_pop=0.025, info_snp_freq=0.95, info_snp_depth=20, 
                 min_info_snp=10, expand_copy_window_by=0.01, max_copy_window_size=2.0,
                 vaf_bin_size=0.05):
        
        self.cohort = cohort
        self.info_snp_pop = info_snp_pop # this defines the 'snp' part of 'info_snp'
        self.info_snp_freq = info_snp_freq # this defines 'info' part of 'info_snp'
        self.info_snp_depth = info_snp_depth # this is just a qc variable for 'info_snp'
        self.min_info_snp = min_info_snp # minimum number of informative snps
        self.expand_copy_window_by = expand_copy_window_by # to increase info_snp count, how quickly should we lower the copy-ratio threshold 
        self.max_copy_window_size = max_copy_window_size # to increase info_snp count, what is the max change in copy-ratio threshold
        self.vaf_bin_size = vaf_bin_size # variant of interest +/- freq, used to calculate probability 
        self.prob_somatic = pd.DataFrame()
    
    def snv(self, sample):
        ''' Read in a fioSeq maf.
        '''
        return pd.read_csv(self.cohort.loc[sample, 'snv'], sep='\t')
    
    def cnv(self, sample):
        ''' Read in a CNVKit call.cns file.
        '''
        
        return pd.read_csv(self.cohort.loc[sample, 'cnv'], sep='\t')
            
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
        self.variants = self.raw[(self.raw['ontology'].isin(['missense', 'frameshift_indel', 'inframe_indel', 'nonsense'])) & \
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
        info_snps['prob_somatic'], info_snps['mu'], info_snps['sigma'] = 0.0, np.nan, np.nan # don't calc stats for info_snps!
        
        # Get the variants we want to classify. Note the stats initialization! If nothing else, is this messing up visualization
        unknown = self.variants.unknown[(self.variants.unknown['chrom'] == chrom) & (self.variants.unknown['seg_start'] == seg_start)]
        unknown['info_snp'] = False
        unknown['prob_somatic'], unknown['mu'], unknown['sigma'] = 0.0, np.nan, np.nan # DEFAULT PROB 0?! Should this always be the case
        
        # Do all this to get the 'known' variants (those we want to classify, but are germline info_snps 
        variants = self.variants[(self.variants['chrom'] == chrom) & (self.variants['seg_start'] == seg_start)]
        # !! Instead do 1 - t_maj_allele !!
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
        segments = get_segments(pd.concat([self.info_snps, self.variants]))
        for segment in segments.index:
            
            chrom, seg_start = segments.loc[segment].values
            info_snps, known, unknown = self._get_segment(chrom, seg_start)
            dfs.append(pd.concat([info_snps, known]).drop_duplicates()) # Nothing else to do with the informative snps. 
            #dfs.append(known) # Nothing left to do becuase we've already said your germline. 
            
            # If there are any unknown variants left in this segment, calc prob somatic
            if unknown.shape[0] and info_snps.shape[0] >= self.min_info_snp:
                
                # Get shape parameters from a maximum liklihood estimate (MLE)
                mu, sigma = norm.fit(info_snps.t_maj_allele)

                # Compute probability if unknown VAFs belonging to a Gaussian defined by the MLE-derived shape paramters
                prob_germ = norm(mu, sigma).cdf(unknown.t_maj_allele + self.vaf_bin_size) - \
                    norm(mu, sigma).cdf(unknown.t_maj_allele - self.vaf_bin_size)
                
                # Convert the above-calulated probability of being germline to probability of being somatic. 
                unknown['prob_somatic'] = 1 - abs(prob_germ) 
                # Store these shape parameters for training machine-learning models
                unknown['mu'] = mu
                unknown['sigma'] = sigma
                
                
                
                dfs.append(unknown)
        
        prob_somatic = pd.concat(dfs)
        prob_somatic.loc[prob_somatic['info_snp'], 'prob_somatic'] = 0.0 # !!! I DONT KNOW HOW THESE ARE SLIPPING THRU !!!
        prob_somatic.loc[prob_somatic['t_maj_allele'] > self.info_snp_freq, 'prob_somatic'] = 0.0
        
        # For now, only using these fields for the prob somatic dataframe object (should I use others? max_cosmic_count?)
        #columns = ['sample', 'chrom', 'pos', 'qual', 'filter', 'ontology',
        #           't_depth', 't_maj_allele', 'max_cosmic_count', 'pop_max', 'fpfilter', 'seg_chr', 'seg_start',
        #           'copy_ratio', 'info_snp','prob_somatic', 'mu', 'sigma', 'window_size', 'support', 'barcode']
        self.prob_somatic = pd.concat([self.prob_somatic, prob_somatic]) #[columns]])
    
    def add_sample(self, sample):
        ''' For a given sample, merge all copy-number and variant info by chromosome and position.
        '''
        
        snv, cnv = self.snv(sample), self.cnv(sample)
        
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

            # Only include variants/mutations that actually fell into a copy-number seg/bin !!! BAD IDEA. COULD THROW OUT VARIANTS!!!
            #snv_chrom = snv_chrom[snv_chrom['pos'].between(snv_chrom['seg_start'], snv_chrom['seg_end'])]
            #dfs.append(snv_chrom[['sample', 'chrom', 'pos', 'ref', 'alt', 'qual', 'filter',
            #               'ontology', 't_alt_count', 't_ref_count', 't_depth', 't_alt_freq',
            #               't_ref_freq', 't_maj_allele', 'max_cosmic_count', 'pop_max', 'fpfilter', 'seg_chr', 
            #               'seg_start', 'seg_end', 'copy_ratio', 'copy_depth', 'probes', 'weight']])
            dfs.append(snv_chrom)

        self.raw = pd.concat(dfs)
        self._format_data()
        self._calc_prob()
        
    def add_samples(self, samples):
        ''' Add multiple samples to the dataset at once (see add_sample for details)
        '''
        
        info_snps = []
        for sample in samples: 
            self.add_sample(sample)
            info_snps.append(self.info_snps)

        self.info_snps = pd.concat(info_snps)

def normal_curve(mu, sigma):
    
    #X,y = np.arange(0, 1.01, 0.01),[]
    X,y = np.arange(mu - (3*sigma), mu + (3*sigma), 0.01), []
    for x in X:
        
        y.append((1/((sigma)*sqrt(2*pi)))*exp(-0.5*((x - mu)/sigma)**2))
        
    
    return pd.DataFrame([X,y]).T.rename(columns={0:'X',1:'y'})

def make_fig(df, index):
    
    segments = get_segments(df[~df['prob_somatic'].isnull()])

    chrom, seg_start = segments.iloc[index][0], segments.iloc[index][1]
    segment = df[(df['chrom'] == chrom) & (df['seg_start'] == seg_start)]
    
    try:
        barcode = segment[~segment.prob_somatic.isnull()]['barcode'].values[0]
        stats = segment[~segment['mu'].isnull()].iloc[0]
        info_snps = df[(df['chrom'] == chrom) & (df['barcode'] == barcode) & (df['info_snp'])].round(3) #.isnull())].round(3)
        unknown = segment[~segment['prob_somatic'].isnull()].round(3)
        positive = unknown[unknown['filter']=='PASS'].round(3)
        negative = unknown[unknown['filter']!='PASS'].round(3)
    except IndexError:
        print('No mutations for chromosome {}, start position {}'.format(chrom, seg_start))
        return None, None
    
    fig = make_subplots(rows=2, cols=1, shared_xaxes=True)
    fig.add_trace(px.histogram(info_snps, x='t_maj_allele', range_x=(0.5, 1.0), nbins=10, color_discrete_sequence=['green']).data[0], row=1, col=1)
    fig.add_trace(px.line(normal_curve(stats.mu, stats.sigma), x='X', y='y', range_x=(0.5, 1.0), color_discrete_sequence=['green']).data[0], row=1, col=1)

    if positive.shape[0]:
        fig.add_trace(px.scatter(positive, x='t_maj_allele', y='copy_ratio', size=positive['prob_somatic']+0.01, color_discrete_sequence=['orange']).data[0], row=2, col=1)
        
    if negative.shape[0]:
        fig.add_trace(px.scatter(negative, x='t_maj_allele', y='copy_ratio', size=negative['prob_somatic']+0.01, color_discrete_sequence=['blue']).data[0], row=2, col=1)
        
    fig.add_trace(px.scatter(info_snps, x='t_maj_allele', y='copy_ratio', color_discrete_sequence=['green'], symbol_sequence=[4]).data[0], row=2, col=1)
    
    # Solving px add_trace issues 
    fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers', marker=dict(size=10, color='orange'), legendgroup='Somatic', showlegend=True, name='Somatic'))
    fig.add_trace(go.Scatter(x=[None], y=[None], mode='markers', marker=dict(size=10, color='blue'), legendgroup='Germline', showlegend=True, name='Germline'))
    fig.data[2]['marker']['sizeref'] = 0.001
    fig.data[3]['marker']['sizeref'] = 0.001
    
    fig.update_layout(title='Chromosome={}, Start position={}, Window size={}, Support={}'.format(chrom, seg_start, stats.window_size, stats.support), yaxis_title='count', xaxis2_title='allele frequency', yaxis2_title='copy ratio')
    fig.update_xaxes(range=[0.5, 1.0])

    return unknown, fig

def performance_plots(df):
    '''
    '''

    # First convert classes to binary
    #df['classes'] = df['filter'].map({'PASS':1,'PASS;tumor_only':0}) # Not working now that we are including germline risk, etc.
    df['classes'] = df['filter'].apply(lambda x: 1 if x == 'PASS' else 0)
    
    # Get the roc/auc data
    roc = pd.DataFrame(roc_curve(df['classes'], df['prob_somatic'])).T.rename(columns={0:'false positive rate',1:'true positive rate (recall)',2:'prob. somatic'})
    roc[roc > 1.0] = 1.0
    roc = roc.round(2)
    auc = round(roc_auc_score(df['classes'], df['prob_somatic']), 3)
    fpr, sensitivity, roc_prob = roc.iloc[((1 - roc['false positive rate']) +  roc['true positive rate (recall)']).idxmax()] 
    specificity = round(1 - fpr, 2)
    balanced_accuracy = round((sensitivity + specificity)/2, 2)
    sensitivity, specificity, balanced_accuracy, roc_prob

    # Get the precision-recall data (prc) and avarage prc
    prc = pd.DataFrame(precision_recall_curve(df['classes'], df['prob_somatic'])).T.rename(columns={0:'precision',1:'recall',2:'prob. somatic'})
    prc[prc > 1.0] = 1.0
    prc = prc.round(2)
    aps = round(average_precision_score(df['classes'], df['prob_somatic']), 3)
    precision, recall, prc_prob = prc.iloc[(prc['precision'] +  prc['recall']).idxmax()] 
    f1_score = round(2*(precision * recall)/(precision + recall), 2)
    precision, recall, f1_score, prc_prob

    fig = make_subplots(rows=1, cols=2, shared_yaxes=True, subplot_titles = \
        ['bal acc = {}; sens = {}, spec = {}, prob = {}'.format(balanced_accuracy, sensitivity, specificity, roc_prob), \
         'f1 = {}; precision = {}, recall = {}, prob = {}'.format(f1_score, precision, recall, prc_prob)])
    for i in fig['layout']['annotations']:
        i['font'] = dict(size=10,color='#ff0000')

    fig.add_trace(px.scatter(roc, x='false positive rate', y='true positive rate (recall)', color='prob. somatic', title='auc = {}'.format(auc)).data[0], row=1, col=1)
    fig.add_trace(px.scatter(prc, x='recall', y='precision', color='prob. somatic', title='average precision score = {}'.format(aps)).data[0], row=1, col=2)
    fig.update_layout(title='continuous-value performance: roc auc = {}; average precision score = {}'.format(auc, aps), \
       xaxis_title='false positive rate', yaxis_title='true positive rate (recall)', \
            xaxis2_title='recall', yaxis2_title='precision', \
               width=1000, height=350)

    return fig, auc, aps

