import pandas as pd


UCSC_SIMPLE_REPEATS_BED = "/unitynfs1/projects/mclaurt/util/UCSC_genomic_tracks/hg38_simpleRepeats.bed"

b = pd.read_csv(UCSC_SIMPLE_REPEATS_BED, sep = '\t', comment='t', header=None)

header = ['chrom', 'start', 'end', 'name'] 

b.columns = header

b = b.dropna(subset = ['chrom', 'start','end'])


b['start'] = b['start'].astype(int)
b['end'] = b['end'].astype(int)

print(b.head())

p = pd.read_csv('pure_tabnet_merge.csv')

repeats = []
for chrom_x in b.chrom.unique():
    print(chrom_x)
    for row in b[b.chrom == chrom_x].itertuples():
        repeats_found = p[(p.chrom == chrom_x) & (p.pos > row.start) & (p.pos <= row.end)]
        if repeats_found.shape[0] > 0:
            repeats.append(repeats_found)

r = pd.concat(repeats)

r.to_csv('repeats.csv')
