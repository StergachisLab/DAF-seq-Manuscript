import csv
import gzip
import pandas as pd
import numpy as np
from itertools import combinations


# Aggregate counts from peak isoforms --------------

iso_counts = 'giab_hg002_isoseq_mapped.collapsed.flnc_count.txt'
count_dict = dict()
with open(iso_counts) as fr:
    reader = csv.reader(fr)
    header = next(reader)
    for line in reader:
        count_dict[line[0]] = int(line[1])

iso_map = 'isoseq_counts_giab_hg002_Peaks_bedmap.bed.gz'
peak_dict = dict()
with gzip.open(iso_map, 'rt') as fr:
    for line in fr:
        peak,iso_delim = line.strip().split('|')
        peak2 = peak.split('\t')
        peak_id = '-'.join(peak2[:3])
        if peak_id not in peak_dict.keys():
            peak_dict[peak_id] = {'chrom':peak2[0], 'start':int(peak2[1]), 'end':int(peak2[2]), 'nreads':0}
        isoforms = iso_delim.split(';')
        for iso in isoforms:
            peak_dict[peak_id]['nreads'] += count_dict[iso]


iso_df = pd.DataFrame.from_dict(peak_dict, orient='index')

iso_df['exp_bin'] = np.log2(iso_df['nreads']).apply(np.ceil)
# collapse bins > 10 into a single bin
idx = iso_df[iso_df['exp_bin'] > 10].index
iso_df.loc[idx, 'exp_bin'] = 11
iso_df['ID'] = iso_df["chrom"].astype(str) + ':' + iso_df["start"].astype(str) + '-' + iso_df["end"].astype(str)

daf_df = pd.read_csv('../scDAF_figures/scDAF_FIRE_actuation_MSP150_withProximal.tsv', delimiter="\t")
daf_df_filt = daf_df[daf_df['is_proximal'] == True]

# Map 'value' from df1 to df2 based on 'key'
daf_df_filt['exp_bin'] = daf_df_filt['ID'].map(iso_df.set_index('ID')['exp_bin'])
daf_df_filt['exp_bin'] = daf_df_filt['exp_bin'].fillna(-1)


# Jaccard Distance by expression bin

def jaccard_similarity(set1, set2):
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return(intersection/union)

def actuation_jaccard_dist_prox_bins(s1, s2, proximal, bin):
    bin_df = daf_df_filt[daf_df_filt['exp_bin'] == float(bin)]
    h1_df = bin_df[(bin_df[f'{s1}_H1'] >= 0) & (bin_df[f'{s2}_H1'] >= 0) & (bin_df['is_proximal'] == proximal)]
    s1_act_H1 = set([f'H1_{idx}' for idx in h1_df[h1_df[f'{s1}_H1'] > 0].index])
    s2_act_H1 = set([f'H1_{idx}' for idx in h1_df[h1_df[f'{s2}_H1'] > 0].index])
    h2_df = bin_df[(bin_df[f'{s1}_H2'] >= 0) & (bin_df[f'{s2}_H2'] >= 0) & (bin_df['is_proximal'] == proximal)]
    s1_act_H2 = set([f'H2_{idx}' for idx in h2_df[h2_df[f'{s1}_H2'] > 0].index])
    s2_act_H2 = set([f'H2_{idx}' for idx in h2_df[h2_df[f'{s2}_H2'] > 0].index])
    # merge datapoint (both haplotypes)
    s1_act_all = s1_act_H1.union(s1_act_H2)
    s2_act_all = s2_act_H1.union(s2_act_H2)
    if len(s1_act_all.union(s2_act_all)) > 0:
        distance = 1 - jaccard_similarity(s1_act_all, s2_act_all)
    else:
        distance = None
    return(distance, len(s1_act_all.union(s2_act_all)))


samples=['PS00718','PS00756','PS00757','PS00758','PS00867','PS00868','PS00869','PS00870','PS00871','PS00872','PS00873','PS00874']

out_rows_bin = []
out_rows_bin.append(['Cell1','Cell2','Jaccard_Dist','Group','Bin','n_items'])
for b in range(-1,12):
    for combo in combinations(samples, 2):
        dist,n = actuation_jaccard_dist_prox_bins(combo[0], combo[1], True, b)
        if dist != None:
            out_rows_bin.append([combo[0], combo[1], dist, 'Proximal', b, n])

with open('jaccard_acc_proximal_binned_by_expression.tsv','w') as fout:
    writer = csv.writer(fout, delimiter="\t")
    for row in out_rows_bin:
        writer.writerow(row)

