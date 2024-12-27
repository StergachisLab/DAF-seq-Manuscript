import pandas as pd
import numpy as np
from itertools import combinations
import csv

fire_df = pd.read_csv('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/msp_analysis/scDAF_FIRE_actuation_MSP150.tsv', sep="\t")

fire_df['ID'] = fire_df['chrom'].astype(str) + ":" + fire_df['start'].astype(str) + "-" + fire_df['end'].astype(str)

# split FIRE peaks into promoter-proximal and promoter-distal
# assignments made previously using "classify_FIRE_peaks_promoter_distal.py"
proximal_peaks = set()
with open('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/msp_analysis/promoter-proximal_FIRE-peaks_HG002_noUnrealiableCov_autosome.txt') as fr:
    for line in fr:
        proximal_peaks.add(line.strip())

def assign_proximal(row):
    if row['ID'] in proximal_peaks:
        return(True)
    else:
        return(False)

fire_df['is_proximal'] = fire_df.apply(assign_proximal, axis=1)

fire_df.to_csv('scDAF_FIRE_actuation_MSP150_withProximal.tsv', sep="\t")

df1 = fire_df[['chrom','start','end','prop_acc','bin','PS00718_H1','PS00718_H2','is_proximal']]
df1 = df1.rename(columns={'PS00718_H1':'H1', 'PS00718_H2':'H2'})
df2 = fire_df[['chrom','start','end','prop_acc','bin','PS00756_H1','PS00756_H2','is_proximal']]
df2 = df2.rename(columns={'PS00756_H1':'H1', 'PS00756_H2':'H2'})
df3 = fire_df[['chrom','start','end','prop_acc','bin','PS00757_H1','PS00757_H2','is_proximal']]
df3 = df3.rename(columns={'PS00757_H1':'H1', 'PS00757_H2':'H2'})
df4 = fire_df[['chrom','start','end','prop_acc','bin','PS00758_H1','PS00758_H2','is_proximal']]
df4 = df4.rename(columns={'PS00758_H1':'H1', 'PS00758_H2':'H2'})

# Jaccard distance between each pair of cells
def jaccard_similarity(set1, set2):
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return(intersection/union)

def actuation_jaccard_dist_overall(s1, s2):
    h1_df = fire_df[(fire_df[f'{s1}_H1'] >= 0) & (fire_df[f'{s2}_H1'] >= 0)]
    s1_act_H1 = set([f'H1_{idx}' for idx in h1_df[h1_df[f'{s1}_H1'] > 0].index])
    s2_act_H1 = set([f'H1_{idx}' for idx in h1_df[h1_df[f'{s2}_H1'] > 0].index])
    h2_df = fire_df[(fire_df[f'{s1}_H2'] >= 0) & (fire_df[f'{s2}_H2'] >= 0)]
    s1_act_H2 = set([f'H2_{idx}' for idx in h2_df[h2_df[f'{s1}_H2'] > 0].index])
    s2_act_H2 = set([f'H2_{idx}' for idx in h2_df[h2_df[f'{s2}_H2'] > 0].index])
    # merge datapoint (both haplotypes)
    s1_act_all = s1_act_H1.union(s1_act_H2)
    s2_act_all = s2_act_H1.union(s2_act_H2)
    distance = 1 - jaccard_similarity(s1_act_all, s2_act_all)
    return(distance)

def actuation_jaccard_dist_prox(s1, s2, proximal):
    h1_df = fire_df[(fire_df[f'{s1}_H1'] >= 0) & (fire_df[f'{s2}_H1'] >= 0) & (fire_df['is_proximal'] == proximal)]
    s1_act_H1 = set([f'H1_{idx}' for idx in h1_df[h1_df[f'{s1}_H1'] > 0].index])
    s2_act_H1 = set([f'H1_{idx}' for idx in h1_df[h1_df[f'{s2}_H1'] > 0].index])
    h2_df = fire_df[(fire_df[f'{s1}_H2'] >= 0) & (fire_df[f'{s2}_H2'] >= 0) & (fire_df['is_proximal'] == proximal)]
    s1_act_H2 = set([f'H2_{idx}' for idx in h2_df[h2_df[f'{s1}_H2'] > 0].index])
    s2_act_H2 = set([f'H2_{idx}' for idx in h2_df[h2_df[f'{s2}_H2'] > 0].index])
    # merge datapoint (both haplotypes)
    s1_act_all = s1_act_H1.union(s1_act_H2)
    s2_act_all = s2_act_H1.union(s2_act_H2)
    distance = 1 - jaccard_similarity(s1_act_all, s2_act_all)
    return(distance)

def num_peaks(s1, s2):
    h1_df = fire_df[(fire_df[f'{s1}_H1'] >= 0) & (fire_df[f'{s2}_H1'] >= 0)]
    h2_df = fire_df[(fire_df[f'{s1}_H2'] >= 0) & (fire_df[f'{s2}_H2'] >= 0)]
    # merge datapoint (both haplotypes)
    return(len(h1_df) + len(h2_df))


samples = ['PS00718','PS00756','PS00757','PS00758']

out_rows = []
out_rows.append(['Cell1','Cell2','Jaccard_Dist','Proximal'])

for combo in combinations(samples, 2):
    dist = actuation_jaccard_dist_prox(combo[0], combo[1], True)
    out_rows.append([combo[0], combo[1], dist, True])
    dist = actuation_jaccard_dist_prox(combo[0], combo[1], False)
    out_rows.append([combo[0], combo[1], dist, False])
    dist = actuation_jaccard_dist_overall(combo[0], combo[1])
    out_rows.append([combo[0], combo[1], dist, 'Overall'])
    print(combo, num_peaks(combo[0], combo[1]))


with open('jaccard_distances_actuation_by_cell.tsv','w') as fout:
    writer = csv.writer(fout, delimiter="\t")
    for row in out_rows:
        writer.writerow(row)


# % actuation by Cell --> proximal vs distal
df1 = df1[(df1['H1'] >=0) & (df1['H2']>=0)]
df2 = df2[(df2['H1'] >=0) & (df2['H2']>=0)]
df3 = df3[(df3['H1'] >=0) & (df3['H2']>=0)]
df4 = df4[(df4['H1'] >=0) & (df4['H2']>=0)]

df1['Cell'] = 'PS00718'
df2['Cell'] = 'PS00756'
df3['Cell'] = 'PS00757'
df4['Cell'] = 'PS00758'

out_act = []
out_act.append(['Cell','is_proximal','prop'])

# prop actuation by group & by cell
grouped = df1.groupby('is_proximal')[['H1','H2']]
sum_df = grouped.sum().sum(axis=1)
count_df = grouped.count().sum(axis=1)
prop_act = sum_df.divide(count_df.values)
out_act.append(['PS00718', 'False', prop_act.loc[False]])
out_act.append(['PS00718', 'True', prop_act.loc[True]])

grouped = df2.groupby('is_proximal')[['H1','H2']]
sum_df = grouped.sum().sum(axis=1)
count_df = grouped.count().sum(axis=1)
prop_act = sum_df.divide(count_df.values)
out_act.append(['PS00756', 'False', prop_act.loc[False]])
out_act.append(['PS00756', 'True', prop_act.loc[True]])

grouped = df3.groupby('is_proximal')[['H1','H2']]
sum_df = grouped.sum().sum(axis=1)
count_df = grouped.count().sum(axis=1)
prop_act = sum_df.divide(count_df.values)
out_act.append(['PS00757', 'False', prop_act.loc[False]])
out_act.append(['PS00757', 'True', prop_act.loc[True]])

grouped = df4.groupby('is_proximal')[['H1','H2']]
sum_df = grouped.sum().sum(axis=1)
count_df = grouped.count().sum(axis=1)
prop_act = sum_df.divide(count_df.values)
out_act.append(['PS00758', 'False', prop_act.loc[False]])
out_act.append(['PS00758', 'True', prop_act.loc[True]])

with open('prop_actuation_proximal_vs_distal_by_cell.tsv','w') as fout:
    writer = csv.writer(fout, delimiter="\t")
    for row in out_act:
        writer.writerow(row)

