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

fire_df.to_csv('scDAF_FIRE_actuation_MSP150_withProximal.tsv', sep="\t", index=False)


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

# Jaccard distance between OPPOSITE haplotypes
def actuation_jaccard_dist_overall_opposite_hap(s1, s2):
    # S1H1 X S2H2
    h12_df = fire_df[(fire_df[f'{s1}_H1'] >= 0) & (fire_df[f'{s2}_H2'] >= 0)]
    s1_act_H12 = set([f'H12_{idx}' for idx in h12_df[h12_df[f'{s1}_H1'] > 0].index])
    s2_act_H12 = set([f'H12_{idx}' for idx in h12_df[h12_df[f'{s2}_H2'] > 0].index])
    # S1H2 X S2H1
    h21_df = fire_df[(fire_df[f'{s1}_H2'] >= 0) & (fire_df[f'{s2}_H1'] >= 0)]
    s1_act_H21 = set([f'H21_{idx}' for idx in h21_df[h21_df[f'{s1}_H2'] > 0].index])
    s2_act_H21 = set([f'H21_{idx}' for idx in h21_df[h21_df[f'{s2}_H1'] > 0].index])
    # merge datapoint (both haplotypes)
    s12_act_all = s1_act_H12.union(s1_act_H21)
    s21_act_all = s2_act_H12.union(s2_act_H21)
    distance = 1 - jaccard_similarity(s12_act_all, s21_act_all)
    return(distance)

# Jaccard distance between OPPOSITE haplotypes
def actuation_jaccard_dist_prox_opposite_hap(s1, s2, proximal):
    # S1H1 X S2H2
    h12_df = fire_df[(fire_df[f'{s1}_H1'] >= 0) & (fire_df[f'{s2}_H2'] >= 0) & (fire_df['is_proximal'] == proximal)]
    s1_act_H12 = set([f'H12_{idx}' for idx in h12_df[h12_df[f'{s1}_H1'] > 0].index])
    s2_act_H12 = set([f'H12_{idx}' for idx in h12_df[h12_df[f'{s2}_H2'] > 0].index])
    # S1H2 X S2H1
    h21_df = fire_df[(fire_df[f'{s1}_H2'] >= 0) & (fire_df[f'{s2}_H1'] >= 0) & (fire_df['is_proximal'] == proximal)]
    s1_act_H21 = set([f'H21_{idx}' for idx in h21_df[h21_df[f'{s1}_H2'] > 0].index])
    s2_act_H21 = set([f'H21_{idx}' for idx in h21_df[h21_df[f'{s2}_H1'] > 0].index])
    # merge datapoint (both haplotypes)
    s1_act_all = s1_act_H12.union(s1_act_H21)
    s2_act_all = s2_act_H12.union(s2_act_H21)
    distance = 1 - jaccard_similarity(s1_act_all, s2_act_all)
    return(distance)

def actuation_jaccard_dist_Within_Cell(s):
    df = fire_df[(fire_df[f'{s}_H1'] >= 0) & (fire_df[f'{s}_H2'] >= 0)]
    s_act_H1 = set([f'{idx}' for idx in df[df[f'{s}_H1'] > 0].index])
    s_act_H2 = set([f'{idx}' for idx in df[df[f'{s}_H2'] > 0].index])
    distance = 1 - jaccard_similarity(s_act_H1, s_act_H2)
    return(distance)

def actuation_jaccard_dist_Within_Cell_prox(s, proximal):
    df = fire_df[(fire_df[f'{s}_H1'] >= 0) & (fire_df[f'{s}_H2'] >= 0) & (fire_df['is_proximal'] == proximal)]
    s_act_H1 = set([f'{idx}' for idx in df[df[f'{s}_H1'] > 0].index])
    s_act_H2 = set([f'{idx}' for idx in df[df[f'{s}_H2'] > 0].index])
    distance = 1 - jaccard_similarity(s_act_H1, s_act_H2)
    return(distance)

def num_peaks(s1, s2):
    h1_df = fire_df[(fire_df[f'{s1}_H1'] >= 0) & (fire_df[f'{s2}_H1'] >= 0)]
    h2_df = fire_df[(fire_df[f'{s1}_H2'] >= 0) & (fire_df[f'{s2}_H2'] >= 0)]
    # merge datapoint (both haplotypes)
    return(len(h1_df) + len(h2_df))


samples=['PS00718','PS00756','PS00757','PS00758','PS00867','PS00868','PS00869','PS00870','PS00871','PS00872','PS00873','PS00874']
out_rows = []
out_rows.append(['Cell1','Cell2','Jaccard_Dist','Group'])

for combo in combinations(samples, 2):
    dist = actuation_jaccard_dist_prox(combo[0], combo[1], True)
    out_rows.append([combo[0], combo[1], dist, 'Proximal'])
    dist = actuation_jaccard_dist_prox(combo[0], combo[1], False)
    out_rows.append([combo[0], combo[1], dist, 'Distal'])
    dist = actuation_jaccard_dist_overall(combo[0], combo[1])
    out_rows.append([combo[0], combo[1], dist, 'Overall'])
    # Opposite Hap
    dist = actuation_jaccard_dist_prox_opposite_hap(combo[0], combo[1], True)
    out_rows.append([combo[0], combo[1], dist, 'Proximal_Opposite_Hap'])
    dist = actuation_jaccard_dist_prox_opposite_hap(combo[0], combo[1], False)
    out_rows.append([combo[0], combo[1], dist, 'Distal_Opposite_Hap'])
    dist = actuation_jaccard_dist_overall_opposite_hap(combo[0], combo[1])
    out_rows.append([combo[0], combo[1], dist, 'Overall_Opposite_Hap'])

# Diff haps WITHIN Cell
for samp in samples:
    dist = actuation_jaccard_dist_Within_Cell(samp)
    out_rows.append([samp, samp, dist, 'Overall_Within_Cell'])
    dist = actuation_jaccard_dist_Within_Cell_prox(samp, True)
    out_rows.append([samp, samp, dist, 'Proximal_Within_Cell'])
    dist = actuation_jaccard_dist_Within_Cell_prox(samp, False)
    out_rows.append([samp, samp, dist, 'Distal_Within_Cell'])

with open('jaccard_distances_actuation_by_cell.tsv','w') as fout:
    writer = csv.writer(fout, delimiter="\t")
    for row in out_rows:
        writer.writerow(row)


# % actuation by Cell --> proximal vs distal ----------------------------------------------------
out_act = []
out_act.append(['Cell','is_proximal','prop'])

for samp in samples:
    temp_df = fire_df[['chrom','start','end','prop_acc','bin',f'{samp}_H1',f'{samp}_H2','is_proximal']]
    temp_df = temp_df.rename(columns={f'{samp}_H1':'H1', f'{samp}_H2':'H2'})
    temp_df = temp_df[(temp_df['H1'] >=0) & (temp_df['H2']>=0)]
    grouped = temp_df.groupby('is_proximal')[['H1','H2']]
    sum_df = grouped.sum().sum(axis=1)
    count_df = grouped.count().sum(axis=1)
    prop_act = sum_df.divide(count_df.values)
    out_act.append([samp, 'False', prop_act.loc[False]])
    out_act.append([samp, 'True', prop_act.loc[True]])

with open('prop_actuation_proximal_vs_distal_by_cell.tsv','w') as fout:
    writer = csv.writer(fout, delimiter="\t")
    for row in out_act:
        writer.writerow(row)


# percent FIRE peaks actuated on at least 1 hap within a cell ----------------------------------------------------
cov_only = 0
cell_act = 0
for samp in samples:
    temp_df = fire_df[['chrom','start','end','prop_acc','bin',f'{samp}_H1',f'{samp}_H2','is_proximal']]
    temp_df['any_hap_act'] = temp_df[[f'{samp}_H1',f'{samp}_H2']].max(axis=1)
    counts = temp_df['any_hap_act'].value_counts()
    cov_only += counts[0]
    cell_act += counts[1]

with open('percent_peaks_actuated_per_cell_anyHap.txt','w') as fw:
    fw.write(f'Inactive,Actuated,Proportion\n')
    fw.write(f'{cov_only},{cell_act},{cell_act/(cell_act+cov_only)}\n')
