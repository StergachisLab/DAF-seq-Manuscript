import csv
import pandas as pd
import numpy as np
from itertools import combinations
from collections import Counter
from scipy.stats import ttest_rel

stats = pd.read_csv('footprint_stats.csv')

# Hap1 ONLY
z_df_H1 = pd.read_csv('zmw_footprint_regions_UBA1_H1.csv', index_col='zmw')

# calculate codependency for each pair of elements
scores_H1 = dict()
rows_H1 = []
for combo in combinations(z_df_H1.columns, 2):
    r1 = combo[0]
    r2 = combo[1]
    counts1 = Counter(z_df_H1[r1])
    counts2 = Counter(z_df_H1[r2])
    if ((counts1[0]+counts1[1]) > 0) and ((counts2[0]+counts2[1]) > 0):
        # prop bound calculated on MSPs that span BOTH regions
        n_both_bound = len(z_df_H1[(z_df_H1[r1] == 1) & (z_df_H1[r2] == 1)])
        n_both_msp = len(z_df_H1[(z_df_H1[r1].isin([0,1])) & (z_df_H1[r2].isin([0,1]))])
        pb1 = len(z_df_H1[(z_df_H1[r1] == 1) & (z_df_H1[r2].isin([0,1]))]) / n_both_msp
        pb2 = len(z_df_H1[(z_df_H1[r2] == 1) & (z_df_H1[r1].isin([0,1]))]) / n_both_msp
        expected = pb1 * pb2
        if n_both_msp > 0:
            observed = n_both_bound / n_both_msp
            rows_H1.append([r1, r2, pb1, pb2, expected, observed, (observed-expected)*4]) # scale from -1 to 1 !!!
            scores_H1[':'.join([r1,r2])] = observed-expected

with open('TF_H1_codep_scores_pairs.csv','w') as fw:
    writer = csv.writer(fw)
    writer.writerow(['reg1','reg2','prop1','prop2','expected','observed','score'])
    for r in rows_H1:
        writer.writerow(r)


with open('TF_H1_Canonical_codep_scores_pairs.csv','w') as fw:
    writer = csv.writer(fw)
    writer.writerow(['reg1','reg2','prop1','prop2','expected','observed','score'])
    for r in rows_H1:
        if int(r[0]) in range(31,37) and int(r[1]) in range(31,37):
            writer.writerow(r)

# Hap2 ONLY
z_df_H2 = pd.read_csv('zmw_footprint_regions_UBA1_H2.csv', index_col='zmw')
scores_H2 = dict()
rows_H2 = []
for combo in combinations(range(31,37), 2):
    r1 = str(combo[0])
    r2 = str(combo[1])
    counts1 = Counter(z_df_H2[r1])
    counts2 = Counter(z_df_H2[r2])
    if ((counts1[0]+counts1[1]) > 0) and ((counts2[0]+counts2[1]) > 0):
        # prop bound calculated on MSPs that span BOTH regions
        n_both_bound = len(z_df_H2[(z_df_H2[r1] == 1) & (z_df_H2[r2] == 1)])
        n_both_msp = len(z_df_H2[(z_df_H2[r1].isin([0,1])) & (z_df_H2[r2].isin([0,1]))])
        pb1 = len(z_df_H2[(z_df_H2[r1] == 1) & (z_df_H2[r2].isin([0,1]))]) / n_both_msp
        pb2 = len(z_df_H2[(z_df_H2[r2] == 1) & (z_df_H2[r1].isin([0,1]))]) / n_both_msp
        expected = pb1 * pb2
        if n_both_msp > 0:
            observed = n_both_bound / n_both_msp
            rows_H2.append([r1, r2, pb1, pb2, expected, observed, (observed-expected)*4]) # scale from -1 to 1 !!!
            scores_H2[':'.join([r1,r2])] = observed-expected

with open('TF_H2_Canonical_codep_scores_pairs.csv','w') as fw:
    writer = csv.writer(fw)
    writer.writerow(['reg1','reg2','prop1','prop2','expected','observed','score'])
    for r in rows_H2:
        writer.writerow(r)


# Paired t-test of UBA1 canonical promoter TF footprint codependency by haplotype ------------------------------------------------------------------------------------------------

# Hap 1 -----------------------------
codep_scores_H1 = pd.read_csv('TF_H1_Canonical_codep_scores_pairs.csv')
codep_scores_H2 = pd.read_csv('TF_H2_Canonical_codep_scores_pairs.csv')

#perform independent two sample t-test
ttest_rel(codep_scores_H1['score'], codep_scores_H2['score'])
