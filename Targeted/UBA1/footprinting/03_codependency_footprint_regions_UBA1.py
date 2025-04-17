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

# FILTER to regions within the canonical UBA1 promoter -------------------------------------------------------------------
# bedtools intersect -wa -a region_footprint_props_UBA1_H1.bed -b canonical_uba1_promoter.bed
# Regions 31-36

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



# Conditional codependency ---------------------------------------

# function for creating codep leave-one-out graph
def loo_codep_graph(omit_node, z_df, peaks_include):
    total_score = 0
    n_pairs = 0
    filt_peaks = [f'{i}' for i in peaks_include if i != omit_node]
    filt_df = z_df[z_df[f'{omit_node}'] == 0]
    for pair in tuple(combinations(filt_peaks, 2)):
        r1 = pair[0]
        r2 = pair[1]
        counts1 = Counter(filt_df[r1])
        counts2 = Counter(filt_df[r2])
        if ((counts1[0]+counts1[1]) > 0) and ((counts2[0]+counts2[1]) > 0):
            # prop bound calculated on MSPs that span BOTH regions
            n_both_bound = len(filt_df[(filt_df[r1] == 1) & (filt_df[r2] == 1)])
            n_both_msp = len(filt_df[(filt_df[r1].isin([0,1])) & (filt_df[r2].isin([0,1]))])
            pb1 = len(filt_df[(filt_df[r1] == 1) & (filt_df[r2].isin([0,1]))]) / n_both_msp
            pb2 = len(filt_df[(filt_df[r2] == 1) & (filt_df[r1].isin([0,1]))]) / n_both_msp
            expected = pb1 * pb2
            if n_both_msp > 0:
                observed = n_both_bound / n_both_msp
                total_score += observed - expected
                n_pairs += 1
                graph_edge_weights[omit_node][':'.join(pair)] = observed - expected
            else:
                graph_edge_weights[omit_node][':'.join(pair)] = 'NA' # NA is either no MSPs overlapping both or only footprints on opposite strands (5/11, GA/CT, both PATZ1)
    if n_pairs > 0:
        total_score /= n_pairs
    else:
        return(None)
    return(total_score)


# H1
# Canonical promoter
peaks_include = range(31,36+1)
# load baseline codep scores
baseline_total = 0
baseline_scores = dict()
for k,v in scores_H1.items():
    r1,r2 = k.split(':')
    if int(r1) in peaks_include and int(r2) in peaks_include:
        baseline_total += v
        baseline_scores[k] = v
# normalize total codep scores by number of comparisons
baseline_total /= len(baseline_scores)
# track edge weights from all graphs
graph_edge_weights = dict()
for i in peaks_include:
    graph_edge_weights[i] = {k:"NA" for k in baseline_scores.keys()}
# leave out each peak
out_lines = []
for i in peaks_include:
    total_score = loo_codep_graph(i, z_df_H1, peaks_include)
    if total_score != None:
        out_lines.append([i, baseline_total/total_score])
# write output
with open(f"conditional_codependency_graphs_raw_H1_{'-'.join([str(i) for i in peaks_include])}.tsv",'w') as fout:
    writer = csv.writer(fout, delimiter="\t")
    header = ['peak_omitted','codep_ratio']
    for k in baseline_scores.keys():
        header.append(k)
    writer.writerow(header)
    base_out = ['baseline','NA']
    for k in graph_edge_weights[peaks_include[0]].keys():
        if k in scores_H1.keys():
            base_out.append(scores_H1[k])
        # else:
        #     base_out.append('NA')
    writer.writerow(base_out)
    for row in out_lines:
        for k in graph_edge_weights[peaks_include[0]].keys():
            if k in scores_H1.keys():
                edge_weight = graph_edge_weights[row[0]][k]
                row.append(edge_weight)
            # else:
            #     base_out.append('NA')
        writer.writerow(row)


# H2
# Canonical promoter
peaks_include = range(31,36+1)
# load baseline codep scores
baseline_total = 0
baseline_scores = dict()
for k,v in scores_H2.items():
    r1,r2 = k.split(':')
    if int(r1) in peaks_include and int(r2) in peaks_include:
        baseline_total += v
        baseline_scores[k] = v
# normalize total codep scores by number of comparisons
baseline_total /= len(baseline_scores)
# track edge weights from all graphs
graph_edge_weights = dict()
for i in peaks_include:
    graph_edge_weights[i] = {k:"NA" for k in baseline_scores.keys()}
# leave out each peak
out_lines = []
for i in peaks_include:
    total_score = loo_codep_graph(i, z_df_H2, peaks_include)
    if total_score != None:
        out_lines.append([i, baseline_total/total_score])
# write output
with open(f"conditional_codependency_graphs_raw_H2_{'-'.join([str(i) for i in peaks_include])}.tsv",'w') as fout:
    writer = csv.writer(fout, delimiter="\t")
    header = ['peak_omitted','codep_ratio']
    for k in baseline_scores.keys():
        header.append(k)
    writer.writerow(header)
    base_out = ['baseline','NA']
    for k in graph_edge_weights[peaks_include[0]].keys():
        if k in scores_H2.keys():
            base_out.append(scores_H2[k])
        # else:
        #     base_out.append('NA')
    writer.writerow(base_out)
    for row in out_lines:
        for k in graph_edge_weights[peaks_include[0]].keys():
            if k in scores_H2.keys():
                edge_weight = graph_edge_weights[row[0]][k]
                row.append(edge_weight)
            # else:
            #     base_out.append('NA')
        writer.writerow(row)



# Paired t-test of UBA1 canonical promoter TF footprint codependency by haplotype ------------------------------------------------------------------------------------------------


# Hap 1 -----------------------------
codep_scores_H1 = pd.read_csv('TF_H1_Canonical_codep_scores_pairs.csv')
codep_scores_H2 = pd.read_csv('TF_H2_Canonical_codep_scores_pairs.csv')

#perform independent two sample t-test
ttest_rel(codep_scores_H1['score'], codep_scores_H2['score'])


