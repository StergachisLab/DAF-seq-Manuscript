import csv
import pandas as pd
import numpy as np
from itertools import combinations
from collections import Counter


z_df = pd.read_csv('zmw_footprint_regions.csv', index_col='zmw')

# calculate codependency for each pair of elements
scores = dict()
rows = []
for combo in combinations(z_df.columns, 2):
    r1 = combo[0]
    r2 = combo[1]
    counts1 = Counter(z_df[r1])
    counts2 = Counter(z_df[r2])
    if ((counts1[0]+counts1[1]) > 0) and ((counts2[0]+counts2[1]) > 0):
        # prop bound calculated on MSPs that span BOTH regions
        n_both_bound = len(z_df[(z_df[r1] == 1) & (z_df[r2] == 1)])
        n_both_msp = len(z_df[(z_df[r1].isin([0,1])) & (z_df[r2].isin([0,1]))])
        pb1 = len(z_df[(z_df[r1] == 1) & (z_df[r2].isin([0,1]))]) / n_both_msp
        pb2 = len(z_df[(z_df[r2] == 1) & (z_df[r1].isin([0,1]))]) / n_both_msp
        expected = pb1 * pb2
    if n_both_msp > 0:
        observed = n_both_bound / n_both_msp
        rows.append([r1, r2, pb1, pb2, expected, observed, (observed-expected)*4]) # scale from -1 to 1 !!!
        scores[':'.join([r1,r2])] = observed-expected


with open('codep_scores_pairs.csv','w') as fw:
    writer = csv.writer(fw)
    writer.writerow(['reg1','reg2','prop1','prop2','expected','observed','score'])
    for r in rows:
        writer.writerow(r)


# Regions ---------------------------------------------------
# {1: ['chr19', '47514936', '47514946'],
#  2: ['chr19', '47514969', '47514985'],
#  3: ['chr19', '47515074', '47515096'],
#  4: ['chr19', '47515095', '47515112'],
#  5: ['chr19', '47515134', '47515152'],
#  6: ['chr19', '47515207', '47515217'],
#  7: ['chr19', '47515268', '47515279'],
#  8: ['chr19', '47515310', '47515324'],
#  9: ['chr19', '47515345', '47515363'],
#  10: ['chr19', '47515383', '47515396'],
#  11: ['chr19', '47515439', '47515456'],
#  12: ['chr19', '47515501', '47515511'],
#  13: ['chr19', '47515578', '47515612']}


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


# 1,2,3 & ALL
graphs = [[1,3], [1,11]]
for g in graphs:
    peaks_include = range(g[0],g[1]+1)
    # load baseline codep scores
    baseline_total = 0
    baseline_scores = dict()
    for k,v in scores.items():
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
        total_score = loo_codep_graph(i, z_df, peaks_include)
        out_lines.append([i, baseline_total/total_score])
    # write output
    with open(f"conditional_codependency_graphs_raw_{'-'.join([str(i) for i in peaks_include])}.tsv",'w') as fout:
        writer = csv.writer(fout, delimiter="\t")
        header = ['peak_omitted','codep_ratio']
        for k in baseline_scores.keys():
            header.append(k)
        writer.writerow(header)
        base_out = ['baseline','NA']
        for k in graph_edge_weights[peaks_include[0]].keys():
            if k in scores.keys():
                base_out.append(scores[k])
        writer.writerow(base_out)
        for row in out_lines:
            for k in graph_edge_weights[peaks_include[0]].keys():
                if k in scores.keys():
                    edge_weight = graph_edge_weights[row[0]][k]
                    row.append(edge_weight)
            writer.writerow(row)

