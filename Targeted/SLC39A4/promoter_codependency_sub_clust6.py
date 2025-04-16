import csv
import os
from glob import glob
import pandas as pd
from collections import Counter
import pysam
from itertools import combinations


# conditional co-dependency analysis with these sub-clusters separately on T and C allele reads from Liver
# to see if the haplotype alters the co-dependency of how this promoter is actuated


# samtools view -b -d HP:1 -N /mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/SLC39A4/clustering/clust_pileups/cluster_6_DA_GA_zmw.txt /mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/SLC39A4/pileups/Liver_SLC39A4_PS00680_nuc.bam > clust6_Liver_SLC39A4_PS00680_nuc_GA_H1.bam
# samtools index clust6_Liver_SLC39A4_PS00680_nuc_GA_H1.bam
# samtools view -b -d HP:2 -N /mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/SLC39A4/clustering/clust_pileups/cluster_6_DA_GA_zmw.txt /mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/SLC39A4/pileups/Liver_SLC39A4_PS00680_nuc.bam > clust6_Liver_SLC39A4_PS00680_nuc_GA_H2.bam
# samtools index clust6_Liver_SLC39A4_PS00680_nuc_GA_H2.bam

# /mmfs1/gscratch/stergachislab/bin/ft extract -x "len(msp)>150" --all clust6_Liver_SLC39A4_PS00680_ft_extract_ALL_H1.bed clust6_Liver_SLC39A4_PS00680_nuc_GA_H1.bam
# /mmfs1/gscratch/stergachislab/bin/ft extract -x "len(msp)>150" --all clust6_Liver_SLC39A4_PS00680_ft_extract_ALL_H2.bed clust6_Liver_SLC39A4_PS00680_nuc_GA_H2.bam


# extract_H1 = 'clust6_Liver_SLC39A4_PS00680_ft_extract_ALL_H1.bed'
# extract_H2 = 'clust6_Liver_SLC39A4_PS00680_ft_extract_ALL_H2.bed'

# name = 'clust6_Liver_SLC39A4_PS00680_H1'
# df_H1 = pd.read_csv(extract_H1, sep = "\t")
# small_msps = [] # <75bp
# large_msps = [] # >150 bp
# for i in range(len(df_H1)):
#     chrom = df_H1.iloc[i,:]['#ct']
#     msps = df_H1.iloc[i,:]['ref_msp_starts'].rstrip(',').split(',')
#     mspl = df_H1.iloc[i,:]['ref_msp_lengths'].rstrip(',').split(',')
#     if len(msps) > 0 and '.' not in msps:
#         for j in range(len(msps)):
#             if int(mspl[j]) < 75 and int(mspl[j]) > 0:
#                 small_msps.append([chrom, int(msps[j]), int(msps[j])+int(mspl[j])+1, df_H1.iloc[i,:]['fiber']])
#             elif int(mspl[j]) > 150:
#                 large_msps.append([chrom, int(msps[j]), int(msps[j])+int(mspl[j])+1, df_H1.iloc[i,:]['fiber']])
# large_msps = pd.DataFrame(large_msps, columns=['chrom','st','en','zmw']).sort_values(['st', 'en'], ascending=[True, True])
# large_bed = f'{name}_large_msp_positions.bed'
# large_msps.to_csv(large_bed, sep="\t", header=False, index=False)

# name = 'clust6_Liver_SLC39A4_PS00680_H2'
# df_H2 = pd.read_csv(extract_H2, sep = "\t")
# small_msps = [] # <75bp
# large_msps = [] # >150 bp
# for i in range(len(df_H2)):
#     chrom = df_H2.iloc[i,:]['#ct']
#     msps = df_H2.iloc[i,:]['ref_msp_starts'].rstrip(',').split(',')
#     mspl = df_H2.iloc[i,:]['ref_msp_lengths'].rstrip(',').split(',')
#     if len(msps) > 0 and '.' not in msps:
#         for j in range(len(msps)):
#             if int(mspl[j]) < 75 and int(mspl[j]) > 0:
#                 small_msps.append([chrom, int(msps[j]), int(msps[j])+int(mspl[j])+1, df_H2.iloc[i,:]['fiber']])
#             elif int(mspl[j]) > 150:
#                 large_msps.append([chrom, int(msps[j]), int(msps[j])+int(mspl[j])+1, df_H2.iloc[i,:]['fiber']])
# large_msps = pd.DataFrame(large_msps, columns=['chrom','st','en','zmw']).sort_values(['st', 'en'], ascending=[True, True])
# large_bed = f'{name}_large_msp_positions.bed'
# large_msps.to_csv(large_bed, sep="\t", header=False, index=False)


# MAP MSPs to modules ----------------------------------------------
# bedmap --ec --fraction-ref 0.5 --echo --echo-map-id promoter_modules_sub_clust6.bed <( sort-bed clust6_Liver_SLC39A4_PS00680_H1_large_msp_positions.bed ) > clust6_large_msps_H1.txt
# bedmap --ec --fraction-ref 0.5 --echo --echo-map-id promoter_modules_sub_clust6.bed <( sort-bed clust6_Liver_SLC39A4_PS00680_H2_large_msp_positions.bed ) > clust6_large_msps_H2.txt


# Hap 1 -------------------------------------------------------------------------------------------------------------------------------
peaks = dict()
peak_regions = dict()
bmap='clust6_large_msps_H1.txt'
with open(bmap) as fr:
    i = 1
    for line in fr:
        bed,zmw = line.strip().split('|')
        peaks[i] = zmw.split(';')
        bed2 = bed.split('\t')
        peak_regions[i] = f"{bed2[0]}:{bed2[1]}-{bed2[2]}"
        i+=1

# total fiber number
fibers = set()
bam_name = 'clust6_Liver_SLC39A4_PS00680_nuc_GA_H1.bam'
bam = pysam.AlignmentFile(bam_name, "rb")
for read in bam.fetch():
    fibers.add(read.qname)
# DF of zmws by region footprint status (not overlapped by fiber (should be filtered out) -> -1, no MSP -> 0, MSP -> 1)
msp_df_H1 = pd.DataFrame(fibers, columns=['zmw'])
for k in peaks.keys():
    new_vals = []
    for i in range(len(msp_df_H1)):
        z = msp_df_H1.iloc[i,:]['zmw']
        if z in peaks[k]:
            new_vals.append(1)
        else:
            new_vals.append(0)
    msp_df_H1[k] = new_vals
msp_df_H1.index = msp_df_H1['zmw']
msp_df_H1 = msp_df_H1.drop(['zmw'], axis=1)
msp_df_H1.to_csv('clust6_msp_df_Liver_H1.csv')

# calculate codependency for each pair of elements
scores_H1 = dict()
rows = []
for combo in combinations(msp_df_H1.columns, 2):
    r1 = combo[0]
    r2 = combo[1]
    counts1 = Counter(msp_df_H1[r1])
    counts2 = Counter(msp_df_H1[r2])
    if ((counts1[0]+counts1[1]) > 0) and ((counts2[0]+counts2[1]) > 0):
        n_both_msp = len(msp_df_H1[(msp_df_H1[r1] == 1) & (msp_df_H1[r2] == 1)])
        n_both_span = len(msp_df_H1[(msp_df_H1[r1].isin([0,1])) & (msp_df_H1[r2].isin([0,1]))])
        pb1 = len(msp_df_H1[(msp_df_H1[r1] == 1) & (msp_df_H1[r2].isin([0,1]))]) / n_both_span
        pb2 = len(msp_df_H1[(msp_df_H1[r2] == 1) & (msp_df_H1[r1].isin([0,1]))]) / n_both_span
        expected = pb1 * pb2
        if n_both_span > 0:
            observed = n_both_msp / n_both_span
            rows.append([r1, r2, pb1, pb2, expected, observed, (observed-expected)*4]) # scale from -1 to 1 !!!
            scores_H1[f'{r1}:{r2}'] = observed-expected
with open('codep_scores_pairs_clust6_Liver_H1.csv','w') as fw:
    writer = csv.writer(fw)
    writer.writerow(['reg1','reg2','prop1','prop2','expected','observed','score'])
    for r in rows:
        writer.writerow(r)


# Hap 2 -------------------------------------------------------------------------------------------------------------------------------
peaks = dict()
peak_regions = dict()
bmap='clust6_large_msps_H2.txt'
with open(bmap) as fr:
    i = 1
    for line in fr:
        bed,zmw = line.strip().split('|')
        peaks[i] = zmw.split(';')
        bed2 = bed.split('\t')
        peak_regions[i] = f"{bed2[0]}:{bed2[1]}-{bed2[2]}"
        i+=1

# total fiber number
fibers = set()
bam_name = 'clust6_Liver_SLC39A4_PS00680_nuc_GA_H2.bam'
bam = pysam.AlignmentFile(bam_name, "rb")
for read in bam.fetch():
    fibers.add(read.qname)
# DF of zmws by region footprint status (not overlapped by fiber (should be filtered out) -> -1, no MSP -> 0, MSP -> 1)
msp_df_H2 = pd.DataFrame(fibers, columns=['zmw'])
for k in peaks.keys():
    new_vals = []
    for i in range(len(msp_df_H2)):
        z = msp_df_H2.iloc[i,:]['zmw']
        if z in peaks[k]:
            new_vals.append(1)
        else:
            new_vals.append(0)
    msp_df_H2[k] = new_vals
msp_df_H2.index = msp_df_H2['zmw']
msp_df_H2 = msp_df_H2.drop(['zmw'], axis=1)
msp_df_H2.to_csv('clust6_msp_df_Liver_H2.csv')

# calculate codependency for each pair of elements
scores_H2 = dict()
rows = []
for combo in combinations(msp_df_H2.columns, 2):
    r1 = combo[0]
    r2 = combo[1]
    counts1 = Counter(msp_df_H2[r1])
    counts2 = Counter(msp_df_H2[r2])
    if ((counts1[0]+counts1[1]) > 0) and ((counts2[0]+counts2[1]) > 0):
        n_both_msp = len(msp_df_H2[(msp_df_H2[r1] == 1) & (msp_df_H2[r2] == 1)])
        n_both_span = len(msp_df_H2[(msp_df_H2[r1].isin([0,1])) & (msp_df_H2[r2].isin([0,1]))])
        pb1 = len(msp_df_H2[(msp_df_H2[r1] == 1) & (msp_df_H2[r2].isin([0,1]))]) / n_both_span
        pb2 = len(msp_df_H2[(msp_df_H2[r2] == 1) & (msp_df_H2[r1].isin([0,1]))]) / n_both_span
        expected = pb1 * pb2
        if n_both_span > 0:
            observed = n_both_msp / n_both_span
            rows.append([r1, r2, pb1, pb2, expected, observed, (observed-expected)*4]) # scale from -1 to 1 !!!
            scores_H2[f'{r1}:{r2}'] = observed-expected
with open('codep_scores_pairs_clust6_Liver_H2.csv','w') as fw:
    writer = csv.writer(fw)
    writer.writerow(['reg1','reg2','prop1','prop2','expected','observed','score'])
    for r in rows:
        writer.writerow(r)



# Conditional codependency ---------------------------------------

# function for creating codep leave-one-out graph
def loo_codep_graph(omit_node, z_df, peaks_include):
    total_score = 0
    n_pairs = 0
    filt_peaks = [i for i in peaks_include if i != omit_node]
    filt_df = z_df[z_df[omit_node] == 0]
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
                graph_edge_weights[omit_node][f'{pair[0]}:{pair[1]}'] = observed - expected
            else:
                graph_edge_weights[omit_node][f'{pair[0]}:{pair[1]}'] = 'NA'
    if n_pairs > 0:
        total_score /= n_pairs
    else:
        return(None)
    return(total_score)


# H1
# Canonical promoter
peaks_include = range(1,6)
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
    total_score = loo_codep_graph(i, msp_df_H1, peaks_include)
    if total_score != None:
        out_lines.append([i, baseline_total/total_score])
# write output
with open(f"conditional_codependency_graphs_raw_clust6_Liver_H1_{'-'.join([str(i) for i in peaks_include])}.tsv",'w') as fout:
    writer = csv.writer(fout, delimiter="\t")
    header = ['peak_omitted','codep_ratio']
    for k in baseline_scores.keys():
        header.append(k)
    writer.writerow(header)
    base_out = ['baseline','NA']
    for k in graph_edge_weights[peaks_include[0]].keys():
        if k in scores_H1.keys():
            base_out.append(scores_H1[k])
    writer.writerow(base_out)
    for row in out_lines:
        for k in graph_edge_weights[peaks_include[0]].keys():
            if k in scores_H1.keys():
                edge_weight = graph_edge_weights[row[0]][k]
                row.append(edge_weight)
        writer.writerow(row)


# -----------------------------------------------------------------------------------------------------------

# H2
# Canonical promoter
peaks_include = range(1,6)
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
    total_score = loo_codep_graph(i, msp_df_H2, peaks_include)
    if total_score != None:
        out_lines.append([i, baseline_total/total_score])
# write output
with open(f"conditional_codependency_graphs_raw_clust6_Liver_H2_{'-'.join([str(i) for i in peaks_include])}.tsv",'w') as fout:
    writer = csv.writer(fout, delimiter="\t")
    header = ['peak_omitted','codep_ratio']
    for k in baseline_scores.keys():
        header.append(k)
    writer.writerow(header)
    base_out = ['baseline','NA']
    for k in graph_edge_weights[peaks_include[0]].keys():
        if k in scores_H2.keys():
            base_out.append(scores_H2[k])
    writer.writerow(base_out)
    for row in out_lines:
        for k in graph_edge_weights[peaks_include[0]].keys():
            if k in scores_H2.keys():
                edge_weight = graph_edge_weights[row[0]][k]
                row.append(edge_weight)
        writer.writerow(row)




