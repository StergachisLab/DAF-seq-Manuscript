import csv
import pandas as pd
import numpy as np
from itertools import combinations
from collections import Counter


# tail -n 1 ../napa_FIRE_peaks.bed > napa_RE2.bed
# bedtools intersect -a ../qc_plots/large_msp_positions.bed -b napa_RE2.bed > re2_large_msp.bed

# bedtools intersect -wa -F 1 -abam ../haplotype_correct/NAPA_PS00626_haplotype_corrected.bam -b napa_RE2.bed -bed | awk '{ print $4 }' > re2_fibers.txt 
spanning_fibers = set()
with open('re2_fibers.txt') as fr:
    for line in fr:
        spanning_fibers.add(line.strip())

msps = pd.read_csv('re2_large_msp.bed', sep="\t", names=['chr','st','en','zmw'])

# fibers with large MSP overlaping RE vs those without
z_df = pd.read_csv('zmw_footprint_regions.csv', index_col='zmw')
z_df = z_df[z_df.index.isin(spanning_fibers)]
regions = z_df.columns
z_df['RE2'] = 0
z_df.loc[z_df.index.isin(msps['zmw']), 'RE2'] = 1


# codep of each region with RE2
re2_counts = Counter(z_df['RE2'])
prop_re2 = re2_counts[1] / (re2_counts[0]+re2_counts[1])
scores = dict()
rows = []
for reg in regions:
    r1 = reg
    r2 = 'RE2'
    counts1 = Counter(z_df[r1])
    pb1 = counts1[1]/(counts1[0]+counts1[1])
    pb2 = prop_re2
    expected = pb1 * pb2
    both_act = len(z_df[(z_df[r1] == 1) & (z_df[r2] == 1)])
    both_spanned = len(z_df[z_df[r1].isin([0,1])])
    if both_spanned > 0:
        observed = both_act / both_spanned
        rows.append([r1, r2, pb1, pb2, expected, observed, (observed-expected)*4]) # scale from -1 to 1 !!!
        scores[':'.join([r1,r2])] = (observed-expected)*4

# No noticeable effect for any region !!!!!!!!!!!!!!!!!!!!!!!


# codep of 3,4,5 with RE2
re2_counts = Counter(z_df['RE2'])
prop_re2 = re2_counts[1] / (re2_counts[0]+re2_counts[1])
scores_all = dict()
rows_all = []
r1 = '3'
r2 = '4'
r3 = '5'
r4 = 'RE2'
pb1 = len(z_df[(z_df[r1] == 1) & (z_df[r2] == 1) & (z_df[r3] == 1)])/len(z_df[(z_df[r1].isin([0,1])) & (z_df[r2].isin([0,1])) & (z_df[r3].isin([0,1]))]) # bound at all 3 / MSP at all 3
pb2 = prop_re2
expected = pb1 * pb2
both_act = len(z_df[(z_df[r1] == 1) & (z_df[r2] == 1) & (z_df[r3] == 1) & (z_df[r4] == 1)])
both_spanned = len(z_df[(z_df[r1].isin([0,1])) & (z_df[r2].isin([0,1])) & (z_df[r3].isin([0,1]))])
if both_spanned > 0:
    observed = both_act / both_spanned
    rows_all.append([r1, r2, r3, r4, pb1, pb2, expected, observed, (observed-expected)*4]) # scale from -1 to 1 !!!
    scores_all[':'.join([r1,r2,r3,r4])] = (observed-expected)*4

# No noticeable effect for bound 3,4,5 with RE2 !!!!!!!!!!!!!!!!!!!!!!!


# +1 nucleosome, taking approx center of downstream nuc peak in pileup (47,514,930)
# bedtools intersect -a ../qc_plots/nuc_positions.bed -b plus_one_nuc.bed > plus_one_nuc_fibers.bed

# Fiber must overlap a range of +/- 146 bp from the estimated +1 Nuc center !!!!!!!!!!
# bedtools intersect -wa -F 1 -abam ../haplotype_correct/NAPA_PS00626_haplotype_corrected.bam -b plus_one_nuc_range.bed -bed | awk '{ print $4 }' > plus_one_nuc_fibers.txt

spanning_fibers_nucs = set()
with open('plus_one_nuc_fibers.txt') as fr:
    for line in fr:
        spanning_fibers_nucs.add(line.strip())

nucs = pd.read_csv('plus_one_nuc_fibers.bed', sep="\t", names=['chr','st','en','zmw'])

# fibers with large MSP overlaping RE vs those without
z_df = pd.read_csv('zmw_footprint_regions.csv', index_col='zmw')
z_df = z_df[z_df.index.isin(spanning_fibers_nucs)]
regions = z_df.columns
z_df['nuc'] = 0
z_df.loc[z_df.index.isin(nucs['zmw']), 'nuc'] = 1

# codep of each region with +1 Nuc
nuc_counts = Counter(z_df['nuc'])
scores_nuc = dict()
rows_nuc = []
for reg in regions:
    r1 = reg
    r2 = 'nuc'
    counts1 = Counter(z_df[r1])
    pb1 = counts1[1]/(counts1[0]+counts1[1])
    pb2 = len(z_df[(z_df[r1].isin([0,1])) & (z_df['nuc'] == 1)]) / len(z_df[(z_df[r1].isin([0,1])) & (z_df['nuc'].isin([0,1]))])
    expected = pb1 * pb2
    both_act = len(z_df[(z_df[r1] == 1) & (z_df[r2] == 1)])
    both_spanned = len(z_df[(z_df[r1].isin([0,1])) & (z_df[r2].isin([0,1]))])
    if both_spanned > 0:
        observed = both_act / both_spanned
        rows_nuc.append([r1, r2, pb1, pb2, expected, observed, (observed-expected)*4]) # scale from -1 to 1 !!!
        scores_nuc[':'.join([r1,r2])] = (observed-expected)*4

# no strong scores for any regions!!!!!!!!!!!!!!!!!!!!!!!!!!


# codep of 3,4,5 with +1 Nuc
enc = 1
nuc_counts = Counter(z_df['nuc'])
scores_nuc_all = dict()
rows_nuc_all = []
r1 = '3'
r2 = '4'
r3 = '5'
r4 = 'nuc'
counts1 = Counter(z_df[r1])
pb1 = len(z_df[(z_df[r1] == 1) & (z_df[r2] == 1) & (z_df[r3] == 1)])/len(z_df[(z_df[r1].isin([0,1])) & (z_df[r2].isin([0,1])) & (z_df[r3].isin([0,1]))]) # bound at all 3 / MSP at all 3
pb2 = len(z_df[(z_df[r1].isin([0,1])) & (z_df[r2].isin([0,1])) & (z_df[r3].isin([0,1])) & (z_df['nuc'] == 1)]) / len(z_df[(z_df[r1].isin([0,1])) & (z_df[r2].isin([0,1])) & (z_df[r3].isin([0,1])) & (z_df['nuc'].isin([0,1]))])
expected = pb1 * pb2
both_act = len(z_df[(z_df[r1] == enc) & (z_df[r2] == enc) & (z_df[r3] == enc) & (z_df[r4] == 1)])
both_spanned = len(z_df[(z_df[r1].isin([0,1])) & (z_df[r2].isin([0,1])) & (z_df[r3].isin([0,1])) & (z_df['nuc'].isin([0,1]))])
if both_spanned > 0:
    observed = both_act / both_spanned
    rows_nuc_all.append([r1, r2, r3, r4, pb1, pb2, expected, observed, (observed-expected)*4]) # scale from -1 to 1 !!!
    scores_nuc_all[':'.join([r1,r2,r3,r4])] = (observed-expected)*4

# {'3:4:5:nuc': -0.017970903611012723} not significant!!!!!!!!!!!!!!!!!!!!!


# BAMs of +/- 3,4,5 & +/- +1 Nuc
unbound_nuc = z_df[(z_df[r1] == 0) & (z_df[r2] == 0) & (z_df[r3] == 0) & (z_df[r4] == 1)]
unbound_NO_nuc = z_df[(z_df[r1] == 0) & (z_df[r2] == 0) & (z_df[r3] == 0) & (z_df[r4] == 0)]

# save ZMWs for visual inspection of +1 nucleosome
bound_nuc = z_df[(z_df[r1] == 1) & (z_df[r2] == 1) & (z_df[r3] == 1) & (z_df[r4] == 1)]
bound_NO_nuc = z_df[(z_df[r1] == 1) & (z_df[r2] == 1) & (z_df[r3] == 1) & (z_df[r4] == 0)]

with open('unbound_345_p1_nuc_zmws.txt','w') as fw:
    for z in unbound_nuc.index:
        fw.write(f'{z}\n')

with open('unbound_345_NO_p1_nuc_zmws.txt','w') as fw:
    for z in unbound_NO_nuc.index:
        fw.write(f'{z}\n')

with open('bound_345_p1_nuc_zmws.txt','w') as fw:
    for z in bound_nuc.index:
        fw.write(f'{z}\n')

with open('bound_345_NO_p1_nuc_zmws.txt','w') as fw:
    for z in bound_NO_nuc.index:
        fw.write(f'{z}\n')

# samtools view -b -N bound_345_p1_nuc_zmws.txt ../haplotype_correct/NAPA_PS00626_haplotype_corrected.bam > bound_345_p1_nuc_zmws.bam


# fibers with called YY2 footprint
yy2_bound = z_df[z_df['1'] == 1]
with open('yy2_bound_zmws.txt','w') as fw:
    for z in yy2_bound.index:
        fw.write(f'{z}\n')




# fibers with large MSP overlaping RE vs those without
z_df = pd.read_csv('zmw_footprint_regions.csv', index_col='zmw')
regions = z_df.columns


# codep of 3,4,5 with 1
enc = 1
scores_one = dict()
rows_one = []
r1 = '3'
r2 = '4'
r3 = '5'
r4 = '1'
counts1 = Counter(z_df[r1])
pb1 = len(z_df[(z_df[r1] == 1) & (z_df[r2] == 1) & (z_df[r3] == 1)])/len(z_df[(z_df[r1].isin([0,1])) & (z_df[r2].isin([0,1])) & (z_df[r3].isin([0,1]))]) # bound at all 3 / MSP at all 3
pb2 = len(z_df[(z_df[r1].isin([0,1])) & (z_df[r2].isin([0,1])) & (z_df[r3].isin([0,1])) & (z_df[r4] == 1)]) / len(z_df[(z_df[r1].isin([0,1])) & (z_df[r2].isin([0,1])) & (z_df[r3].isin([0,1])) & (z_df[r4].isin([0,1]))])
expected = pb1 * pb2
both_act = len(z_df[(z_df[r1] == enc) & (z_df[r2] == enc) & (z_df[r3] == enc) & (z_df[r4] == 1)])
both_spanned = len(z_df[(z_df[r1].isin([0,1])) & (z_df[r2].isin([0,1])) & (z_df[r3].isin([0,1])) & (z_df[r4].isin([0,1]))])
if both_spanned > 0:
    observed = both_act / both_spanned
    rows_one.append([r1, r2, r3, r4, pb1, pb2, expected, observed, (observed-expected)*4]) # scale from -1 to 1 !!!
    scores_one[':'.join([r1,r2,r3,r4])] = (observed-expected)*4

# not significant {'3:4:5:1': 0.031335695653056855}





