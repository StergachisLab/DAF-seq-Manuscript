import csv
import os
from glob import glob
import pandas as pd
from collections import Counter
import pysam
from itertools import combinations

# TSSs as defined by Mas Iso-seq (first 4 of 5 are within FIRE peaks)
# bedtools intersect -a /gscratch/stergachislab/swansoe/projects/GM12878_scATAC/phased-fdr-and-peaks/XCI/isoseq/GM12878_mas_collapsed_TSS.bed -b region_UBA1.bed | head -n 4 > tss_uba1_masSeq.bed
# bedmap --ec --echo --echo-map-id tss_uba1_masSeq.bed <( sort-bed pileups/GM12878_UBA1_PS00685_H1_large_msp_positions.bed ) > large_msps_H1_fire_map.txt

# ALL 7 PEAKS
# bedmap --ec --fraction-ref 0.8 --echo --echo-map-id FIRE_narrowPeaks_GM12878.bed <( sort-bed pileups/GM12878_UBA1_PS00685_H1_large_msp_positions.bed ) > large_msps_H1_ALL_fire_map.txt



# Hap 1 -------------------------------------------------------------------------------------------------------------------------------

peaks = dict()
peak_regions = dict()
bmap='large_msps_H1_fire_map.txt'
with open(bmap) as fr:
    i = 1
    for line in fr:
        bed,zmw = line.strip().split('|')
        peaks[i] = zmw.split(';')
        bed2 = bed.split('\t')
        peak_regions[i] = f"{bed2[0]}:{bed2[1]}-{bed2[2]}"
        i+=1

# samtools view pileups/GM12878_UBA1_PS00685_nuc_H1.bam | awk '{ print $1 }' > H1_GA_GM12878_UBA1_zmws.txt

# total fiber number
# chrX:47190561-47194939 UBA1
fibers = set()
reg_chrom = 'chrX'
reg_start = 47190561
reg_end = 47194939
bam_name = 'pileups/GM12878_UBA1_PS00685_nuc_H1.bam'
bam = pysam.AlignmentFile(bam_name, "rb")
for read in bam.fetch():
    if read.is_secondary == False and read.is_supplementary == False and read.reference_start <= (reg_start+20) and read.reference_end >= (reg_end-20):
        fibers.add(read.qname)


# DF of zmws by region footprint status (not overlapped by fiber (should be filtered out) -> -1, no MSP -> 0, MSP -> 1)
msp_df = pd.DataFrame(fibers, columns=['zmw'])
for k in peaks.keys():
    new_vals = []
    for i in range(len(msp_df)):
        z = msp_df.iloc[i,:]['zmw']
        if z in peaks[k]:
            new_vals.append(1)
        else:
            new_vals.append(0)
    msp_df[k] = new_vals

msp_df.index = msp_df['zmw']
msp_df = msp_df.drop(['zmw'], axis=1)
msp_df.to_csv('uba1_msp_df_H1.csv')

# calculate codependency for each pair of elements
scores = dict()
rows = []
for combo in combinations(msp_df.columns, 2):
    r1 = combo[0]
    r2 = combo[1]
    counts1 = Counter(msp_df[r1])
    counts2 = Counter(msp_df[r2])
    if ((counts1[0]+counts1[1]) > 0) and ((counts2[0]+counts2[1]) > 0):
        n_both_msp = len(msp_df[(msp_df[r1] == 1) & (msp_df[r2] == 1)])
        n_both_span = len(msp_df[(msp_df[r1].isin([0,1])) & (msp_df[r2].isin([0,1]))])
        pb1 = len(msp_df[(msp_df[r1] == 1) & (msp_df[r2].isin([0,1]))]) / n_both_span
        pb2 = len(msp_df[(msp_df[r2] == 1) & (msp_df[r1].isin([0,1]))]) / n_both_span
        expected = pb1 * pb2
        if n_both_span > 0:
            observed = n_both_msp / n_both_span
            rows.append([r1, r2, pb1, pb2, expected, observed, (observed-expected)*4]) # scale from -1 to 1 !!!
            scores[f'{r1}:{r2}'] = observed-expected

with open('codep_scores_pairs_H1.csv','w') as fw:
    writer = csv.writer(fw)
    writer.writerow(['reg1','reg2','prop1','prop2','expected','observed','score'])
    for r in rows:
        writer.writerow(r)


# Hap 2 -------------------------------------------------------------------------------------------------------------------------------

# bedmap --ec --echo --echo-map-id tss_uba1_masSeq.bed <( sort-bed pileups/GM12878_UBA1_PS00685_H2_large_msp_positions.bed ) > large_msps_H2_fire_map.txt

peaks = dict()
peak_regions = dict()
bmap='large_msps_H2_fire_map.txt'
with open(bmap) as fr:
    i = 1
    for line in fr:
        bed,zmw = line.strip().split('|')
        peaks[i] = zmw.split(';')
        bed2 = bed.split('\t')
        peak_regions[i] = f"{bed2[0]}:{bed2[1]}-{bed2[2]}"
        i+=1


# samtools view pileups/GM12878_UBA1_PS00685_nuc_H1.bam | awk '{ print $1 }' > H1_GA_GM12878_UBA1_zmws.txt

# total fiber number
# chrX:47190561-47194939 UBA1
fibers = set()
reg_chrom = 'chrX'
reg_start = 47190561
reg_end = 47194939
bam_name = 'pileups/GM12878_UBA1_PS00685_nuc_H2.bam'
bam = pysam.AlignmentFile(bam_name, "rb")
for read in bam.fetch():
    if read.is_secondary == False and read.is_supplementary == False and read.reference_start <= (reg_start+20) and read.reference_end >= (reg_end-20):
        fibers.add(read.qname)


# DF of zmws by region footprint status (not overlapped by fiber (should be filtered out) -> -1, no MSP -> 0, MSP -> 1)
msp_df = pd.DataFrame(fibers, columns=['zmw'])
for k in peaks.keys():
    new_vals = []
    for i in range(len(msp_df)):
        z = msp_df.iloc[i,:]['zmw']
        if z in peaks[k]:
            new_vals.append(1)
        else:
            new_vals.append(0)
    msp_df[k] = new_vals

msp_df.index = msp_df['zmw']
msp_df = msp_df.drop(['zmw'], axis=1)
msp_df.to_csv('uba1_msp_df_H2.csv')

# calculate codependency for each pair of elements
scores = dict()
rows = []
for combo in combinations(msp_df.columns, 2):
    r1 = combo[0]
    r2 = combo[1]
    counts1 = Counter(msp_df[r1])
    counts2 = Counter(msp_df[r2])
    if ((counts1[0]+counts1[1]) > 0) and ((counts2[0]+counts2[1]) > 0):
        n_both_msp = len(msp_df[(msp_df[r1] == 1) & (msp_df[r2] == 1)])
        n_both_span = len(msp_df[(msp_df[r1].isin([0,1])) & (msp_df[r2].isin([0,1]))])
        pb1 = len(msp_df[(msp_df[r1] == 1) & (msp_df[r2].isin([0,1]))]) / n_both_span
        pb2 = len(msp_df[(msp_df[r2] == 1) & (msp_df[r1].isin([0,1]))]) / n_both_span
        expected = pb1 * pb2
        if n_both_span > 0:
            observed = n_both_msp / n_both_span
            rows.append([r1, r2, pb1, pb2, expected, observed, (observed-expected)*4]) # scale from -1 to 1 !!!
            scores[f'{r1}:{r2}'] = observed-expected

with open('codep_scores_pairs_H2.csv','w') as fw:
    writer = csv.writer(fw)
    writer.writerow(['reg1','reg2','prop1','prop2','expected','observed','score'])
    for r in rows:
        writer.writerow(r)


# Combined codependency (Merged H1 & H2) -------------------------------------------------------------------------------------------------------------------------------

msp_df_h1 = pd.read_csv('uba1_msp_df_H1.csv', index_col='zmw')
msp_df_h2 = pd.read_csv('uba1_msp_df_H2.csv', index_col='zmw')
both_haps = pd.concat([msp_df_h1,msp_df_h2], axis=0)

# calculate codependency for each pair of elements
scores = dict()
rows = []
for combo in combinations(both_haps.columns, 2):
    r1 = combo[0]
    r2 = combo[1]
    counts1 = Counter(both_haps[r1])
    counts2 = Counter(both_haps[r2])
    if ((counts1[0]+counts1[1]) > 0) and ((counts2[0]+counts2[1]) > 0):
        n_both_msp = len(both_haps[(both_haps[r1] == 1) & (both_haps[r2] == 1)])
        n_both_span = len(both_haps[(both_haps[r1].isin([0,1])) & (both_haps[r2].isin([0,1]))])
        pb1 = len(both_haps[(both_haps[r1] == 1) & (both_haps[r2].isin([0,1]))]) / n_both_span
        pb2 = len(both_haps[(both_haps[r2] == 1) & (both_haps[r1].isin([0,1]))]) / n_both_span
        expected = pb1 * pb2
        if n_both_span > 0:
            observed = n_both_msp / n_both_span
            rows.append([r1, r2, pb1, pb2, expected, observed, (observed-expected)*4]) # scale from -1 to 1 !!!
            scores[f'{r1}:{r2}'] = observed-expected

with open('codep_scores_pairs_merged_bothHaps.csv','w') as fw:
    writer = csv.writer(fw)
    writer.writerow(['reg1','reg2','prop1','prop2','expected','observed','score'])
    for r in rows:
        writer.writerow(r)
