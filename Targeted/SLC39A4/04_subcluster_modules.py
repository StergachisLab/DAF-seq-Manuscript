import csv
import os
from glob import glob
import pandas as pd
from collections import Counter
import pysam
from itertools import combinations



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

