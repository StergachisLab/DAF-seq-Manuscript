""" Identify heterozygous variants from DAF-seq data based on the % base at each position in the target region.
"""

import pysam
import csv
import pandas as pd
import numpy as np


# At each base, plot the proportion of each base on top vs bottom strands.


# length from the end of the region to consider (priming regions)
buffer_len = 25

samples = [{'name':'NAPA', 'reg_chrom':'chr19', 'reg_start':47514488, 'reg_end':47518985, 'bam':'../data/PS00626.m84046_240619_124816_s1.bc2072.ft.map-pb_corrected_realigned.bam'},
           {'name':'WASF1', 'reg_chrom':'chr6', 'reg_start':110176895, 'reg_end':110181442, 'bam':'../data/PS00626.m84046_240619_124816_s1.bc2072.ft.map-pb_corrected_realigned.bam'},
           {'name':'UBA1', 'reg_chrom':'chrX', 'reg_start':47190561, 'reg_end':47194939, 'bam':'../data/PS00685_GM12878_UBA1_m84046_240802_231812_s1.hifi_reads.bc2093.map-pb_corrected_realigned.bam'},
           {'name':'GM12878_SLC', 'reg_chrom':'chr8', 'reg_start':144415793, 'reg_end':144417939, 'bam':'../data/GM12878_SLC39A4_PS00686_m84046_240802_231812_s1.hifi_reads.bc2094.map-pb_corrected_realigned.bam'},
           {'name':'Liver_SLC', 'reg_chrom':'chr8', 'reg_start':144415793, 'reg_end':144417939, 'bam':'../data/Liver_SLC39A4_PS00680_m84046_240802_231812_s1.hifi_reads.bc2084.map-pb_corrected_realigned.bam'},
           {'name':'Heart_SLC', 'reg_chrom':'chr8', 'reg_start':144415793, 'reg_end':144417939, 'bam':'../data/Heart_SLC39A4_PS00681_m84046_240802_231812_s1.hifi_reads.bc2085.map-pb_corrected_realigned.bam'},
           {'name':'Colon_SLC', 'reg_chrom':'chr8', 'reg_start':144415793, 'reg_end':144417939, 'bam':'../data/Colon_SLC39A4_PS00682_m84046_240802_231812_s1.hifi_reads.bc2086.map-pb_corrected_realigned.bam'},
           {'name':'COLO_Mix', 'reg_chrom':'chr17', 'reg_start':19446516, 'reg_end':19450322, 'bam':'../data/PS00719_COLO_Region2_m84055_240822_221245_s3.hifi_reads.bc2069.map-pb_corrected_realigned.bam'}]

for s in samples:
    print(f"Getting strand % DA on: {s['name']}.....")
    bam = pysam.AlignmentFile(s['bam'], "rb")
    ct_pos = dict() # basecalls by CT & GA strands
    ga_pos = dict() # basecalls by CT & GA strands
    for read in bam.fetch(s['reg_chrom'], s['reg_start'], s['reg_end']):
        if read.is_secondary == False and read.is_supplementary == False:
            seq = read.seq
            pair = read.get_aligned_pairs(matches_only=False, with_seq=True)
            strand = read.get_tag('ST')
            if strand == 'CT' and 'R' not in seq:
                for pos in pair:
                    if pos[0] != None and pos[1] != None: # ignore insertions & deletions
                        mcoord = pos[0]
                        ref_coord = pos[1] + 1
                        if ref_coord not in ct_pos.keys():
                            ct_pos[ref_coord] = {'A':0, 'C':0, 'G':0, 'T':0}
                        if seq[mcoord] == 'Y':
                            ct_pos[ref_coord]['T'] += 1
                        else:
                            ct_pos[ref_coord][seq[mcoord]] += 1
            elif strand == 'GA' and 'Y' not in seq:
                for pos in pair:
                    if pos[0] != None and pos[1] != None: # ignore insertions & deletions
                        mcoord = pos[0]
                        ref_coord = pos[1] + 1
                        if ref_coord not in ga_pos.keys():
                            ga_pos[ref_coord] = {'A':0, 'C':0, 'G':0, 'T':0}
                        if seq[mcoord] == 'R':
                            ga_pos[ref_coord]['A'] += 1
                        else:
                            ga_pos[ref_coord][seq[mcoord]] += 1
    shared_keys = set(ct_pos.keys()).intersection(set(ga_pos.keys()))
    bases = ['A','C','G','T']
    for b in bases:
        keys = []
        top_prop = []
        bottom_prop = []
        for key in shared_keys:
            if key >= s['reg_start']+buffer_len and key <= s['reg_end']-buffer_len:
                keys.append(key)
                top_prop.append(ct_pos[key][b]/sum(ct_pos[key].values())) # CT is top strand
                bottom_prop.append(ga_pos[key][b]/sum(ga_pos[key].values())) # GA is bottom strand
        df = pd.DataFrame({'top_prop':top_prop, 'bottom_prop':bottom_prop}, index = keys)
        df.to_csv(f"prop_by_strand/{s['name']}_prop_{b}_by_strand.csv", index_label='position')


