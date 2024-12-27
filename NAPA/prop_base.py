""" At each base, plot the proportion of each base on top vs bottom strands. """

import pysam
import csv
import pandas as pd
import numpy as np

bam_name = '/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/PCR_Dev/Napa_WASF1/PS00626.m84046_240619_124816_s1.bc2072.ft.map-pb_corrected_realigned.bam'
bam = pysam.AlignmentFile(bam_name, "rb")

# hap basecalls
basecall_tsv = 'NAPA_PS00626_haplotype_corrected_hap_basecalls.tsv'
hap_basecalls = dict()
with open(basecall_tsv) as fr:
    reader = csv.reader(fr, delimiter="\t")
    for line in reader:
        if line[0] != 'Position':
            hap_basecalls[int(line[0])-1] = [line[1], line[2]]


# chr19:47,514,458-47,519,061 NAPA
reg_chrom = 'chr19'
reg_start = 47514458
reg_end = 47519061

# basecalls by CT & GA strands
ct_pos = dict()
ga_pos = dict()

# ref coordinates are 1-based !!!!!!!!!!
for read in bam.fetch(reg_chrom, reg_start, reg_end):
    if read.is_secondary == False and read.is_supplementary == False:
        seq = read.seq
        pair = read.get_aligned_pairs(matches_only=False, with_seq=True)
        strand = read.get_tag('ST')
        if strand == 'CT':
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
        elif strand == 'GA':
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

# convert to DataFrames
# CT is top strand, GA is bottom strand
bases = ['A','C','G','T']
for b in bases:
    keys = []
    top_prop = []
    bottom_prop = []
    for key in shared_keys:
        keys.append(key)
        top_prop.append(ct_pos[key][b]/sum(ct_pos[key].values()))
        bottom_prop.append(ga_pos[key][b]/sum(ga_pos[key].values()))
    df = pd.DataFrame({'top_prop':top_prop, 'bottom_prop':bottom_prop}, index = keys)
    df.to_csv(f'prop_{b}_by_strand.csv', index_label='position')


# # reload DFs
# df_a = pd.read_csv('prop_A_by_strand.csv', index_col='position')
# df_c = pd.read_csv('prop_C_by_strand.csv', index_col='position')
# df_g = pd.read_csv('prop_G_by_strand.csv', index_col='position')
# df_t = pd.read_csv('prop_T_by_strand.csv', index_col='position')


# # plot prop top strand on Y and prop bottom strand on X for each base

# import altair as alt
# chart = alt.Chart(df_c)

# alt.Chart(df_c).mark_point()
# chart.save('chart.png')

# alt.renderers.enable("html")

# alt.renderers.enable('png')
# alt.renderers.enable('image/png')
