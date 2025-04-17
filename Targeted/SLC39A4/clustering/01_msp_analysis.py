import pandas as pd
import numpy as np
import pysam
from collections import Counter


bam_dir = '../'
bams = ['GM12878_SLC39A4_PS00686_haplotype_corrected.bam','Liver_SLC39A4_PS00680_haplotype_corrected.bam']

# extract counts for EACH position of fibers
# GA and CT separately

# chr8:144,415,793-144,417,939 SLC39A4
reg_chrom = 'chr8'
reg_start = 144415793
reg_end = 144417939

n_reads = 10000
prop_cutoff = 0.9
ct_pos_by_t = []
ga_pos_by_t = []
for b in bams:
    bam = pysam.AlignmentFile(f'{bam_dir}/{b}', "rb")
    ct_pos = dict() # basecalls by CT & GA strands
    ga_pos = dict() # basecalls by CT & GA strands
    n_ct = 0
    n_ga = 0
    for read in bam.fetch(reg_chrom, reg_start, reg_end):
        if read.is_secondary == False and read.is_supplementary == False:
            seq = read.seq
            pair = read.get_aligned_pairs(matches_only=False, with_seq=True)
            strand = read.get_tag('ST')
            if strand == 'CT' and 'R' not in seq and n_ct < n_reads:
                for pos in pair:
                    if pos[0] != None and pos[1] != None: # ignore insertions & deletions
                        mcoord = pos[0]
                        ref_coord = pos[1] + 1
                        if ref_coord not in ct_pos.keys():
                            ct_pos[ref_coord] = {'A':0, 'C':0, 'G':0, 'T':0, 'Y':0}
                        ct_pos[ref_coord][seq[mcoord]] += 1
                n_ct += 1
            elif strand == 'GA' and 'Y' not in seq and n_ga < n_reads:
                for pos in pair:
                    if pos[0] != None and pos[1] != None: # ignore insertions & deletions
                        mcoord = pos[0]
                        ref_coord = pos[1] + 1
                        if ref_coord not in ga_pos.keys():
                            ga_pos[ref_coord] = {'A':0, 'C':0, 'G':0, 'T':0, 'R':0}
                        ga_pos[ref_coord][seq[mcoord]] += 1
                n_ga += 1
        if n_ct >= n_reads and n_ga >= n_reads:
            break
    ct_set = set()
    for k,v in ct_pos.items():
        if ((v['C']+v['Y']) / sum(v.values())) > prop_cutoff and (sum(v.values())/n_reads) > prop_cutoff: # ensure high counts for position, remove edges
            ct_set.add(k)
    ct_pos_by_t.append(ct_set)
    ga_set = set()
    for k,v in ga_pos.items():
        if ((v['G']+v['R']) / sum(v.values())) > prop_cutoff and (sum(v.values())/n_reads) > prop_cutoff: # ensure high counts for position, remove edges
            ga_set.add(k)
    ga_pos_by_t.append(ga_set)

shared_CT_keys = ct_pos_by_t[0].intersection(*ct_pos_by_t)
shared_GA_keys = ga_pos_by_t[0].intersection(*ga_pos_by_t)

with open('shared_keys_CT.txt','w') as fw:
    for k in sorted(shared_CT_keys):
        fw.write(f'{k}\n')
with open('shared_keys_GA.txt','w') as fw:
    for k in sorted(shared_GA_keys):
        fw.write(f'{k}\n')

# read postions
shared_CT_keys = set()
with open('shared_keys_CT.txt') as fr:
    for line in fr:
        shared_CT_keys.add(int(line.strip()))
shared_GA_keys = set()
with open('shared_keys_GA.txt') as fr:
    for line in fr:
        shared_GA_keys.add(int(line.strip()))


ct_dfs = []
ga_dfs = []
for b in bams:
    tissue = b.split('_')[0]
    bam = pysam.AlignmentFile(f'{bam_dir}/{b}', "rb")
    ct_fibers = dict() # basecalls by CT & GA strands
    ga_fibers = dict() # basecalls by CT & GA strands
    for read in bam.fetch(reg_chrom, reg_start, reg_end):
        if read.is_secondary == False and read.is_supplementary == False:
            zmw = read.qname
            seq = read.seq
            pair = read.get_aligned_pairs(matches_only=False, with_seq=True)
            strand = read.get_tag('ST')
            if read.has_tag('HP'):
                hap = read.get_tag('HP')
                if strand == 'CT' and 'R' not in seq:
                    ct_fibers[zmw] = {p:-1 for p in sorted(shared_CT_keys)}
                    ct_fibers[zmw]['hap'] = hap
                    for pos in pair:
                        if pos[0] != None and pos[1] != None: # ignore insertions & deletions
                            mcoord = pos[0]
                            ref_coord = pos[1] + 1
                            if ref_coord in shared_CT_keys:
                                if seq[mcoord] == 'Y':
                                    ct_fibers[zmw][ref_coord] = 1
                                elif seq[mcoord] == 'C': # NEW
                                    ct_fibers[zmw][ref_coord] = 0 # NEW
                elif strand == 'GA' and 'Y' not in seq:
                    ga_fibers[zmw] = {p:0 for p in sorted(shared_GA_keys)}
                    ga_fibers[zmw]['hap'] = hap
                    for pos in pair:
                        if pos[0] != None and pos[1] != None: # ignore insertions & deletions
                            mcoord = pos[0]
                            ref_coord = pos[1] + 1
                            if ref_coord in shared_GA_keys:
                                if seq[mcoord] == 'R':
                                    ga_fibers[zmw][ref_coord] = 1
                                if seq[mcoord] == 'G': # NEW
                                    ga_fibers[zmw][ref_coord] = 0 # NEW
    c_df = pd.DataFrame(ct_fibers).T
    c_df['tissue'] = tissue
    g_df = pd.DataFrame(ga_fibers).T
    g_df['tissue'] = tissue
    ct_dfs.append(c_df)
    ga_dfs.append(g_df)
    print(f'DONE for {b}')

all_ct = pd.concat(ct_dfs, axis=0)
all_ga = pd.concat(ga_dfs, axis=0)
all_ct['zmw'] = all_ct.index
all_ga['zmw'] = all_ga.index
all_ct.to_csv('tissue_DA_features_CT.csv', index=False)
all_ga.to_csv('tissue_DA_features_GA.csv', index=False)
