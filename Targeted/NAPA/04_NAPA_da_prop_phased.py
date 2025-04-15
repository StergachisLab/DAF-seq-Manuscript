import pysam
import numpy as np
import pandas as pd
import csv
import os

bam_name = 'NAPA_PS00626_haplotype_corrected.bam'
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

h1_pos_da = {}
h2_pos_da = {}

# Allow for non-ref C/G to be modified. % R/G|A & % Y/C|T
for read in bam.fetch(reg_chrom, reg_start, reg_end):
    if read.is_secondary == False and read.is_supplementary == False and read.has_tag('HP'):
        strand = read.get_tag('ST')
        hap = read.get_tag('HP')
        seq = read.seq
        pair = read.get_aligned_pairs(matches_only=False, with_seq=True)
        for pos in pair:
            mcoord = pos[0]
            ref_coord = pos[1]
            if mcoord != None and pos[1] != None:
                if strand == 'CT':
                    if hap == 1:
                        if ref_coord in hap_basecalls.keys():
                            base = hap_basecalls[ref_coord][0]
                        else:
                            base = pos[2].upper()
                        if base == 'C':                       
                            if ref_coord not in h1_pos_da.keys():
                                h1_pos_da[ref_coord] = {'total':0, 'da':0}
                            h1_pos_da[ref_coord]['total'] += 1
                            if seq[mcoord] == 'Y':
                                h1_pos_da[ref_coord]['da'] += 1
                    elif hap == 2:
                        if ref_coord in hap_basecalls.keys():
                            base = hap_basecalls[ref_coord][1]
                        else:
                            base = pos[2].upper()
                        if base == 'C':                       
                            if ref_coord not in h2_pos_da.keys():
                                h2_pos_da[ref_coord] = {'total':0, 'da':0}
                            h2_pos_da[ref_coord]['total'] += 1
                            if seq[mcoord] == 'Y':
                                h2_pos_da[ref_coord]['da'] += 1
                elif strand == 'GA':
                    if hap == 1:
                        if ref_coord in hap_basecalls.keys():
                            base = hap_basecalls[ref_coord][0]
                        else:
                            base = pos[2].upper()
                        if base == 'G':                       
                            if ref_coord not in h1_pos_da.keys():
                                h1_pos_da[ref_coord] = {'total':0, 'da':0}
                            h1_pos_da[ref_coord]['total'] += 1
                            if seq[mcoord] == 'R':
                                h1_pos_da[ref_coord]['da'] += 1
                    elif hap == 2:
                        if ref_coord in hap_basecalls.keys():
                            base = hap_basecalls[ref_coord][1]
                        else:
                            base = pos[2].upper()
                        if base == 'G':                       
                            if ref_coord not in h2_pos_da.keys():
                                h2_pos_da[ref_coord] = {'total':0, 'da':0}
                            h2_pos_da[ref_coord]['total'] += 1
                            if seq[mcoord] == 'R':
                                h2_pos_da[ref_coord]['da'] += 1


# out bedgraph of DA proportions (both haps)
out_bg = 'NAPA_PS00626_DA_density.bg'
with open(out_bg,'w') as fout:
    writer = csv.writer(fout, delimiter="\t")
    for k in h1_pos_da.keys():
        if k in h2_pos_da.keys():
            total = h1_pos_da[k]['total'] + h2_pos_da[k]['total']
            da = h1_pos_da[k]['da'] + h2_pos_da[k]['da']
        else:
            total = h1_pos_da[k]['total']
            da = h1_pos_da[k]['da']
        if da > 0: # filter 0 positions
            writer.writerow([reg_chrom, k, k+1, da/total])
    for k in h2_pos_da.keys():
        if k not in h1_pos_da.keys():
            total = h2_pos_da[k]['total']
            da = h2_pos_da[k]['da']
            if da > 0: # filter 0 positions
                writer.writerow([reg_chrom, k, k+1, da/total])

# Hap 1 only
out_bg = 'NAPA_PS00626_Hap1_DA_density.bg'
with open(out_bg,'w') as fout:
    writer = csv.writer(fout, delimiter="\t")
    for k,v in h1_pos_da.items():
        if v['da'] > 0: # filter 0 positions
            writer.writerow([reg_chrom, k, k+1, v['da']/v['total']])

# Hap 2 only
out_bg = 'NAPA_PS00626_Hap2_DA_density.bg'
with open(out_bg,'w') as fout:
    writer = csv.writer(fout, delimiter="\t")
    for k,v in h2_pos_da.items():
        if v['da'] > 0: # filter 0 positions
            writer.writerow([reg_chrom, k, k+1, v['da']/v['total']])


# sort -k1,1 -k2,2n NAPA_PS00626_DA_density.bg > temp.bg
# mv temp.bg NAPA_PS00626_DA_density.bg
# /gscratch/stergachislab/install_dir/bedGraphToBigWig NAPA_PS00626_DA_density.bg /mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.chrom.sizes NAPA_PS00626_DA_density.bw

# sort -k1,1 -k2,2n NAPA_PS00626_Hap1_DA_density.bg > temp.bg
# mv temp.bg NAPA_PS00626_Hap1_DA_density.bg
# /gscratch/stergachislab/install_dir/bedGraphToBigWig NAPA_PS00626_Hap1_DA_density.bg /mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.chrom.sizes NAPA_PS00626_Hap1_DA_density.bw

# sort -k1,1 -k2,2n NAPA_PS00626_Hap2_DA_density.bg > temp.bg
# mv temp.bg NAPA_PS00626_Hap2_DA_density.bg
# /gscratch/stergachislab/install_dir/bedGraphToBigWig NAPA_PS00626_Hap2_DA_density.bg /mmfs1/gscratch/stergachislab/assemblies/hg38.analysisSet.chrom.sizes NAPA_PS00626_Hap2_DA_density.bw

