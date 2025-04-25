import pandas as pd
import numpy as np
import pysam
from glob import glob
import os


""" Create dataframe of per-cell and per-haplotype CTCF occupancy """


# Get single-cell read / haplotype-strand associations
fiber_dict = dict()

bams = glob("../msp_analysis/fibertools_msp/*_consensus_BothStrands_HG38_corrected.haplotagged.m6A_nuc.bam")
for b in bams:
    bam = pysam.AlignmentFile(b, "rb")
    sname = os.path.basename(b).split('_')[0]
    for read in bam.fetch():
        if read.is_secondary == False and read.is_supplementary == False:
            rname = read.qname
            strand = read.get_tag('ST')
            if read.has_tag('HP'):
                hap = read.get_tag('HP')
            else:
                hap = 'UNK'
            if strand == 'CT':
                fiber_dict[rname] = {'hap':hap, 'strand':'Top'}
            elif strand == 'GA':
                fiber_dict[rname] = {'cell':sname, 'hap':hap, 'strand':'Bottom'}


df_rows = []
ft_files = glob("ft_out/ft_CTCF_PS*.bed")
for file in ft_files:
    sname = file.rstrip('.bed').split('_')[-1]
    ft_out = pd.read_csv(file, sep="\t")
    ft_out['motif_ID'] = ft_out.index + 1
    ft_dict = {motif:{'1_Top':-1, '2_Top':-1, '1_Bottom':-1, '2_Bottom':-1, 'UNK_Top':-1, 'UNK_Bottom':-1} for motif in ft_out['motif_ID']} # -1 for no overlap, 0 for overlap only, 1 for MSP, 2 for footprinted CTCF
    for row in ft_out.iterrows():
        if row[1].isna().any(): # no overlapping fibers
            continue
        else:
            fiber_names = row[1]['fiber_names'].split(',')
            footprint_codes = [int(fc) for fc in row[1]['footprint_codes'].split(',')]
            for i in range(len(fiber_names)):
                fname = fiber_names[i]
                if fname in fiber_dict.keys():
                    hap = fiber_dict[fname]['hap']
                    strand = fiber_dict[fname]['strand']
                    if (footprint_codes[i] & 1) > 0:
                        msp = True
                    else:
                        msp = False
                    if (footprint_codes[i] & (1 << 3)) > 0 and (footprint_codes[i] & (1 << 4)) > 0 and (footprint_codes[i] & 1) > 0: # mods 2&3 bound and spanning MSP
                        bound = True
                    else:
                        bound = False
                    update_key = f'{hap}_{strand}'
                    if not msp:
                        ft_dict[row[1]['motif_ID']][update_key] = 0
                    else:
                        if not bound:
                            ft_dict[row[1]['motif_ID']][update_key] = 1
                        else:
                            ft_dict[row[1]['motif_ID']][update_key] = 2
    for motif,values in ft_dict.items():
        chrom = ft_out.loc[ft_out['motif_ID'] == motif, '#chrom'].tolist()[0]
        mstart = ft_out.loc[ft_out['motif_ID'] == motif, 'start'].tolist()[0]
        mend = ft_out.loc[ft_out['motif_ID'] == motif, 'end'].tolist()[0]
        mstrand = ft_out.loc[ft_out['motif_ID'] == motif, 'strand'].tolist()[0]
        for hs,code in values.items():
            hap,strand = hs.split('_')
            if hap in ['1','2']:
                hap = f'H{hap}'
            df_rows.append([motif, chrom, mstart, mend, mstrand, sname, hap, strand, code])


out_df = pd.DataFrame(df_rows, columns = ['motif_ID','chrom','start','end','motif_strand','cell','hap','strand','code'])
out_df.to_csv('ctcf_footprint_codes.csv', index=False)


# quantify Fiber-seq footprints for HG002 and GM12878
HG_ft = pd.read_csv("ft_out/ft_CTCF_FS_HG002.bed", sep="\t")
HG_ft['motif_ID'] = HG_ft.index + 1
HG_ft['sample'] = "HG002"

n_overlap_l = []
n_msp_l = []
n_bound_l = []

for row in HG_ft.iterrows():
    n_overlap = 0
    n_msp = 0
    n_bound = 0
    if not row[1].isna().any():
        fiber_names = row[1]['fiber_names'].split(',')
        footprint_codes = [int(fc) for fc in row[1]['footprint_codes'].split(',')]
        for i in range(len(fiber_names)):
            n_overlap += 1
            fname = fiber_names[i]
            if (footprint_codes[i] & 1) > 0:
                n_msp += 1
            if (footprint_codes[i] & (1 << 3)) > 0 and (footprint_codes[i] & (1 << 4)) > 0 and (footprint_codes[i] & 1) > 0: # mods 2&3 bound and spanning MSP
                n_bound += 1
    n_overlap_l.append(n_overlap)
    n_msp_l.append(n_msp)
    n_bound_l.append(n_bound)
        
HG_ft['n_overlap'] = n_overlap_l
HG_ft['n_msp'] = n_msp_l
HG_ft['n_bound'] = n_bound_l
HG_ft['prop_bound'] = HG_ft['n_bound'] / HG_ft['n_msp']
HG_ft['prop_act'] = HG_ft['n_msp'] / HG_ft['n_overlap']


GM_ft = pd.read_csv("ft_out/ft_CTCF_FS_GM12878.bed", sep="\t")
GM_ft['motif_ID'] = GM_ft.index + 1
GM_ft['sample'] = "GM12878"

n_overlap_l = []
n_msp_l = []
n_bound_l = []

for row in GM_ft.iterrows():
    n_overlap = 0
    n_msp = 0
    n_bound = 0
    if not row[1].isna().any():
        fiber_names = row[1]['fiber_names'].split(',')
        footprint_codes = [int(fc) for fc in row[1]['footprint_codes'].split(',')]
        for i in range(len(fiber_names)):
            n_overlap += 1
            fname = fiber_names[i]
            if (footprint_codes[i] & 1) > 0:
                n_msp += 1
            if (footprint_codes[i] & (1 << 3)) > 0 and (footprint_codes[i] & (1 << 4)) > 0 and (footprint_codes[i] & 1) > 0: # mods 2&3 bound and spanning MSP
                n_bound += 1
    n_overlap_l.append(n_overlap)
    n_msp_l.append(n_msp)
    n_bound_l.append(n_bound)
        
GM_ft['n_overlap'] = n_overlap_l
GM_ft['n_msp'] = n_msp_l
GM_ft['n_bound'] = n_bound_l
GM_ft['prop_bound'] = GM_ft['n_bound'] / GM_ft['n_msp']
GM_ft['prop_act'] = GM_ft['n_msp'] / GM_ft['n_overlap']


fs_merged = pd.concat([HG_ft, GM_ft], ignore_index = True)
fs_merged = fs_merged[['motif_ID','#chrom','start','end','strand','sample','n_overlap','n_msp','n_bound','prop_bound','prop_act']]
fs_merged = fs_merged.rename(columns={'#chrom':'chrom'})

fs_merged.to_csv('FiberSeq_merged_ctcf_footprint.csv', index=False)
