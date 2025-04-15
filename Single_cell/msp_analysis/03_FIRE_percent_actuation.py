import pandas as pd
import numpy as np
import pysam


# track single-cell MSP overlaps as 1, non-MSP overlap with coverage 0, no coverage -1

# initialize peak DF
peaks='FIRE-peaks_HG002_noUnreliableCoverage.AutosomesOnly.bed.gz'
peak_df = pd.read_csv(peaks, sep="\t", names=['chrom','start','end','name','prop_acc','strand'])
daf_df = peak_df[['chrom','start','end','prop_acc']]
daf_df.index = daf_df["chrom"].astype(str) + ':' + daf_df["start"].astype(str) + ':' + daf_df["end"].astype(str)

# Bin FIRE peaks by 10% actuation
daf_df['bin'] = daf_df['prop_acc'] // .1

samples=['PS00718','PS00756','PS00757','PS00758','PS00867','PS00868','PS00869','PS00870','PS00871','PS00872','PS00873','PS00874']
msp_lengths = [150]

for LEN in msp_lengths:
    # read in CT and GA BED overlaps as separate files
    for PS in samples:
        fiber_phase = dict()
        # Track fiber Hap info
        bam_file = f"fibertools_msp/{PS}_consensus_BothStrands_HG38_corrected.haplotagged.m6A_nuc.bam"
        bam = pysam.AlignmentFile(bam_file, "rb")
        for read in bam.fetch():
            if read.is_secondary == False and read.is_supplementary == False:
                try:
                    hap = read.get_tag("HP")
                except:
                    hap = "UNK"
                fiber_phase[read.qname] = hap
        # Fiber overlaps
        fiber_int_file=f'fibertools_msp/{PS}_collapsed_haplotagged_FIRE_read_intersect.bed.gz'
        fiber_int = pd.read_csv(fiber_int_file, sep="\t", names=['chrom','start','end','name','score','strand','rgb','v1','v2','v3','v4','v5'])
        fiber_int = fiber_int[['chrom','start','end','name']]
        fiber_int['name'] = fiber_int['name'].astype(str)
        fiber_int.index = fiber_int["chrom"].astype(str) + ':' + fiber_int["start"].astype(str) + ':' + fiber_int["end"].astype(str)
        # MSPs
        msp_int_file=f'fibertools_msp/{PS}_collapsed_haplotagged_FIRE_MSP_{LEN}_intersect.bed.gz'
        msp_int = pd.read_csv(msp_int_file, sep="\t", names=['chrom','start','end','name','score','strand'])
        msp_int = msp_int[['chrom','start','end','name']]
        msp_int['name'] = msp_int['name'].astype(str)
        msp_int.index = msp_int["chrom"].astype(str) + ':' + msp_int["start"].astype(str) + ':' + msp_int["end"].astype(str)
        # fill out FIRE peak DF for sample
        calls_all = [-1 for i in daf_df.index]
        fiber_index_check = set(fiber_int.index)
        for i in range(len(daf_df.index)):
            if daf_df.index[i] in fiber_index_check:
                calls_all[i] = 0
        msp_index_check = set(msp_int.index)
        for i in range(len(daf_df.index)):
            if daf_df.index[i] in msp_index_check:
                calls_all[i] = 1
        daf_df[f'{PS}_ALL'] = calls_all
        # by Hap
        calls_H1 = [-1 for i in daf_df.index]
        calls_H2 = [-1 for i in daf_df.index]
        daf_index_check = list(daf_df.index)
        fiber_int_idx = list(fiber_int.index)
        for i in range(len(fiber_int)):
            fiber_name = fiber_int.iloc[i]['name']
            f_peak = fiber_int_idx[i]
            if fiber_name in fiber_phase.keys():
                if fiber_phase[fiber_name] == 1:
                    calls_H1[daf_index_check.index(f_peak)] = 0
                elif fiber_phase[fiber_name] == 2:
                    calls_H2[daf_index_check.index(f_peak)] = 0
        msp_int_idx = list(msp_int.index)
        for i in range(len(msp_int)):
            fiber_name = msp_int.iloc[i]['name']
            f_peak = msp_int_idx[i]
            if fiber_name in fiber_phase.keys():
                if fiber_phase[fiber_name] == 1:
                    calls_H1[daf_index_check.index(f_peak)] = 1
                elif fiber_phase[fiber_name] == 2:
                    calls_H2[daf_index_check.index(f_peak)] = 1
        daf_df[f'{PS}_H1'] = calls_H1
        daf_df[f'{PS}_H2'] = calls_H2
        # prop accessible of haps covered
        daf_df[f'{PS}_nCov'] = daf_df[[f'{PS}_H1',f'{PS}_H2']].apply(lambda row: (row >= 0).sum(), axis=1)
        daf_df[f'{PS}_nMSP'] = daf_df[[f'{PS}_H1',f'{PS}_H2']].apply(lambda row: (row >= 1).sum(), axis=1)
    daf_df.to_csv(f'scDAF_FIRE_actuation_MSP{LEN}.tsv', sep="\t", index=False)

