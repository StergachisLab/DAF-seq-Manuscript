""" Aggregate collapsed read stats: read length, coverage, phasing. """

import pandas as pd
import numpy as np
from glob import glob
import pysam


samples = ['PS00718','PS00756','PS00757','PS00758']

out_file = 'collapsed_read_stats.tsv'


def calc_N50(len_list):
    total_bases = sum(len_list)
    tracked_bases = 0
    for l in sorted(len_list, reverse=True):
        tracked_bases += l
        if (tracked_bases/total_bases) >= 0.5:
            return(l)
    return(None)


# raw collapsed read stats ------------------------------------------------------------------
sample_stats = {s:{} for s in samples}
for s in samples:
    fastas = glob(f'/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/collapse/consensus_seqs/{s}_*.fa')
    lengths = []
    for f in fastas:
        with open(f) as fr:
            for line in fr:
                if line.startswith(">") == False:
                    lengths.append(len(line.strip()))
    sample_stats[s]['N50'] = calc_N50(lengths)
    a,b,c,d = np.percentile(lengths, [75,50,25,90])
    sample_stats[s]['Q90'] = d
    sample_stats[s]['Q75'] = a
    sample_stats[s]['Q50'] = b
    sample_stats[s]['Q25'] = c
    sample_stats[s]['mean'] = np.mean(lengths)
    sample_stats[s]['longest'] = max(lengths)
    sample_stats[s]['num_100Kb'] = sum(i >= 100000 for i in lengths)
    sample_stats[s]['ultra_long_bases'] = sum([i for i in lengths if i >= 100000])


# Calculated read counts and seq bases using Samtools

# PS00718_DddA_A2_HighYield_HG002_m84055_240822_201326_s2.hifi_reads.bc2068.bam, PS00718_DddA_A2_Seq2_m84046_241004_204600_s1.hifi_reads.bc2068.bam
# PS00756_DddA_C2_Seq2_m84046_241004_224513_s2.hifi_reads.bc2001.bam
# PS00757_DddA_E2_Seq2_m84046_241005_004435_s3.hifi_reads.bc2002.bam
# PS00758_DddA_F2_Seq2_m84046_241005_024355_s4.hifi_reads.bc2003.bam, PS00758_DddA_F2_Seq2_m84046_241014_163432_s2.hifi_reads.bc2003.bam, PS00758_DddA_F2_Seq2_m84046_241021_233822_s1.hifi_reads.bc2003.bam, PS00758_DddA_F2_Seq2_m84046_241104_232153_s1.hifi_reads.bc2003.bam

# add raw read counts -----------------------------
# total raw reads
# PS00718 6548122
# PS00756 6028101
# PS00757 5732285
# PS00758 6784775 + 8298680 + 5645855 + 7005740 = 27735050

raw_read_counts = {'PS00718':6548122 ,'PS00756':6028101,'PS00757':5732285,'PS00758':27735050}
for k,v in raw_read_counts.items():
    sample_stats[k]['raw_reads'] = v

# total sequenced bases -----------------------------
# samtools stats BAM | grep ^SN | cut -f 2- | grep total &

# IGNORES CLIPPING!
# PS00718 5325652896 + 17495394154
# PS00756 29137666533
# PS00757 22096376656
# PS00758 22217877534 + 25719673247 + 30400003589 + 22938518498

raw_read_bases = {'PS00718':22821047050 ,'PS00756':29137666533,'PS00757':22096376656,'PS00758':101276072868}
for k,v in raw_read_bases.items():
    sample_stats[k]['bases'] = v


# coverage of the MAPPABLE Genome (hg38) -----------------------------------------------------

# Precomputed inuputs with phasing/mappable_coverage_precompute.sh !!!

cov_dir='/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/phasing'
cov_bed = f'{cov_dir}/mappable_coverage_precompute.tsv'
with open(cov_bed) as fr:
    for line in fr:
        sname,group,count = line.strip().split()
        sample_stats[sname][group] = int(count)

# total mappable bases    
with open(f'{cov_dir}/mappable_hg38_precompute.txt') as fr:
    header = next(fr)
    tot_map_bp = int(next(fr).strip())

for s in samples:
    sample_stats[s]['mappable_ref_bp'] = tot_map_bp


# (Same peaks and coverage regardless of MSP length)
fire_df = pd.read_csv('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/msp_analysis/scDAF_FIRE_actuation_MSP150.tsv', sep="\t") # 150 bp MSP

for s in samples:
    # FIRE peaks covered
    vc = fire_df[f'{s}_ALL'].value_counts()
    sample_stats[s]['FIRE_cov_ALL'] = vc[0]+vc[1]
    vc = fire_df[f'{s}_H1'].value_counts()
    sample_stats[s]['FIRE_cov_H1'] = vc[0]+vc[1]
    vc = fire_df[f'{s}_H2'].value_counts()
    sample_stats[s]['FIRE_cov_H2'] = vc[0]+vc[1]
    sample_stats[s]['FIRE_cov_BOTH_haps'] = len(fire_df[(fire_df[f'{s}_H1'] >= 0) & (fire_df[f'{s}_H2'] >= 0)])
    sample_stats[s]['mappable_peaks'] = len(fire_df)

stat_df = pd.DataFrame(sample_stats).T
stat_df.to_csv(out_file, sep="\t", float_format='%.0f', index_label='Cell')
