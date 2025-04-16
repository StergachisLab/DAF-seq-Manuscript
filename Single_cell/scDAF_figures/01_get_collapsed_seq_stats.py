""" Aggregate collapsed read stats: read length, coverage, phasing. """

import pandas as pd
import numpy as np
from glob import glob
import csv
import os
import pysam


samples=['PS00718','PS00756','PS00757','PS00758','PS00867','PS00868','PS00869','PS00870','PS00871','PS00872','PS00873','PS00874']

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


# Ran "python /mmfs1/gscratch/stergachislab/swansoe/projects/DddA/seq_stats.py $(ls ./*.bam) -t 30 > read_stats.tbl" for raw seq stats

read_stats = pd.read_csv("/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/data/revision_raw_bams/read_stats.tbl", sep="\t")

read_stats['sample'] = read_stats['file'].str.replace('^./','', regex=True).str.split('_').str[0]
read_stats = read_stats[read_stats['sample'].isin(samples)]

raw_read_counts = read_stats.groupby('sample')['nSeqs'].agg(['sum'])
for i in range(len(raw_read_counts)):
    sample_stats[raw_read_counts.index[i]]['raw_reads'] = raw_read_counts['sum'][i]

raw_read_bases = read_stats.groupby('sample')['totalBp'].agg(['sum'])
for i in range(len(raw_read_bases)):
    sample_stats[raw_read_bases.index[i]]['bases'] = raw_read_bases['sum'][i]


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


# average phasing rate from the four cells
haplotagged_dir = '/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/phasing/haplotagged'
tag_lists = glob(f'{haplotagged_dir}/*HG38_corrected.haplotag_list.txt.gz')

total = 0
phased = 0
for tl in tag_lists:
    df = pd.read_csv(tl, sep="\t")
    df_c = df['haplotype'].value_counts()
    total += sum(df_c)
    phased += df_c['H1'] + df_c['H2']

with open('avg_phased_read_prop.txt','w') as fout:
    fout.write(f'Avg Collapsed Read Phasing Proportion (all cells): {phased/total}\n')


# FIRE stats
fire_df = pd.read_csv('/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/msp_analysis/scDAF_FIRE_actuation_MSP150.tsv', sep="\t")
freq_df = fire_df[['chrom','start','end','prop_acc','bin']]

def check_columns(row):
    if row['H1'] == 0 and row['H2'] == 0:
        return('none')
    elif row['H1'] == 1 and row['H2'] == 0:
        return('H1')
    elif row['H1'] == 0 and row['H2'] == 1:
        return('H2')
    else:
        return('both')

for samp in samples:
    temp_df = fire_df[['chrom','start','end','prop_acc','bin',f'{samp}_H1',f'{samp}_H2']]
    temp_df = temp_df.rename(columns={f'{samp}_H1':'H1', f'{samp}_H2':'H2'})
    filt = temp_df[(temp_df['H1']>=0) & (temp_df['H2']>=0)]
    filt['result'] = filt.apply(check_columns, axis=1)
    freq_df[f'{samp}'] = filt['result']

freq_df['neither'] = freq_df.apply(lambda row: (row == 'none').sum(), axis=1)
freq_df['H1_only'] = freq_df.apply(lambda row: (row == 'H1').sum(), axis=1)
freq_df['H2_only'] = freq_df.apply(lambda row: (row == 'H2').sum(), axis=1)
freq_df['both'] = freq_df.apply(lambda row: (row == 'both').sum(), axis=1)

freq_df.to_csv('hap_actuation_bins.csv', index=False)


# Consensus read lengths
seq_dir = '/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/collapse/consensus_seqs'
fastas = glob(f'{seq_dir}/PS*.fa')

out_rows = []
for file in fastas:
    name = os.path.basename(file).split('_')[0]
    with open(file) as fr:
        for line in fr:
            if not line.startswith('>'):
                out_rows.append((len(line.strip()), name))

with open('consensus_read_lengths.tsv','w') as fw:
    writer = csv.writer(fw, delimiter="\t")
    for row in out_rows:
        writer.writerow(row)


# Deamination rate by cell
bam_dir='/mmfs1/gscratch/stergachislab/swansoe/projects/DddA/single_cell/expt1_HG002_FACS/collapse/consensus_bams'
BAMS = [f'{bam_dir}/{s}_consensus_BothStrands_HG38_corrected.bam' for s in samples]

def da_rate(bam_name):
    da_bp = 0
    total_C = 0
    total_Y = 0
    total_G = 0
    total_R = 0
    bam = pysam.AlignmentFile(bam_name, "rb")
    for read in bam.fetch():
        strand = read.get_tag('ST')
        if strand == 'CT':
            da_bp += read.seq.count('Y')
            total_C += read.seq.count('C')
            total_Y += read.seq.count('Y')
        elif strand == 'GA':
            da_bp += read.seq.count('R')
            total_G += read.seq.count('G')
            total_R += read.seq.count('R')
            pass
    return(da_bp, total_C, total_Y, total_G, total_R)


out_rows = [['Cell','DA_bp','Tot_bp','prop_CT_da','prop_GA_da','prop_da_both_strands']]

for i in range(len(samples)):
    da_bp,total_C,total_Y,total_G,total_R = da_rate(BAMS[i])
    out_rows.append([samples[i], da_bp, da_bp+total_C+total_G, total_Y/(total_Y+total_C), total_R/(total_R+total_G), da_bp/(total_Y+total_R+total_C+total_G)])

with open('deamination_rate_by_cell.csv','w') as fw:
    writer = csv.writer(fw)
    for line in out_rows:
        writer.writerow(line)

